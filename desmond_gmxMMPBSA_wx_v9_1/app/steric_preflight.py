#!/usr/bin/env python3
"""Topology-aware steric preflight for prepared gmx_MMPBSA trajectories.

The program uses the actual Amber COM.prmtop Lennard-Jones parameters and
exclusion graph to find severe nonbonded overlaps in every prepared frame.
Only synthetic ACE/NME atoms may be moved automatically.  Source protein and
ligand coordinates are never changed by this tool.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import mdtraj as md
import numpy as np
from scipy.spatial import cKDTree

try:
    import parmed as pmd
except ImportError as exc:  # pragma: no cover - dependency check
    raise SystemExit(
        "ParmEd is required for steric preflight. Activate the gmx_MMPBSA environment."
    ) from exc

VERSION = "1.1.0"


@dataclass(frozen=True)
class PairParameter:
    i: int
    j: int
    interaction_class: str
    scnb: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a prepared complex with the actual Amber topology for "
            "Lennard-Jones overlaps before gmx_MMPBSA."
        )
    )
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--pdb", default=None)
    parser.add_argument("--xtc", default=None)
    parser.add_argument("--prmtop", default="COM.prmtop")
    parser.add_argument("--repair-script", default=None)
    parser.add_argument("--cutoff-nm", type=float, default=0.25)
    parser.add_argument("--pair-limit-kcal", type=float, default=1.0e4)
    parser.add_argument("--frame-limit-kcal", type=float, default=1.0e5)
    parser.add_argument("--hard-distance-nm", type=float, default=0.055)
    parser.add_argument("--top-pairs", type=int, default=100)
    parser.add_argument(
        "--no-repair", action="store_true", help="Diagnose only; never rotate synthetic caps."
    )
    return parser.parse_args()


def pair_key(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def load_report(workdir: Path) -> dict[str, object]:
    path = workdir / "preparation_report.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    report = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    if not isinstance(report, dict):
        raise ValueError("preparation_report.json must contain an object")
    return report


def resolve_input_files(
    workdir: Path, report: dict[str, object], pdb_arg: str | None, xtc_arg: str | None
) -> tuple[Path, Path]:
    files = report.get("files", {})
    if not isinstance(files, dict):
        files = {}
    pdb_name = pdb_arg or str(files.get("complex_pdb", "complex_amber_order.pdb"))
    xtc_name = xtc_arg or str(files.get("complex_xtc", "complex_stride10.xtc"))
    pdb = (workdir / pdb_name).resolve() if not Path(pdb_name).is_absolute() else Path(pdb_name)
    xtc = (workdir / xtc_name).resolve() if not Path(xtc_name).is_absolute() else Path(xtc_name)
    if not pdb.is_file():
        raise FileNotFoundError(pdb)
    if not xtc.is_file():
        raise FileNotFoundError(xtc)
    return pdb, xtc


def synthetic_residue_numbers(workdir: Path) -> set[int]:
    path = workdir / "protein_templates.tsv"
    numbers: set[int] = set()
    if not path.is_file():
        return numbers
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("source_type", "")).strip() != "synthetic_cap":
                continue
            try:
                numbers.add(int(str(row.get("output_residue", "")).strip()))
            except ValueError:
                continue
    return numbers


def topology_sets(structure: object) -> tuple[set[tuple[int, int]], dict[tuple[int, int], float]]:
    """Return excluded 1-2/1-3 pairs and scaled 1-4 pairs."""

    excluded: set[tuple[int, int]] = set()
    one_four: dict[tuple[int, int], float] = {}
    for bond in structure.bonds:
        excluded.add(pair_key(int(bond.atom1.idx), int(bond.atom2.idx)))
    for angle in structure.angles:
        excluded.add(pair_key(int(angle.atom1.idx), int(angle.atom3.idx)))
    # ParmEd exposes the Amber exclusion list on atoms.  This captures any
    # extra exclusions that are not reducible to the simple bond/angle graph
    # (for example, special points in a converted topology).  Amber 1-4 pairs
    # also occur in that list, so they are restored below as scaled pairs.
    for atom in structure.atoms:
        for partner in getattr(atom, "exclusion_partners", ()):
            excluded.add(pair_key(int(atom.idx), int(partner.idx)))
    for dihedral in structure.dihedrals:
        if bool(getattr(dihedral, "improper", False)):
            continue
        key = pair_key(int(dihedral.atom1.idx), int(dihedral.atom4.idx))
        dtype = getattr(dihedral, "type", None)
        types: Iterable[object]
        if dtype is None:
            types = ()
        elif isinstance(dtype, (list, tuple)):
            types = dtype
        else:
            try:
                # DihedralTypeList is iterable, ordinary DihedralType is not.
                types = list(dtype) if dtype.__class__.__name__.endswith("TypeList") else (dtype,)
            except TypeError:
                types = (dtype,)
        scnb_values = [
            float(getattr(item, "scnb"))
            for item in types
            if getattr(item, "scnb", None) not in (None, 0)
        ]
        scnb = scnb_values[0] if scnb_values else 2.0
        one_four[key] = scnb
    # 1-4 pairs are absent from the ordinary nonbonded list but are evaluated
    # separately with SCNB scaling; do not discard them as generic exclusions.
    excluded.difference_update(one_four)
    return excluded, one_four


def atom_lj(atom: object) -> tuple[float, float]:
    epsilon = float(getattr(atom, "epsilon", 0.0) or 0.0)
    rmin_half = float(getattr(atom, "rmin", 0.0) or 0.0)
    return abs(epsilon), max(0.0, rmin_half)


def lj_energy_kcal(
    atom_i: object,
    atom_j: object,
    distance_nm: float,
    scnb: float = 1.0,
) -> float:
    epsilon_i, rmin_i = atom_lj(atom_i)
    epsilon_j, rmin_j = atom_lj(atom_j)
    if epsilon_i <= 0.0 or epsilon_j <= 0.0 or rmin_i <= 0.0 or rmin_j <= 0.0:
        return 0.0
    distance_a = max(distance_nm * 10.0, 1.0e-8)
    rmin_pair = rmin_i + rmin_j
    ratio = rmin_pair / distance_a
    if ratio > 100.0:
        return 1.0e300
    ratio6 = ratio**6
    energy = math.sqrt(epsilon_i * epsilon_j) * (ratio6 * ratio6 - 2.0 * ratio6)
    return energy / max(scnb, 1.0)


def atom_label(topology: md.Topology, index: int) -> str:
    atom = topology.atom(index)
    chain_index = atom.residue.chain.index
    return f"{index + 1}:{chain_index}:{atom.residue.name}{atom.residue.resSeq}:{atom.name}"


def scan(
    workdir: Path,
    pdb: Path,
    xtc: Path,
    prmtop: Path,
    report: dict[str, object],
    cutoff_nm: float,
    pair_limit_kcal: float,
    frame_limit_kcal: float,
    hard_distance_nm: float,
    top_pairs_n: int,
) -> dict[str, object]:
    structure = pmd.load_file(str(prmtop))
    trajectory = md.load_xtc(str(xtc), top=str(pdb))
    if len(structure.atoms) != trajectory.n_atoms:
        raise RuntimeError(
            f"COM.prmtop has {len(structure.atoms)} atoms but trajectory has {trajectory.n_atoms}"
        )

    excluded, one_four = topology_sets(structure)
    synthetic_resnums = synthetic_residue_numbers(workdir)
    synthetic_atoms = {
        atom.index
        for atom in trajectory.topology.atoms
        if atom.residue.resSeq in synthetic_resnums
        and atom.residue.name.upper() in {"ACE", "NME"}
    }
    source_frames = report.get("source_frame_indices", [])
    if not isinstance(source_frames, list) or len(source_frames) != trajectory.n_frames:
        source_frames = list(range(trajectory.n_frames))
    protein_atoms = int(report.get("protein_atoms", trajectory.n_atoms))

    all_records: list[dict[str, object]] = []
    frame_summaries: list[dict[str, object]] = []
    unsafe_frames: set[int] = set()
    unsafe_source_pair = False
    unsafe_synthetic_pair = False
    global_min = (math.inf, -1, -1, -1)
    global_max = (-math.inf, -1, -1, -1)

    for frame_index, xyz in enumerate(np.asarray(trajectory.xyz, dtype=np.float64)):
        candidate_pairs = cKDTree(xyz).query_pairs(cutoff_nm, output_type="ndarray")
        frame_positive = 0.0
        frame_max = 0.0
        frame_min = math.inf
        frame_records: list[dict[str, object]] = []
        for raw_i, raw_j in candidate_pairs:
            i, j = int(raw_i), int(raw_j)
            key = pair_key(i, j)
            if key in excluded:
                continue
            distance = float(np.linalg.norm(xyz[i] - xyz[j]))
            scnb = one_four.get(key, 1.0)
            interaction_class = "1-4" if key in one_four else "nonbonded"
            epsilon_i, rmin_i = atom_lj(structure.atoms[i])
            epsilon_j, rmin_j = atom_lj(structure.atoms[j])
            has_lj = (
                epsilon_i > 0.0 and epsilon_j > 0.0
                and rmin_i > 0.0 and rmin_j > 0.0
            )
            energy = lj_energy_kcal(structure.atoms[i], structure.atoms[j], distance, scnb)
            positive = max(0.0, energy)
            frame_positive += positive
            frame_max = max(frame_max, positive)
            frame_min = min(frame_min, distance)
            synthetic = i in synthetic_atoms or j in synthetic_atoms
            severe = (
                positive >= pair_limit_kcal
                or (has_lj and distance < hard_distance_nm)
            )
            if severe:
                unsafe_frames.add(frame_index)
                unsafe_synthetic_pair |= synthetic
                unsafe_source_pair |= not synthetic
            if severe or positive >= 10.0 or distance < 0.12:
                record = {
                    "prepared_frame_0based": frame_index,
                    "prepared_frame_1based": frame_index + 1,
                    "source_frame_0based": int(source_frames[frame_index]),
                    "atom1_index_1based": i + 1,
                    "atom1": atom_label(trajectory.topology, i),
                    "atom2_index_1based": j + 1,
                    "atom2": atom_label(trajectory.topology, j),
                    "distance_nm": distance,
                    "interaction_class": interaction_class,
                    "scnb": scnb,
                    "lj_kcal_mol": energy,
                    "positive_lj_kcal_mol": positive,
                    "synthetic_cap_involved": synthetic,
                    "atom1_region": "protein" if i < protein_atoms else "ligand",
                    "atom2_region": "protein" if j < protein_atoms else "ligand",
                    "severe": severe,
                }
                frame_records.append(record)
                all_records.append(record)
            if distance < global_min[0]:
                global_min = (distance, frame_index, i, j)
            if positive > global_max[0]:
                global_max = (positive, frame_index, i, j)

        if frame_positive >= frame_limit_kcal:
            unsafe_frames.add(frame_index)
            # Frame-wide excess cannot safely be attributed to a cap unless all
            # substantial repulsive records involve synthetic atoms.
            substantial = [item for item in frame_records if item["positive_lj_kcal_mol"] >= 10.0]
            if substantial and all(bool(item["synthetic_cap_involved"]) for item in substantial):
                unsafe_synthetic_pair = True
            else:
                unsafe_source_pair = True
        frame_summaries.append(
            {
                "prepared_frame_0based": frame_index,
                "prepared_frame_1based": frame_index + 1,
                "source_frame_0based": int(source_frames[frame_index]),
                "minimum_interacting_distance_nm": None if math.isinf(frame_min) else frame_min,
                "maximum_positive_pair_lj_kcal_mol": frame_max,
                "summed_positive_short_range_lj_kcal_mol": frame_positive,
                "unsafe": frame_index in unsafe_frames,
            }
        )

    all_records.sort(key=lambda item: float(item["positive_lj_kcal_mol"]), reverse=True)
    top_records = all_records[:top_pairs_n]
    passed = not unsafe_frames
    return {
        "version": VERSION,
        "pdb": str(pdb),
        "xtc": str(xtc),
        "prmtop": str(prmtop),
        "frames": trajectory.n_frames,
        "atoms": trajectory.n_atoms,
        "cutoff_nm": cutoff_nm,
        "pair_limit_kcal_mol": pair_limit_kcal,
        "frame_limit_kcal_mol": frame_limit_kcal,
        "hard_distance_nm": hard_distance_nm,
        "synthetic_cap_residue_numbers": sorted(synthetic_resnums),
        "synthetic_cap_atom_count": len(synthetic_atoms),
        "validation_passed": passed,
        "unsafe_prepared_frames_0based": sorted(unsafe_frames),
        "unsafe_source_frames_0based": [int(source_frames[i]) for i in sorted(unsafe_frames)],
        "unsafe_source_pair_present": unsafe_source_pair,
        "unsafe_synthetic_cap_pair_present": unsafe_synthetic_pair,
        "global_minimum_interacting_pair": {
            "distance_nm": global_min[0],
            "prepared_frame_0based": global_min[1],
            "source_frame_0based": int(source_frames[global_min[1]]) if global_min[1] >= 0 else None,
            "atom1": atom_label(trajectory.topology, global_min[2]) if global_min[2] >= 0 else None,
            "atom2": atom_label(trajectory.topology, global_min[3]) if global_min[3] >= 0 else None,
        },
        "global_maximum_positive_pair_lj": {
            "lj_kcal_mol": global_max[0],
            "prepared_frame_0based": global_max[1],
            "source_frame_0based": int(source_frames[global_max[1]]) if global_max[1] >= 0 else None,
            "atom1": atom_label(trajectory.topology, global_max[2]) if global_max[2] >= 0 else None,
            "atom2": atom_label(trajectory.topology, global_max[3]) if global_max[3] >= 0 else None,
        },
        "top_pairs": top_records,
        "frame_summaries": frame_summaries,
    }


def write_outputs(workdir: Path, result: dict[str, object], suffix: str = "") -> None:
    json_name = f"steric_preflight{suffix}.json"
    tsv_name = f"steric_clashes{suffix}.tsv"
    (workdir / json_name).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    records = result.get("top_pairs", [])
    if not isinstance(records, list):
        records = []
    columns = [
        "prepared_frame_0based", "prepared_frame_1based", "source_frame_0based",
        "atom1_index_1based", "atom1", "atom2_index_1based", "atom2",
        "distance_nm", "interaction_class", "scnb", "lj_kcal_mol",
        "positive_lj_kcal_mol", "synthetic_cap_involved", "atom1_region",
        "atom2_region", "severe",
    ]
    with (workdir / tsv_name).open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def run_cap_repair(
    workdir: Path,
    repair_script: Path,
    pdb: Path,
    xtc: Path,
    pass_number: int,
) -> None:
    output_pdb = f"{pdb.stem}_stericfix{pass_number}{pdb.suffix}"
    output_xtc = f"{xtc.stem}_stericfix{pass_number}{xtc.suffix}"
    if pass_number == 1:
        hard, trigger, tilt, angle, tilt_step, azimuth = (
            "0.16", "0.20", "60", "2.0", "5.0", "30.0"
        )
    else:
        # The second pass uses the same physically bounded 60-degree tilt but
        # a finer angular mesh and a larger trigger radius.
        hard, trigger, tilt, angle, tilt_step, azimuth = (
            "0.18", "0.22", "60", "1.0", "2.0", "15.0"
        )
    command = [
        sys.executable,
        str(repair_script),
        "--workdir", str(workdir),
        "--pdb", pdb.name,
        "--xtc", xtc.name,
        "--output-pdb", output_pdb,
        "--output-xtc", output_xtc,
        "--allow-no-caps",
        "--backup-and-replace",
        "--trigger-clash-nm", trigger,
        "--hard-clash-nm", hard,
        "--maximum-tilt-deg", tilt,
        "--angle-step-deg", angle,
        "--tilt-step-deg", tilt_step,
        "--tilt-azimuth-step-deg", azimuth,
        "--environment-radius-nm", "1.20",
    ]
    completed = subprocess.run(command, cwd=str(workdir))
    if completed.returncode != 0:
        raise RuntimeError(
            f"Synthetic-cap steric repair pass {pass_number} failed with exit status "
            f"{completed.returncode}"
        )


def main() -> int:
    args = parse_args()
    workdir = Path(args.workdir).expanduser().resolve()
    if not workdir.is_dir():
        raise NotADirectoryError(workdir)
    report = load_report(workdir)
    pdb, xtc = resolve_input_files(workdir, report, args.pdb, args.xtc)
    prmtop = (workdir / args.prmtop).resolve()
    if not prmtop.is_file():
        raise FileNotFoundError(prmtop)

    result = scan(
        workdir, pdb, xtc, prmtop, report,
        args.cutoff_nm, args.pair_limit_kcal, args.frame_limit_kcal,
        args.hard_distance_nm, args.top_pairs,
    )
    write_outputs(workdir, result, "_before")

    repair_history: list[dict[str, object]] = []
    if (
        not bool(result["validation_passed"])
        and not args.no_repair
        and bool(result["unsafe_synthetic_cap_pair_present"])
    ):
        # Repair synthetic-cap pairs first, even when the frame-wide heuristic
        # also flags source atoms.  The repair script is restricted to atoms
        # explicitly marked as synthetic caps.  After each pass, rescan the
        # complete system; any genuine source-only clash remains and causes a
        # hard stop below.
        repair_script = (
            Path(args.repair_script).expanduser().resolve()
            if args.repair_script
            else Path(__file__).with_name("repair_synthetic_caps.py")
        )
        if not repair_script.is_file():
            raise FileNotFoundError(repair_script)
        for pass_number in (1, 2):
            run_cap_repair(workdir, repair_script, pdb, xtc, pass_number)
            repaired = scan(
                workdir, pdb, xtc, prmtop, report,
                args.cutoff_nm, args.pair_limit_kcal, args.frame_limit_kcal,
                args.hard_distance_nm, args.top_pairs,
            )
            repair_history.append(
                {
                    "pass": pass_number,
                    "validation_passed": repaired["validation_passed"],
                    "global_maximum_positive_pair_lj": repaired[
                        "global_maximum_positive_pair_lj"
                    ],
                    "global_minimum_interacting_pair": repaired[
                        "global_minimum_interacting_pair"
                    ],
                    "unsafe_source_pair_present": repaired[
                        "unsafe_source_pair_present"
                    ],
                    "unsafe_synthetic_cap_pair_present": repaired[
                        "unsafe_synthetic_cap_pair_present"
                    ],
                }
            )
            result = repaired
            if bool(result["validation_passed"]):
                break
            if not bool(result["unsafe_synthetic_cap_pair_present"]):
                break

    result["repair_history"] = repair_history
    result["automatic_coordinate_changes_limited_to_synthetic_caps"] = True
    write_outputs(workdir, result, "")
    print(json.dumps(result, indent=2))

    if not bool(result["validation_passed"]):
        if bool(result["unsafe_source_pair_present"]):
            raise RuntimeError(
                "Topology-aware steric preflight found an unsafe pair involving "
                "source protein or ligand atoms. No source coordinates were changed. "
                "See steric_preflight.json and steric_clashes.tsv."
            )
        raise RuntimeError(
            "Synthetic-cap repair did not reduce all topology-aware Lennard-Jones "
            "overlaps below the safety limits. See steric_preflight.json."
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
