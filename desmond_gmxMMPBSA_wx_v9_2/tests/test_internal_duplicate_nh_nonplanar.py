#!/usr/bin/env python3
"""Regression test: non-planar duplicate N-H on an internal peptide residue.

The M4 (5DSG) production system carries a build artifact at MET345: the
residue is peptide-bonded to GLN344 (real C-N bond) yet its backbone N bears
two hydrogens because the source nitrogen is sp3-like.  Across replicas the
freely rotatable CA-C-N-CA torsion samples very different values:

* job_1: omega 164 deg, both O-C-N-H dihedrals non-planar (-134/+108 deg);
* job_2: omega 108 deg (far from amide-like), one H nearly planar (165 deg);
* job_3: omega 175 deg, non-planar hydrogens (-121/+126 deg).

The strict amide-geometry gates previously rejected some of these frames,
the duplicate H was then reinterpreted as a neutral N-terminal cap, and
segment-boundary validation aborted with "Duplicate-H cap on internal
residue ... expected a segment N terminus".  The engine must instead keep
the better-scoring hydrogen and omit the excess one as an internal
correction for every frame.

Coordinates below are the GLN344/MET345 atoms from the failing GROs.
"""
from __future__ import annotations

import importlib.util
import sys
import tempfile
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "app" / "desmond_gmxmmpbsa.py"

# The parser path under test needs only numpy; stub the heavier imports.
md_stub = types.ModuleType("mdtraj")
md_stub.open = lambda *_a, **_k: None
md_stub.iterload = lambda *_a, **_k: iter(())
formats_stub = types.ModuleType("mdtraj.formats")
md_stub.formats = formats_stub
sys.modules.setdefault("mdtraj", md_stub)
sys.modules.setdefault("mdtraj.formats", formats_stub)
nx_stub = types.ModuleType("networkx")
nx_stub.Graph = object
sys.modules.setdefault("networkx", nx_stub)

spec = importlib.util.spec_from_file_location("workflow_engine_nonplanar_v914", ENGINE)
assert spec and spec.loader
engine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = engine
spec.loader.exec_module(engine)
assert engine.__version__ == "9.1.4"


def atom_line(resid: int, resname: str, name: str, number: int, x: float, y: float, z: float) -> str:
    return f"{resid:5d}{resname:<5.5s}{name:<5.5s}{number:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def write_gro(path: Path, lines: list[str]) -> None:
    path.write_text(
        "regression non-planar internal duplicate-H test\n"
        f"{len(lines)}\n" + "".join(line + "\n" for line in lines)
        + "  10.00000  10.00000  10.00000\n",
        encoding="utf-8",
    )


# job_1 frame: omega 164.4 deg, O-C-N-H1 -133.8 deg, O-C-N-H2 +108.0 deg.
ATOMS_JOB1 = [
    atom_line(344, "GLN", "N", 1, -0.033, -1.341, -3.050),
    atom_line(344, "GLN", "CA", 2, 0.059, -1.320, -2.936),
    atom_line(344, "GLN", "C", 3, 0.011, -1.199, -2.858),
    atom_line(344, "GLN", "O", 4, -0.009, -1.090, -2.911),
    atom_line(344, "GLN", "CB", 5, 0.208, -1.299, -2.984),
    atom_line(345, "MET", "N", 6, -0.026, -1.216, -2.724),
    atom_line(345, "MET", "CA", 7, -0.038, -1.094, -2.641),
    atom_line(345, "MET", "C", 8, 0.103, -1.030, -2.626),
    atom_line(345, "MET", "O", 9, 0.202, -1.098, -2.635),
    atom_line(345, "MET", "CB", 10, -0.092, -1.124, -2.501),
    atom_line(345, "MET", "H1", 11, 0.046, -1.261, -2.669),
    atom_line(345, "MET", "H2", 12, -0.107, -1.276, -2.726),
]

# job_2 frame: omega 107.8 deg (non-amide-like), O-C-N-H1 +165.1 deg,
# O-C-N-H2 +57.1 deg.  The old omega gate rejected this frame outright.
ATOMS_JOB2 = [
    atom_line(344, "GLN", "N", 1, -0.261, -1.165, -3.233),
    atom_line(344, "GLN", "CA", 2, -0.195, -1.102, -3.117),
    atom_line(344, "GLN", "C", 3, -0.299, -1.030, -3.030),
    atom_line(344, "GLN", "O", 4, -0.360, -0.932, -3.077),
    atom_line(344, "GLN", "CB", 5, -0.089, -1.007, -3.163),
    atom_line(345, "MET", "N", 6, -0.323, -1.071, -2.901),
    atom_line(345, "MET", "CA", 7, -0.273, -0.983, -2.789),
    atom_line(345, "MET", "C", 8, -0.124, -0.965, -2.794),
    atom_line(345, "MET", "O", 9, -0.049, -1.058, -2.811),
    atom_line(345, "MET", "CB", 10, -0.315, -1.030, -2.649),
    atom_line(345, "MET", "H1", 11, -0.281, -1.163, -2.901),
    atom_line(345, "MET", "H2", 12, -0.422, -1.091, -2.892),
]


def run_case(label: str, atoms: list[str], expected_removed: str, expected_retained: str) -> None:
    with tempfile.TemporaryDirectory() as directory:
        gro = Path(directory) / f"{label}.gro"
        write_gro(gro, atoms)
        _title, all_atoms, box = engine.parse_gro(gro)
        raw_blocks = engine.contiguous_residue_blocks(all_atoms, engine.AA_NAMES)
        residues, caps, corrections, dropped = engine.separate_noncontiguous_hydrogen_caps(
            raw_blocks, box
        )
        assert dropped == []
        assert [r.source_resid for r in residues] == [344, 345]
        # MET345 must be an internal correction, not a terminal N_H cap.
        assert (345, "MET") not in caps, (label, caps)
        assert len(corrections) == 1, (label, corrections)
        correction = corrections[0]
        assert correction.residue_key == (345, "MET"), label
        assert correction.previous_residue_key == (344, "GLN"), label
        assert correction.removed_atom.source_name == expected_removed, label
        assert correction.retained_atom is not None, label
        assert correction.retained_atom.source_name == expected_retained, label
        assert correction.retained_o_c_n_h_dihedral_deg is not None, label
        # The kept atom set no longer contains the excess hydrogen.
        met = residues[1]
        kept_names = [atom.source_name for atom in met.atoms]
        assert expected_removed not in kept_names, label
        assert expected_retained in kept_names, label
        # Segment detection sees one peptide-connected segment and cap
        # assignment must not raise the former segment-N-terminus error.
        segments = engine.detect_chain_segments(residues, box)
        assert segments == [[0, 1]], (label, segments)
        assignments = engine.assign_hydrogen_caps_to_segment_boundaries(
            residues, segments, caps
        )
        assert assignments == [], label


# job_1 and job_3: neither hydrogen is planar; H1 scores better in job_1.
run_case("job1", ATOMS_JOB1, "H2", "H1")
# job_2: omega is far from amide-like; H1 (165 deg) is retained over H2 (57 deg).
run_case("job2", ATOMS_JOB2, "H2", "H1")

print("test_internal_duplicate_nh_nonplanar: OK")
