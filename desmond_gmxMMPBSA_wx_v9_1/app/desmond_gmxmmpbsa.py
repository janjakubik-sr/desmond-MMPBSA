#!/usr/bin/env python3
"""Prepare and run an open-source gmx_MMPBSA rescoring workflow.

The miniapp accepts a full-system GRO/XTC pair converted from a Desmond
trajectory plus an Antechamber-compatible ligand MOL2. It rebuilds a dry
protein-ligand complex in Amber ff19SB/GAFF2 atom order, removes periodic
imaging, fits the trajectory, retains every Nth input frame, generates Amber
and GROMACS topologies through AmberTools/ParmEd, and runs gmx_MMPBSA.

Subcommands
-----------
prepare  Create a topology-ordered dry PDB/GRO/XTC, index, and build inputs.
run      Build/validate topologies and execute gmx_MMPBSA in a prepared folder.
all      Perform prepare followed by run.
summarize Parse an existing FINAL_RESULTS_MMPBSA.dat file.

Preparation dependencies: numpy, mdtraj, networkx.
Run dependencies: AmberTools (parmchk2, tleap), ParmEd, GROMACS, gmx_MMPBSA.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import mdtraj as md
import networkx as nx
import numpy as np

__version__ = "9.1.0"

AA_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN",
    "GLU", "GLH", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU",
    "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR",
    "VAL", "ACE", "NME", "NMA",
}
VARIANT_NAMES = AA_NAMES | {"HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN"}

# Desmond/Schrodinger residue aliases that describe standard Amber cap
# chemistry but use different residue or atom nomenclature.  NMA is the
# methylamide cap represented by ff19SB as NME.
SOURCE_TEMPLATE_ALIASES = {"NMA": "NME"}
CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
COVALENT_RADII_A = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20,
    "I": 1.39, "SI": 1.11,
}
ATOMIC_NUMBERS = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "SI": 14,
    "P": 15, "S": 16, "CL": 17, "BR": 35, "I": 53,
}


@dataclass
class GroAtom:
    source_index: int
    source_resid: int
    source_resname: str
    source_name: str
    atom_number: int
    xyz_nm: np.ndarray
    element: str


@dataclass
class ResidueBlock:
    source_resid: int
    source_resname: str
    atoms: list[GroAtom]
    chain_index: int = -1
    output_resid: int = -1
    amber_template: str = ""
    output_resname: str = ""


@dataclass
class OutputAtom:
    output_index: int
    source_index: int
    atom_name: str
    residue_name: str
    residue_number: int
    chain_id: str
    element: str
    record_name: str
    is_protein: bool


@dataclass
class HydrogenCap:
    """A Desmond neutral hydrogen used to cap a protein chain break."""

    kind: str
    residue_key: tuple[int, str]
    source_atom: GroAtom
    anchor_atom: GroAtom
    plane_atom: GroAtom
    distance_nm: float


@dataclass
class HydrogenCapAssignment:
    """One Desmond terminal hydrogen assigned to a protein segment boundary.

    C-terminal carbonyl-H and N-terminal duplicate-H caps are handled independently.
    This supports mixed prepared structures in which one side of a chain break
    is already represented by a standard ACE or NME residue.
    """

    residue_index: int
    segment_index: int
    boundary: str
    replacement_template: str
    cap: HydrogenCap


@dataclass
class InternalBackboneHydrogenCorrection:
    """One excess N-bound hydrogen removed from an internal peptide residue.

    Some Desmond/VMD GRO exports contain a residue with two backbone-N
    hydrogens even though its N is covalently joined to the preceding peptide
    carbonyl.  Such an atom set cannot be represented by a standard internal
    ff19SB residue.  The chemically planar amide hydrogen is retained and the
    other source hydrogen is omitted from the Amber rescoring model.
    """

    residue_key: tuple[int, str]
    previous_residue_key: tuple[int, str]
    removed_atom: GroAtom
    retained_atom: GroAtom | None
    c_n_distance_nm: float
    ca_c_n_angle_deg: float
    c_n_ca_angle_deg: float
    omega_deg: float
    retained_o_c_n_h_dihedral_deg: float | None
    removed_o_c_n_h_dihedral_deg: float


@dataclass
class SyntheticCapSpec:
    """Instructions for regenerating one standard ACE/NME cap per frame."""

    template_name: str
    output_residue: int
    chain_id: str
    output_indices: dict[str, int]
    origin_source_index: int
    direction_source_index: int
    plane_source_index: int
    replaced_source_index: int
    source_resid: int
    source_resname: str


# Local ff19SB-compatible cap coordinates in nm. They were measured from
# explicit ACE/NME caps in a successfully mapped Desmond receptor system and
# are expressed in a residue-local orthonormal frame. Per-frame placement is
# controlled by the original Desmond terminal hydrogen direction, so the caps
# follow the source trajectory rather than remaining fixed in laboratory space.
NME_LOCAL_COORDS_NM: dict[str, np.ndarray] = {
    "N": np.array([0.122593, 0.000000, 0.000000]),
    "H": np.array([0.172008, -0.086993, -0.000607]),
    "C": np.array([0.198030, 0.124054, -0.000011]),
    "H1": np.array([0.296216, 0.104813, -0.000381]),
    "H2": np.array([0.177694, 0.178215, 0.081439]),
    "H3": np.array([0.177050, 0.179141, -0.080628]),
}

ACE_LOCAL_COORDS_NM: dict[str, np.ndarray] = {
    "H1": np.array([0.292250, -0.137356, 0.000225]),
    "CH3": np.array([0.192487, -0.140138, -0.000484]),
    "H2": np.array([0.162955, -0.191340, -0.081424]),
    "H3": np.array([0.163171, -0.191824, 0.080554]),
    "C": np.array([0.134227, 0.000000, 0.000000]),
    "O": np.array([0.194834, 0.092222, 0.000275]),
}


def natural_key(text: str) -> list[object]:
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", text)]


def infer_element(atom_name: str, residue_name: str = "") -> str:
    """Infer a chemical element from a PDB/GRO/MOL2 atom name.

    Protein atom CA is carbon, while a residue/atom named CL is chlorine.
    MOL2 force-field type names are deliberately not used because GAFF types
    such as ca denote aromatic carbon, not calcium.
    """
    name = re.sub(r"^[0-9]+", "", atom_name.strip().upper())
    residue = residue_name.strip().upper()
    if residue == "CL" and name.startswith("CL"):
        return "CL"
    if residue == "BR" and name.startswith("BR"):
        return "BR"
    for symbol in ("CL", "BR", "SI"):
        if name.startswith(symbol) and name not in {"CA", "CB", "CG", "CD", "CE", "CZ"}:
            return symbol
    return name[0] if name else "X"


def parse_gro(path: Path) -> tuple[str, list[GroAtom], np.ndarray]:
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        title = handle.readline().rstrip("\n")
        count_line = handle.readline()
        try:
            n_atoms = int(count_line.strip())
        except ValueError as exc:
            raise ValueError(f"Invalid GRO atom-count line in {path}: {count_line!r}") from exc
        atoms: list[GroAtom] = []
        for index in range(n_atoms):
            line = handle.readline().rstrip("\n")
            if len(line) < 44:
                raise ValueError(f"Malformed GRO atom line {index + 3}: {line!r}")
            resid = int(line[0:5])
            resname = line[5:10].strip().upper()
            atom_name = line[10:15].strip()
            atom_number = int(line[15:20])
            xyz = np.array(
                [float(line[20:28]), float(line[28:36]), float(line[36:44])],
                dtype=float,
            )
            atoms.append(
                GroAtom(
                    source_index=index,
                    source_resid=resid,
                    source_resname=resname,
                    source_name=atom_name,
                    atom_number=atom_number,
                    xyz_nm=xyz,
                    element=infer_element(atom_name, resname),
                )
            )
        box_fields = [float(value) for value in handle.readline().split()]
    if len(box_fields) == 3:
        box = np.diag(box_fields)
    elif len(box_fields) == 9:
        # GRO order: v1x v2y v3z v1y v1z v2x v2z v3x v3y
        box = np.array(
            [
                [box_fields[0], box_fields[3], box_fields[4]],
                [box_fields[5], box_fields[1], box_fields[6]],
                [box_fields[7], box_fields[8], box_fields[2]],
            ],
            dtype=float,
        )
    else:
        raise ValueError(f"Expected 3 or 9 GRO box values; found {len(box_fields)}")
    return title, atoms, box


def minimum_image(vector: np.ndarray, box: np.ndarray) -> np.ndarray:
    inverse = np.linalg.inv(box)
    fractional = vector @ inverse
    fractional -= np.round(fractional)
    return fractional @ box


def contiguous_residue_blocks(atoms: Sequence[GroAtom], allowed: set[str]) -> list[ResidueBlock]:
    blocks: list[ResidueBlock] = []
    current_key: tuple[int, str] | None = None
    current: ResidueBlock | None = None
    for atom in atoms:
        if atom.source_resname not in allowed:
            continue
        key = (atom.source_resid, atom.source_resname)
        if current_key != key:
            current = ResidueBlock(key[0], key[1], [])
            blocks.append(current)
            current_key = key
        assert current is not None
        current.atoms.append(atom)
    return blocks


def backbone_n_hydrogens(
    residue: ResidueBlock,
    box: np.ndarray,
    minimum_nm: float = 0.07,
    maximum_nm: float = 0.14,
) -> list[GroAtom]:
    """Return hydrogens covalently attached to the peptide-backbone N.

    Desmond/GRO exports use several naming schemes for terminal amine
    hydrogens (``H``, ``1H``/``2H``/``3H``, or ``H1``/``H2``/``H3``).
    Geometry is therefore more reliable than atom names.  The distance test is
    evaluated with minimum-image vectors so wrapped terminal atoms are handled
    correctly.
    """

    nitrogen = next(
        (atom for atom in residue.atoms if atom.source_name.upper() == "N"),
        None,
    )
    if nitrogen is None:
        return []
    result: list[GroAtom] = []
    for atom in residue.atoms:
        if atom.element.upper() != "H" or atom is nitrogen:
            continue
        distance = float(
            np.linalg.norm(minimum_image(atom.xyz_nm - nitrogen.xyz_nm, box))
        )
        if minimum_nm <= distance <= maximum_nm:
            result.append(atom)
    result.sort(key=lambda atom: atom.source_index)
    return result


def angle_degrees(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Return the A-B-C angle in degrees."""

    left = a - b
    right = c - b
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator == 0.0:
        raise ValueError("Cannot calculate an angle from coincident atoms")
    cosine = float(np.clip(np.dot(left, right) / denominator, -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def dihedral_degrees(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    d: np.ndarray,
) -> float:
    """Return the signed A-B-C-D dihedral in degrees."""

    b0 = -(b - a)
    b1 = c - b
    b2 = d - c
    norm_b1 = float(np.linalg.norm(b1))
    if norm_b1 == 0.0:
        raise ValueError("Cannot calculate a dihedral across a zero-length bond")
    b1 = b1 / norm_b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if float(np.linalg.norm(v)) == 0.0 or float(np.linalg.norm(w)) == 0.0:
        raise ValueError("Cannot calculate a dihedral from collinear atoms")
    return math.degrees(math.atan2(float(np.dot(np.cross(b1, v), w)), float(np.dot(v, w))))


def internal_duplicate_nh_correction(
    previous: ResidueBlock | None,
    current: ResidueBlock,
    n_hydrogens: Sequence[GroAtom],
    ordinary_count: int,
    box: np.ndarray,
) -> InternalBackboneHydrogenCorrection | None:
    """Classify an extra backbone-N hydrogen on an internal peptide residue.

    A two-hydrogen backbone N is normally interpreted as a neutral N terminus.
    That interpretation is wrong when the residue is already joined to the
    preceding carbonyl.  We require peptide-like C-N distance, backbone bond
    angles, and a trans omega dihedral.  For an ordinary amino acid, the
    retained hydrogen is the one whose O-C-N-H dihedral is closest to 180
    degrees, i.e. the planar amide hydrogen.  Ambiguous cases are left for the
    terminal-cap logic or a strict mapping error rather than guessed.
    """

    if previous is None or len(n_hydrogens) != ordinary_count + 1:
        return None

    # The provisional current residue may still contain two atoms with the
    # same literal name (commonly two ``H`` records).  Build a lookup only for
    # the unique heavy atoms needed for peptide geometry rather than calling
    # atom_lookup(), whose duplicate-name guard is intentionally stricter.
    previous_lookup = {
        atom.source_name: atom
        for atom in previous.atoms
        if atom.source_name in {"CA", "C", "O", "OXT"}
    }
    current_lookup = {
        atom.source_name: atom
        for atom in current.atoms
        if atom.source_name in {"N", "CA"}
    }
    if "OXT" in previous_lookup:
        return None
    required_previous = {"CA", "C", "O"}
    required_current = {"N", "CA"}
    if not required_previous.issubset(previous_lookup) or not required_current.issubset(current_lookup):
        return None

    carbon = previous_lookup["C"].xyz_nm
    previous_ca = carbon + minimum_image(previous_lookup["CA"].xyz_nm - carbon, box)
    oxygen = carbon + minimum_image(previous_lookup["O"].xyz_nm - carbon, box)
    nitrogen = carbon + minimum_image(current_lookup["N"].xyz_nm - carbon, box)
    current_ca = nitrogen + minimum_image(current_lookup["CA"].xyz_nm - nitrogen, box)

    c_n_distance = float(np.linalg.norm(nitrogen - carbon))
    if not 0.115 <= c_n_distance <= 0.170:
        return None

    ca_c_n_angle = angle_degrees(previous_ca, carbon, nitrogen)
    c_n_ca_angle = angle_degrees(carbon, nitrogen, current_ca)
    omega = dihedral_degrees(previous_ca, carbon, nitrogen, current_ca)
    omega_deviation = abs(180.0 - abs(omega))
    if not (90.0 <= ca_c_n_angle <= 145.0):
        return None
    if not (95.0 <= c_n_ca_angle <= 145.0):
        return None
    if omega_deviation > 45.0:
        return None

    if ordinary_count == 0:
        # Internal proline has no N-H.  One source N-H on a peptide-connected
        # proline is therefore the excess atom.
        removed = n_hydrogens[0]
        removed_xyz = nitrogen + minimum_image(removed.xyz_nm - nitrogen, box)
        removed_dihedral = dihedral_degrees(oxygen, carbon, nitrogen, removed_xyz)
        return InternalBackboneHydrogenCorrection(
            residue_key=(current.source_resid, current.source_resname),
            previous_residue_key=(previous.source_resid, previous.source_resname),
            removed_atom=removed,
            retained_atom=None,
            c_n_distance_nm=c_n_distance,
            ca_c_n_angle_deg=ca_c_n_angle,
            c_n_ca_angle_deg=c_n_ca_angle,
            omega_deg=omega,
            retained_o_c_n_h_dihedral_deg=None,
            removed_o_c_n_h_dihedral_deg=removed_dihedral,
        )

    if ordinary_count != 1 or len(n_hydrogens) != 2:
        return None

    scored: list[tuple[float, float, GroAtom]] = []
    for hydrogen in n_hydrogens:
        hydrogen_xyz = nitrogen + minimum_image(hydrogen.xyz_nm - nitrogen, box)
        dihedral = dihedral_degrees(oxygen, carbon, nitrogen, hydrogen_xyz)
        score = abs(180.0 - abs(dihedral))
        scored.append((score, dihedral, hydrogen))
    scored.sort(key=lambda item: (item[0], item[2].source_index))

    best_score, retained_dihedral, retained = scored[0]
    second_score, removed_dihedral, removed = scored[1]
    # A standard peptide amide hydrogen is nearly trans to the carbonyl O.
    # Require both a good retained geometry and a clear distinction between
    # the two candidates.  This prevents arbitrary deletion at ambiguous
    # termini or malformed structures.
    if best_score > 35.0 or second_score - best_score < 25.0:
        return None

    return InternalBackboneHydrogenCorrection(
        residue_key=(current.source_resid, current.source_resname),
        previous_residue_key=(previous.source_resid, previous.source_resname),
        removed_atom=removed,
        retained_atom=retained,
        c_n_distance_nm=c_n_distance,
        ca_c_n_angle_deg=ca_c_n_angle,
        c_n_ca_angle_deg=c_n_ca_angle,
        omega_deg=omega,
        retained_o_c_n_h_dihedral_deg=retained_dihedral,
        removed_o_c_n_h_dihedral_deg=removed_dihedral,
    )


def separate_noncontiguous_hydrogen_caps(
    blocks: Sequence[ResidueBlock], box: np.ndarray
) -> tuple[
    list[ResidueBlock],
    dict[tuple[int, str], list[HydrogenCap]],
    list[InternalBackboneHydrogenCorrection],
]:
    """Extract neutral chain-break hydrogens from protein residues.

    Maestro/VMD GRO exports may place terminal hydrogens either inside the
    main residue block or in a later, non-contiguous block with the same
    residue number/name.  Two chemically equivalent naming conventions are
    supported for an N-terminal neutral amine: two atoms both named ``H`` or
    two N-bound hydrogens named, for example, ``1H`` and ``2H``.

    A neutral terminal amine is converted to an ACE-capped peptide terminus by
    retaining the normal peptide N-H (none for internal PRO) and using the
    other N-H direction as the per-frame ACE placement anchor.  A fully
    protonated free N terminus (three N-H atoms for ordinary residues, two for
    PRO) is retained and mapped to the ff19SB N-terminal template.
    """

    atoms_by_key: dict[tuple[int, str], list[GroAtom]] = defaultdict(list)
    key_order: list[tuple[int, str]] = []
    for block in blocks:
        key = (block.source_resid, block.source_resname)
        if key not in atoms_by_key:
            key_order.append(key)
        atoms_by_key[key].extend(block.atoms)

    normalized: list[ResidueBlock] = []
    caps: dict[tuple[int, str], list[HydrogenCap]] = defaultdict(list)
    internal_corrections: list[InternalBackboneHydrogenCorrection] = []

    for key_position, key in enumerate(key_order):
        residue_atoms = list(atoms_by_key[key])

        # Desmond protein-preparation exports use both HXT and HC for a
        # neutral hydrogen attached directly to a terminal backbone carbonyl
        # carbon.  Detect this site by chemistry rather than by one literal
        # atom name.  Ordinary peptide residues have no C-bound hydrogen at
        # their backbone carbonyl; ACE/NME are excluded because their carbon
        # naming describes different cap chemistry.
        carbonyl = next(
            (atom for atom in residue_atoms if atom.source_name.upper() == "C"),
            None,
        )
        hxt_atoms: list[GroAtom] = []
        if carbonyl is not None and SOURCE_TEMPLATE_ALIASES.get(key[1].upper(), key[1].upper()) not in {"ACE", "NME"}:
            for atom in residue_atoms:
                if atom.element.upper() != "H":
                    continue
                distance = float(
                    np.linalg.norm(minimum_image(atom.xyz_nm - carbonyl.xyz_nm, box))
                )
                if 0.07 <= distance <= 0.14:
                    hxt_atoms.append(atom)
        if len(hxt_atoms) > 1:
            raise ValueError(
                f"Multiple hydrogens bonded to backbone carbonyl C in residue {key}: "
                f"{[atom.source_name for atom in hxt_atoms]}"
            )
        if hxt_atoms:
            residue_atoms.remove(hxt_atoms[0])

        provisional = ResidueBlock(key[0], key[1], list(residue_atoms))
        n_hydrogens = backbone_n_hydrogens(provisional, box)
        is_proline = key[1].upper() == "PRO"
        ordinary_count = 0 if is_proline else 1
        protonated_count = 2 if is_proline else 3
        extra_n_h: GroAtom | None = None

        if len(n_hydrogens) == ordinary_count + 1:
            previous_block = None
            if key_position > 0:
                previous_key = key_order[key_position - 1]
                previous_block = ResidueBlock(
                    previous_key[0],
                    previous_key[1],
                    list(atoms_by_key[previous_key]),
                )
            internal_correction = internal_duplicate_nh_correction(
                previous_block,
                provisional,
                n_hydrogens,
                ordinary_count,
                box,
            )
            if internal_correction is not None:
                # The residue is peptide-connected to the previous residue,
                # so this is not an N-terminal ACE placement vector.  Retain
                # the chemically planar amide H (if any) and omit only the
                # excess internal N-H from the ff19SB rescoring model.
                residue_atoms.remove(internal_correction.removed_atom)
                internal_corrections.append(internal_correction)
            else:
                # Neutral free NH/NH2 terminus.  Use one hydrogen as the
                # direction of the missing acyl C and retain the ordinary
                # peptide N-H when the residue is not proline.  Prefer a
                # conventional atom named H; otherwise retain the earlier
                # source atom.  At a true terminus the hydrogens are
                # chemically equivalent before ACE placement.
                if ordinary_count == 0:
                    extra_n_h = n_hydrogens[0]
                else:
                    conventional = [
                        atom for atom in n_hydrogens if atom.source_name.upper() == "H"
                    ]
                    retained = conventional[0] if conventional else n_hydrogens[0]
                    extra_n_h = next(
                        atom for atom in n_hydrogens if atom is not retained
                    )
                residue_atoms.remove(extra_n_h)
        elif len(n_hydrogens) in {ordinary_count, protonated_count}:
            pass
        elif len(n_hydrogens) > protonated_count:
            labels = [atom.source_name for atom in n_hydrogens]
            raise ValueError(
                f"Too many hydrogens ({labels}) bonded to backbone N in residue {key}"
            )
        # Counts below the ordinary template expectation are left for the
        # strict template mapper/topology reconciliation to diagnose.

        # All remaining names must be unique for deterministic template
        # mapping.  Different terminal-hydrogen names (1H/2H etc.) are unique
        # and are normalized later by map_protein_names().
        seen_names: set[str] = set()
        duplicates: list[str] = []
        for atom in residue_atoms:
            if atom.source_name in seen_names:
                duplicates.append(atom.source_name)
            seen_names.add(atom.source_name)
        if duplicates:
            raise ValueError(
                f"Unsupported duplicate atom names in residue {key}: "
                f"{sorted(set(duplicates))}"
            )

        main = ResidueBlock(key[0], key[1], residue_atoms)
        lookup = atom_lookup(main)

        if hxt_atoms:
            atom = hxt_atoms[0]
            if atom.element != "H" or "C" not in lookup or "O" not in lookup:
                raise ValueError(f"Carbonyl-H residue {key} lacks valid H/C/O cap geometry")
            anchor = lookup["C"]
            plane = lookup["O"]
            distance = float(
                np.linalg.norm(minimum_image(atom.xyz_nm - anchor.xyz_nm, box))
            )
            if not 0.07 <= distance <= 0.14:
                raise ValueError(
                    f"Carbonyl-bound H in residue {key} is {distance:.4f} nm from C; "
                    "expected a covalent C-H distance"
                )
            caps[key].append(
                HydrogenCap(
                    kind="C_HXT",
                    residue_key=key,
                    source_atom=atom,
                    anchor_atom=anchor,
                    plane_atom=plane,
                    distance_nm=distance,
                )
            )

        if extra_n_h is not None:
            if extra_n_h.element != "H" or "N" not in lookup or "CA" not in lookup:
                raise ValueError(
                    f"Neutral N-terminal residue {key} lacks valid H/N/CA geometry"
                )
            anchor = lookup["N"]
            plane = lookup["CA"]
            distance = float(
                np.linalg.norm(minimum_image(extra_n_h.xyz_nm - anchor.xyz_nm, box))
            )
            if not 0.07 <= distance <= 0.14:
                raise ValueError(
                    f"Extra N-H in residue {key} is {distance:.4f} nm from N; "
                    "expected a covalent N-H distance"
                )
            caps[key].append(
                HydrogenCap(
                    kind="N_H",
                    residue_key=key,
                    source_atom=extra_n_h,
                    anchor_atom=anchor,
                    plane_atom=plane,
                    distance_nm=distance,
                )
            )

        normalized.append(main)

    return normalized, dict(caps), internal_corrections

def assign_hydrogen_caps_to_segment_boundaries(
    residues: Sequence[ResidueBlock],
    segments: Sequence[Sequence[int]],
    caps: dict[tuple[int, str], list[HydrogenCap]],
) -> list[HydrogenCapAssignment]:
    """Assign terminal hydrogens independently to segment boundaries.

    Earlier versions required every chain break to contain a matched
    C-terminal carbonyl-H *and* N-terminal duplicate H.  Desmond-prepared structures
    can legitimately mix representations, for example an explicit NME at the
    end of one segment followed by a duplicate-H capped amino acid at the
    beginning of the next.  Requiring a symmetric pair rejects such systems.

    Each detected cap is therefore validated against its own segment boundary:

    * ``C_HXT`` (source HXT/HC) is allowed only on the last residue of a segment and is
      replaced by a standard ff19SB NME residue.
    * ``N_H`` is allowed only on the first residue of a segment and is
      replaced by a standard ff19SB ACE residue.

    Existing source ACE/NME residues are retained unchanged.  A segment may
    consequently need zero, one, or two synthetic cap replacements, and the
    first/last protein segment is handled in exactly the same way as an
    internal chain break.
    """

    key_to_index: dict[tuple[int, str], int] = {}
    for index, residue in enumerate(residues):
        key = (residue.source_resid, residue.source_resname)
        if key in key_to_index:
            raise ValueError(
                f"Protein residue key {key} occurs more than once after normalization"
            )
        key_to_index[key] = index

    residue_to_segment: dict[int, int] = {}
    segment_starts: set[int] = set()
    segment_ends: set[int] = set()
    for segment_index, segment in enumerate(segments):
        if not segment:
            raise ValueError(f"Protein segment {segment_index + 1} is empty")
        segment_starts.add(segment[0])
        segment_ends.add(segment[-1])
        for residue_index in segment:
            if residue_index in residue_to_segment:
                raise ValueError(
                    f"Protein residue index {residue_index} belongs to multiple segments"
                )
            residue_to_segment[residue_index] = segment_index

    assignments: list[HydrogenCapAssignment] = []
    assigned_boundaries: set[tuple[int, str]] = set()
    all_caps = [cap for residue_caps in caps.values() for cap in residue_caps]

    for cap in all_caps:
        if cap.residue_key not in key_to_index:
            raise ValueError(
                f"Hydrogen cap {cap.kind} refers to unknown residue {cap.residue_key}"
            )
        residue_index = key_to_index[cap.residue_key]
        segment_index = residue_to_segment[residue_index]
        residue = residues[residue_index]

        if cap.kind == "C_HXT":
            boundary = "C"
            replacement = "NME"
            if residue_index not in segment_ends:
                raise ValueError(
                    f"Carbonyl-H cap on internal residue {cap.residue_key}; expected a segment C terminus"
                )
            if SOURCE_TEMPLATE_ALIASES.get(residue.source_resname, residue.source_resname) == "NME":
                raise ValueError(
                    f"Residue {cap.residue_key} is already NME but also contains a carbonyl-bound H"
                )
        elif cap.kind == "N_H":
            boundary = "N"
            replacement = "ACE"
            if residue_index not in segment_starts:
                raise ValueError(
                    f"Duplicate-H cap on internal residue {cap.residue_key}; "
                    "expected a segment N terminus"
                )
            if SOURCE_TEMPLATE_ALIASES.get(residue.source_resname, residue.source_resname) == "ACE":
                raise ValueError(
                    f"Residue {cap.residue_key} is already ACE but also contains a duplicate N-H"
                )
        else:
            raise ValueError(f"Unsupported hydrogen-cap kind {cap.kind!r}")

        boundary_key = (residue_index, boundary)
        if boundary_key in assigned_boundaries:
            raise ValueError(
                f"Multiple hydrogen caps assigned to residue {cap.residue_key} "
                f"at its {boundary} terminus"
            )
        assigned_boundaries.add(boundary_key)
        assignments.append(
            HydrogenCapAssignment(
                residue_index=residue_index,
                segment_index=segment_index,
                boundary=boundary,
                replacement_template=replacement,
                cap=cap,
            )
        )

    assignments.sort(
        key=lambda assignment: (
            assignment.segment_index,
            0 if assignment.boundary == "N" else 1,
            assignment.cap.source_atom.source_index,
        )
    )
    return assignments


def atom_lookup(residue: ResidueBlock) -> dict[str, GroAtom]:
    lookup: dict[str, GroAtom] = {}
    for atom in residue.atoms:
        if atom.source_name in lookup:
            raise ValueError(
                f"Duplicate atom name {atom.source_name!r} in residue "
                f"{residue.source_resid} {residue.source_resname}"
            )
        lookup[atom.source_name] = atom
    return lookup


def detect_chain_segments(
    residues: list[ResidueBlock], box: np.ndarray, cutoff_nm: float = 0.20
) -> list[list[int]]:
    segments: list[list[int]] = []
    for index, residue in enumerate(residues):
        if index == 0:
            segments.append([index])
            continue
        previous = residues[index - 1]
        prev_atoms = atom_lookup(previous)
        curr_atoms = atom_lookup(residue)
        connected = False
        if "C" in prev_atoms and "N" in curr_atoms:
            delta = minimum_image(curr_atoms["N"].xyz_nm - prev_atoms["C"].xyz_nm, box)
            connected = float(np.linalg.norm(delta)) <= cutoff_nm
        if connected:
            segments[-1].append(index)
        else:
            segments.append([index])
    if len(segments) > len(CHAIN_IDS):
        raise ValueError(f"Too many protein chains ({len(segments)}); maximum {len(CHAIN_IDS)}")
    for chain_index, segment in enumerate(segments):
        for residue_index in segment:
            residues[residue_index].chain_index = chain_index
    return segments


def load_ffxml(path: Path) -> tuple[dict[str, list[str]], dict[str, list[tuple[str, str]]]]:
    root = ET.parse(path).getroot()
    residue_root = root.find("Residues")
    if residue_root is None:
        raise ValueError(f"No <Residues> section in {path}")
    atoms: dict[str, list[str]] = {}
    bonds: dict[str, list[tuple[str, str]]] = {}
    for residue in residue_root.findall("Residue"):
        name = residue.attrib["name"]
        atoms[name] = [atom.attrib["name"] for atom in residue.findall("Atom")]
        bonds[name] = [
            (bond.attrib["atomName1"], bond.attrib["atomName2"])
            for bond in residue.findall("Bond")
        ]
    return atoms, bonds


def detect_disulfides(
    residues: list[ResidueBlock], box: np.ndarray, cutoff_nm: float = 0.23
) -> list[tuple[int, int, float]]:
    candidates: list[tuple[int, GroAtom]] = []
    for index, residue in enumerate(residues):
        if residue.source_resname not in {"CYS", "CYX", "CYM"}:
            continue
        lookup = atom_lookup(residue)
        if "SG" in lookup and "HG" not in lookup:
            candidates.append((index, lookup["SG"]))
    edges: list[tuple[float, int, int]] = []
    for left in range(len(candidates)):
        index_left, atom_left = candidates[left]
        for right in range(left + 1, len(candidates)):
            index_right, atom_right = candidates[right]
            delta = minimum_image(atom_right.xyz_nm - atom_left.xyz_nm, box)
            distance = float(np.linalg.norm(delta))
            if distance <= cutoff_nm:
                edges.append((distance, index_left, index_right))
    edges.sort()
    used: set[int] = set()
    pairs: list[tuple[int, int, float]] = []
    for distance, left, right in edges:
        if left in used or right in used:
            continue
        used.update((left, right))
        pairs.append((left, right, distance))
    unmatched = [index for index, _atom in candidates if index not in used]
    if unmatched:
        labels = [residues[index].source_resid for index in unmatched]
        raise ValueError(
            "Cysteines without HG could not be paired into disulfides: " + ", ".join(map(str, labels))
        )
    return pairs


def base_variant(residue: ResidueBlock, cyx_indices: set[int], residue_index: int) -> str:
    name = SOURCE_TEMPLATE_ALIASES.get(residue.source_resname, residue.source_resname)
    names = {atom.source_name for atom in residue.atoms}
    if name in {"HIS", "HID", "HIE", "HIP"}:
        has_hd1 = "HD1" in names
        has_he2 = "HE2" in names
        if has_hd1 and has_he2:
            return "HIP"
        if has_hd1:
            return "HID"
        if has_he2:
            return "HIE"
        raise ValueError(
            f"Cannot assign histidine protonation for source residue {residue.source_resid}: "
            "neither HD1 nor HE2 is present"
        )
    if name in {"CYS", "CYX", "CYM"}:
        if residue_index in cyx_indices:
            return "CYX"
        return "CYS" if "HG" in names else "CYM"
    if name in {"LYS", "LYN"}:
        hz = {"1HZ", "2HZ", "3HZ", "HZ1", "HZ2", "HZ3"} & names
        return "LYS" if len(hz) >= 3 else "LYN"
    if name in {"ASP", "ASH"}:
        return "ASH" if {"HD1", "HD2"} & names else "ASP"
    if name in {"GLU", "GLH"}:
        return "GLH" if {"HE1", "HE2"} & names else "GLU"
    return name


def choose_template(
    residue_index: int,
    residue: ResidueBlock,
    segments: list[list[int]],
    cyx_indices: set[int],
    templates: dict[str, list[str]],
    box: np.ndarray,
) -> str:
    variant = base_variant(residue, cyx_indices, residue_index)
    segment = next(segment for segment in segments if residue_index in segment)
    source_names = {atom.source_name for atom in residue.atoms}
    is_first = residue_index == segment[0]
    is_last = residue_index == segment[-1]
    if variant not in {"ACE", "NME"} and is_first:
        candidate = "N" + variant
        expected_n_h = 2 if variant == "PRO" else 3
        if candidate in templates and len(backbone_n_hydrogens(residue, box)) == expected_n_h:
            return candidate
    if variant not in {"ACE", "NME"} and is_last:
        candidate = "C" + variant
        if candidate in templates and "OXT" in source_names:
            return candidate
    if variant not in templates:
        raise ValueError(f"No ff19SB template {variant!r} for source residue {residue.source_resid}")
    return variant

def pdb_resname_from_template(template_name: str) -> str:
    if template_name in {"ACE", "NME"}:
        return template_name
    if template_name.startswith(("N", "C")) and template_name[1:] in VARIANT_NAMES:
        return template_name[1:]
    return template_name


def map_protein_names(
    source_names: list[str], template_names: list[str], template_name: str
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    remaining = set(template_names)
    if template_name in {"ACE", "NME"}:
        # Schrodinger/Desmond uses several equivalent methyl-cap naming
        # schemes.  Examples observed in VMD GRO exports are:
        #   ACE: CH3, HH31, HH32, HH33
        #   NME: N, CH3, H, HH31, HH32, HH33
        #   NMA: N, CA, H, HA1, HA2, HA3
        # ff19SB names the methyl hydrogens H1/H2/H3 and the NME methyl
        # carbon C.  Resolve these aliases explicitly rather than relying on
        # one historical 1HH3/2HH3/3HH3 convention.
        def cap_target(source_name: str) -> str:
            upper = source_name.upper()
            if template_name == "NME" and upper in {"CH3", "CA"}:
                return "C"
            if template_name == "ACE" and upper == "CH3":
                return "CH3"
            for pattern in (r"HH3([123])", r"([123])HH3", r"HA([123])", r"([123])HA"):
                match = re.fullmatch(pattern, upper)
                if match:
                    return "H" + match.group(1)
            return source_name

        for source in source_names:
            target = cap_target(source)
            if target not in remaining:
                raise ValueError(
                    f"Cannot map {source} to {target} in {template_name}; "
                    f"remaining={sorted(remaining)}"
                )
            mapping[source] = target
            remaining.remove(target)
    else:
        # Preserve exact names first.
        for source in source_names:
            if source in remaining:
                mapping[source] = source
                remaining.remove(source)

        # Normalize bare backbone/terminal amine hydrogen conventions.  This
        # handles H + 1H + 2H, 1H/2H/3H, H1/H2/H3 and a retained single 1H or
        # 2H that must map to the ordinary peptide atom H after ACE capping.
        generic_sources = [
            source
            for source in source_names
            if source not in mapping
            and re.fullmatch(r"(?:H|[123]H|H[123])", source.upper())
        ]
        generic_targets = [
            name for name in ("H", "H1", "H2", "H3") if name in remaining
        ]
        if generic_sources and len(generic_sources) == len(generic_targets):
            def source_h_key(name: str) -> tuple[int, int]:
                upper = name.upper()
                if upper == "H":
                    return (0, source_names.index(name))
                match = re.search(r"([123])", upper)
                return (int(match.group(1)) if match else 9, source_names.index(name))

            generic_sources.sort(key=source_h_key)
            generic_targets.sort(key=natural_key)
            for source, target in zip(generic_sources, generic_targets):
                mapping[source] = target
                remaining.remove(target)

        # Old PDB-style atom names, e.g. 1HB/2HB -> HB2/HB3.  Candidate names
        # are paired in natural order because many force-field templates use
        # suffixes 2/3 rather than 1/2.
        grouped: dict[str, list[tuple[int, str]]] = defaultdict(list)
        for source in source_names:
            if source in mapping:
                continue
            match = re.fullmatch(r"([123])(H[A-Z0-9']*)", source)
            if match:
                grouped[match.group(2)].append((int(match.group(1)), source))
        for base, source_group in grouped.items():
            candidates = sorted(
                [name for name in remaining if name.startswith(base)],
                key=natural_key,
            )
            source_group.sort()
            if len(candidates) != len(source_group):
                continue
            for (_number, source), target in zip(source_group, candidates):
                mapping[source] = target
                remaining.remove(target)

    unresolved = [name for name in source_names if name not in mapping]
    if unresolved or remaining or len(set(mapping.values())) != len(mapping):
        raise ValueError(
            f"Protein atom mapping failed for {template_name}: "
            f"unresolved={unresolved}, unused={sorted(remaining)}"
        )
    return mapping

def _template_prefix(template_name: str) -> str:
    """Return the Amber terminal prefix (N/C/empty) for an AA template."""

    if (
        len(template_name) > 1
        and template_name[0] in {"N", "C"}
        and template_name[1:] in VARIANT_NAMES
    ):
        return template_name[0]
    return ""


def _mapped_template_bonds_are_geometric(
    residue: ResidueBlock,
    mapping: dict[str, str],
    template_bonds: Sequence[tuple[str, str]],
    box: np.ndarray,
) -> bool:
    """Validate that every candidate-template bond exists geometrically."""

    source_by_name = atom_lookup(residue)
    target_to_source = {target: source for source, target in mapping.items()}
    expected_source_bonds: set[frozenset[str]] = set()
    for target_left, target_right in template_bonds:
        if target_left not in target_to_source or target_right not in target_to_source:
            return False
        source_left = target_to_source[target_left]
        source_right = target_to_source[target_right]
        expected_source_bonds.add(frozenset((source_left, source_right)))
        left = source_by_name[source_left]
        right = source_by_name[source_right]
        delta = minimum_image(right.xyz_nm - left.xyz_nm, box)
        distance = float(np.linalg.norm(delta))
        if left.element == "H" or right.element == "H":
            if not 0.065 <= distance <= 0.140:
                return False
        elif not 0.090 <= distance <= 0.190:
            return False

    # For the deliberately small ALA/GLY surrogate set, require the complete
    # inferred intrarezidue covalent graph to equal the template graph.  This
    # prevents a mere atom-name coincidence from triggering substitution.
    inferred = {
        frozenset((residue.atoms[left].source_name, residue.atoms[right].source_name))
        for left, right in infer_coordinate_bonds(residue.atoms, box)
    }
    return inferred == expected_source_bonds


def find_exact_truncated_residue_surrogate(
    residue: ResidueBlock,
    expected_template: str,
    templates: dict[str, list[str]],
    template_bonds: dict[str, list[tuple[str, str]]],
    box: np.ndarray,
) -> tuple[str, dict[str, str], str] | None:
    """Find a conservative ff19SB surrogate for a truncated side chain.

    Some prepared crystal constructs contain a residue labelled with its
    sequence identity even though unresolved side-chain atoms were deleted and
    the remaining C-beta was hydrogen-capped.  The resulting atom set can be
    exactly alanine-like (or, less commonly, glycine-like) and cannot be loaded
    as the original ff19SB residue because LEaP would add an unsampled side
    chain.

    Automatic substitution is deliberately restricted to exact ALA/GLY atom
    sets with valid covalent geometry and the same terminal state as the
    expected template.  No other side-chain reconstruction or residue guessing
    is performed.
    """

    expected_base = expected_template[1:] if _template_prefix(expected_template) else expected_template
    if expected_base in {"ALA", "GLY", "ACE", "NME"}:
        return None
    prefix = _template_prefix(expected_template)
    source_names = list(atom_lookup(residue))

    for base in ("ALA", "GLY"):
        candidate = prefix + base
        if candidate not in templates:
            continue
        if len(source_names) != len(templates[candidate]):
            continue
        try:
            mapping = map_protein_names(source_names, templates[candidate], candidate)
        except ValueError:
            continue
        if not _mapped_template_bonds_are_geometric(
            residue, mapping, template_bonds[candidate], box
        ):
            continue
        reason = (
            f"source {residue.source_resname}{residue.source_resid} is labelled "
            f"{expected_base} but its complete retained atom set and covalent "
            f"geometry exactly match {candidate}; use the sampled {base}-like "
            "chemistry rather than allowing LEaP to invent missing side-chain atoms"
        )
        return candidate, mapping, reason
    return None


def parse_mol2(path: Path) -> tuple[list[dict[str, object]], list[tuple[int, int, str]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms: list[dict[str, object]] = []
    bonds: list[tuple[int, int, str]] = []
    section = ""
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        if not line.strip():
            continue
        if section == "@<TRIPOS>ATOM":
            fields = line.split()
            if len(fields) < 9:
                raise ValueError(f"Malformed MOL2 atom line: {line}")
            atoms.append(
                {
                    "id": int(fields[0]),
                    "name": fields[1],
                    "xyz_a": np.array([float(fields[2]), float(fields[3]), float(fields[4])]),
                    "type": fields[5],
                    "resid": int(fields[6]),
                    "resname": fields[7],
                    "charge": float(fields[8]),
                    "element": infer_element(fields[1], fields[7]),
                }
            )
        elif section == "@<TRIPOS>BOND":
            fields = line.split()
            if len(fields) >= 4:
                bonds.append((int(fields[1]) - 1, int(fields[2]) - 1, fields[3]))
    if not atoms or not bonds:
        raise ValueError(f"MOL2 must contain ATOM and BOND sections: {path}")
    if [int(atom["id"]) for atom in atoms] != list(range(1, len(atoms) + 1)):
        raise ValueError("MOL2 atom IDs must be contiguous and start at 1")
    return atoms, bonds


def assess_mol2_parameterization(
    mol2_atoms: Sequence[dict[str, object]],
    requested_charge: int | None,
    preserve_charges: bool,
) -> dict[str, object]:
    """Decide whether a MOL2 can be used directly as GAFF2 input.

    Desmond/VMD commonly exports Tripos/SYBYL atom types (for example C.2,
    N.4) and force-field charges that do not sum to the molecular formal
    charge.  Those files are chemically useful for connectivity and atom
    mapping, but they are not GAFF2 parameter files.  Such inputs are routed
    through Antechamber/AM1-BCC during the ``run`` stage.
    """

    atom_types = [str(atom["type"]) for atom in mol2_atoms]
    charge_sum = float(sum(float(atom["charge"]) for atom in mol2_atoms))
    gaff2_like = all(re.fullmatch(r"[a-z][a-z0-9]*", atom_type) for atom_type in atom_types)

    if requested_charge is None:
        nearest = int(round(charge_sum))
        if not gaff2_like or abs(charge_sum - nearest) > 0.10:
            raise ValueError(
                "The ligand MOL2 is not a charge-consistent GAFF/GAFF2 file. "
                "Specify --ligand-charge so Antechamber can generate GAFF2/AM1-BCC "
                "parameters."
            )
        target_charge = nearest
    else:
        target_charge = int(requested_charge)

    reasons: list[str] = []
    if not gaff2_like:
        reasons.append("non-GAFF/SYBYL atom types")
    if abs(charge_sum - target_charge) > 0.10:
        reasons.append(
            f"input charge sum {charge_sum:.6f} differs from target {target_charge:+d}"
        )
    requires_antechamber = bool(reasons)
    if requires_antechamber and preserve_charges:
        raise ValueError(
            "--preserve-ligand-charges cannot be combined with a MOL2 that requires "
            "Antechamber reparameterization. Supply a charge-consistent GAFF2 MOL2 "
            "or omit --preserve-ligand-charges."
        )

    return {
        "input_atom_types": sorted(set(atom_types)),
        "input_charge_sum_e": charge_sum,
        "target_charge_e": target_charge,
        "gaff2_like_input_types": gaff2_like,
        "requires_antechamber": requires_antechamber,
        "method": "GAFF2/AM1-BCC via Antechamber" if requires_antechamber else "use supplied GAFF2 types/charges",
        "reasons": reasons,
    }


def infer_coordinate_bonds(
    atoms: Sequence[GroAtom], box: np.ndarray, scale: float = 1.25
) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    for left, atom_left in enumerate(atoms):
        radius_left = COVALENT_RADII_A.get(atom_left.element)
        if radius_left is None:
            raise ValueError(f"No covalent radius for ligand element {atom_left.element}")
        for right in range(left + 1, len(atoms)):
            atom_right = atoms[right]
            if atom_left.element == atom_right.element == "H":
                continue
            radius_right = COVALENT_RADII_A.get(atom_right.element)
            if radius_right is None:
                raise ValueError(f"No covalent radius for ligand element {atom_right.element}")
            delta_nm = minimum_image(atom_right.xyz_nm - atom_left.xyz_nm, box)
            distance_a = float(np.linalg.norm(delta_nm) * 10.0)
            if distance_a <= scale * (radius_left + radius_right):
                bonds.append((left, right))
    return bonds


def kabsch_rmsd(source: np.ndarray, target: np.ndarray) -> float:
    source_centered = source - source.mean(axis=0)
    target_centered = target - target.mean(axis=0)
    covariance = source_centered.T @ target_centered
    u, _singular, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    aligned = source_centered @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - target_centered) ** 2, axis=1))))


def map_ligand_to_mol2(
    source_atoms: list[GroAtom],
    mol2_atoms: list[dict[str, object]],
    mol2_bonds: list[tuple[int, int, str]],
    box: np.ndarray,
) -> tuple[list[int], list[dict[str, object]], list[tuple[int, int]], dict[str, object]]:
    if len(source_atoms) != len(mol2_atoms):
        raise ValueError(
            f"Ligand atom-count mismatch: GRO={len(source_atoms)}, MOL2={len(mol2_atoms)}"
        )
    source_bonds = infer_coordinate_bonds(source_atoms, box)
    source_graph = nx.Graph()
    target_graph = nx.Graph()
    for index, atom in enumerate(source_atoms):
        source_graph.add_node(index, element=atom.element)
    source_graph.add_edges_from(source_bonds)
    for index, atom in enumerate(mol2_atoms):
        target_graph.add_node(index, element=str(atom["element"]))
    target_graph.add_edges_from((left, right) for left, right, _order in mol2_bonds)
    if source_graph.number_of_edges() != target_graph.number_of_edges():
        raise ValueError(
            f"Ligand bond-count mismatch: inferred GRO={source_graph.number_of_edges()}, "
            f"MOL2={target_graph.number_of_edges()}"
        )

    source_heavy = [node for node in source_graph if source_graph.nodes[node]["element"] != "H"]
    target_heavy = [node for node in target_graph if target_graph.nodes[node]["element"] != "H"]
    source_h_graph = source_graph.subgraph(source_heavy).copy()
    target_h_graph = target_graph.subgraph(target_heavy).copy()
    for heavy_graph, full_graph in ((source_h_graph, source_graph), (target_h_graph, target_graph)):
        for node in heavy_graph:
            heavy_graph.nodes[node]["signature"] = (
                full_graph.nodes[node]["element"],
                heavy_graph.degree[node],
                sum(full_graph.nodes[neighbor]["element"] == "H" for neighbor in full_graph.neighbors(node)),
                tuple(sorted(full_graph.nodes[neighbor]["element"] for neighbor in heavy_graph.neighbors(node))),
            )
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        source_h_graph,
        target_h_graph,
        node_match=lambda left, right: left["signature"] == right["signature"],
    )
    source_xyz = np.array([source_atoms[index].xyz_nm * 10.0 for index in source_heavy])
    best_mapping: dict[int, int] | None = None
    best_score: tuple[float, int] | None = None
    tested = 0
    for mapping in matcher.isomorphisms_iter():
        tested += 1
        target_xyz = np.array([mol2_atoms[mapping[index]]["xyz_a"] for index in source_heavy])
        rmsd = kabsch_rmsd(source_xyz, target_xyz)
        exact_names = sum(
            source_atoms[index].source_name == str(mol2_atoms[mapping[index]]["name"])
            for index in source_heavy
        )
        score = (rmsd, -exact_names)
        if best_score is None or score < best_score:
            best_score = score
            best_mapping = dict(mapping)
        if tested >= 50000:
            break
    if best_mapping is None:
        raise ValueError("No element/connectivity-preserving ligand heavy-atom mapping was found")

    mapping = dict(best_mapping)
    for source_heavy_index, target_heavy_index in best_mapping.items():
        source_h = [
            node for node in source_graph.neighbors(source_heavy_index)
            if source_graph.nodes[node]["element"] == "H"
        ]
        target_h = [
            node for node in target_graph.neighbors(target_heavy_index)
            if target_graph.nodes[node]["element"] == "H"
        ]
        if len(source_h) != len(target_h):
            raise ValueError(
                "Hydrogen-count mismatch around mapped ligand atoms "
                f"{source_atoms[source_heavy_index].source_name} -> "
                f"{mol2_atoms[target_heavy_index]['name']}"
            )
        source_by_name = {source_atoms[index].source_name: index for index in source_h}
        target_by_name = {str(mol2_atoms[index]["name"]): index for index in target_h}
        exact = sorted(set(source_by_name) & set(target_by_name), key=natural_key)
        for name in exact:
            mapping[source_by_name[name]] = target_by_name[name]
        remaining_source = sorted(
            [index for index in source_h if index not in mapping],
            key=lambda index: natural_key(source_atoms[index].source_name),
        )
        used_target = set(mapping.values())
        remaining_target = sorted(
            [index for index in target_h if index not in used_target],
            key=lambda index: natural_key(str(mol2_atoms[index]["name"])),
        )
        mapping.update(zip(remaining_source, remaining_target))

    if set(mapping) != set(range(len(source_atoms))) or set(mapping.values()) != set(range(len(mol2_atoms))):
        raise ValueError("Ligand mapping is not a complete one-to-one mapping")
    for left, right in source_bonds:
        if not target_graph.has_edge(mapping[left], mapping[right]):
            raise ValueError("Ligand mapping does not preserve all covalent bonds")

    target_to_source_global: list[int | None] = [None] * len(mol2_atoms)
    for source_local, target_local in mapping.items():
        target_to_source_global[target_local] = source_atoms[source_local].source_index
    if any(index is None for index in target_to_source_global):
        raise RuntimeError("Internal ligand mapping error")

    source_global_to_local = {atom.source_index: index for index, atom in enumerate(source_atoms)}
    report: list[dict[str, object]] = []
    for target_index, source_global in enumerate(target_to_source_global):
        assert source_global is not None
        source_local = source_global_to_local[source_global]
        report.append(
            {
                "output_position_1based": target_index + 1,
                "source_index_1based": source_global + 1,
                "source_name": source_atoms[source_local].source_name,
                "mol2_name": mol2_atoms[target_index]["name"],
                "element": mol2_atoms[target_index]["element"],
            }
        )
    assert best_score is not None
    quality = {
        "heavy_atom_kabsch_rmsd_angstrom": float(best_score[0]),
        "exact_heavy_atom_name_matches": int(-best_score[1]),
        "heavy_atoms": len(source_heavy),
        "graph_isomorphisms_tested": tested,
        "source_bonds": source_graph.number_of_edges(),
        "mol2_bonds": target_graph.number_of_edges(),
    }
    return [int(index) for index in target_to_source_global], report, source_bonds, quality


def normalize_mol2_charges(
    source: Path,
    destination: Path,
    target_charge: int | None,
    preserve: bool,
    coordinates_angstrom: np.ndarray | None = None,
) -> dict[str, float | bool]:
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    atom_line_indices: list[int] = []
    charges: list[float] = []
    section = ""
    for index, line in enumerate(lines):
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        if section == "@<TRIPOS>ATOM" and line.strip():
            fields = line.split()
            if len(fields) >= 9:
                atom_line_indices.append(index)
                charges.append(float(fields[8]))
    if not charges:
        raise ValueError(f"No MOL2 charges found in {source}")
    replacement_coordinates: np.ndarray | None = None
    if coordinates_angstrom is not None:
        replacement_coordinates = np.asarray(coordinates_angstrom, dtype=float)
        expected_shape = (len(charges), 3)
        if replacement_coordinates.shape != expected_shape:
            raise ValueError(
                f"Replacement MOL2 coordinates must have shape {expected_shape}; "
                f"found {replacement_coordinates.shape}"
            )
        if not np.isfinite(replacement_coordinates).all():
            raise ValueError("Replacement MOL2 coordinates contain NaN or infinity")
    original = float(sum(charges))
    inferred = int(round(original)) if target_charge is None else int(target_charge)
    if target_charge is None and abs(original - inferred) > 0.10:
        raise ValueError(
            f"MOL2 charge sum {original:.6f} is not within 0.10 e of an integer. "
            "Specify --ligand-charge explicitly or repair the MOL2."
        )
    correction = 0.0 if preserve else (inferred - original) / len(charges)
    for atom_position, (line_index, old_charge) in enumerate(zip(atom_line_indices, charges)):
        fields = lines[line_index].split()
        new_charge = old_charge + correction
        fields[8] = f"{new_charge:.8f}"
        if replacement_coordinates is None:
            x, y, z = map(float, fields[2:5])
        else:
            x, y, z = map(float, replacement_coordinates[atom_position])
        lines[line_index] = (
            f"{int(fields[0]):7d} {fields[1]:<8s} "
            f"{x:12.4f} {y:12.4f} {z:12.4f} "
            f"{fields[5]:<8s} {int(fields[6]):4d} {fields[7]:<8s} {new_charge:12.8f}"
        )
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    normalized = float(sum(charge + correction for charge in charges))
    return {
        "original_charge": original,
        "target_charge": float(inferred),
        "per_atom_correction": correction,
        "output_charge": normalized,
        "charges_preserved": preserve,
        "coordinates_replaced_with_first_trajectory_frame": replacement_coordinates is not None,
    }


def finalize_antechamber_mol2(
    reference_path: Path,
    generated_path: Path,
    destination_path: Path,
    ligand_resname: str,
    target_charge: int,
    coordinate_tolerance_angstrom: float = 0.05,
) -> dict[str, object]:
    """Preserve prepared trajectory atom order after Antechamber typing.

    Antechamber normally preserves atom order and coordinates, but the MM/PBSA
    topology must match the prepared XTC exactly.  This function verifies that
    order, restores the reference atom names/coordinates, retains generated
    GAFF2 types and AM1-BCC charges, and checks the final charge.
    """

    reference_atoms, reference_bonds = parse_mol2(reference_path)
    generated_atoms, generated_bonds = parse_mol2(generated_path)
    if len(reference_atoms) != len(generated_atoms):
        raise ValueError(
            f"Antechamber changed the ligand atom count: reference={len(reference_atoms)}, "
            f"generated={len(generated_atoms)}"
        )
    reference_elements = [str(atom["element"]) for atom in reference_atoms]
    generated_elements = [str(atom["element"]) for atom in generated_atoms]
    if reference_elements != generated_elements:
        mismatches = [
            index + 1
            for index, (left, right) in enumerate(zip(reference_elements, generated_elements))
            if left != right
        ]
        raise ValueError(
            "Antechamber reordered ligand atoms or changed their elements; first "
            f"mismatching positions: {mismatches[:20]}"
        )
    coordinate_differences = np.linalg.norm(
        np.asarray([atom["xyz_a"] for atom in generated_atoms], dtype=float)
        - np.asarray([atom["xyz_a"] for atom in reference_atoms], dtype=float),
        axis=1,
    )
    maximum_coordinate_change = float(coordinate_differences.max(initial=0.0))
    if maximum_coordinate_change > coordinate_tolerance_angstrom:
        raise ValueError(
            "Antechamber altered/reordered ligand coordinates by up to "
            f"{maximum_coordinate_change:.4f} A; trajectory/topology order cannot be "
            "certified automatically."
        )
    if {
        tuple(sorted((left, right))) for left, right, _order in reference_bonds
    } != {
        tuple(sorted((left, right))) for left, right, _order in generated_bonds
    }:
        raise ValueError("Antechamber changed ligand bond connectivity")

    lines = generated_path.read_text(encoding="utf-8", errors="replace").splitlines()
    atom_line_indices: list[int] = []
    section = ""
    for line_index, line in enumerate(lines):
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        if section == "@<TRIPOS>ATOM" and line.strip():
            atom_line_indices.append(line_index)
    if len(atom_line_indices) != len(reference_atoms):
        raise ValueError("Could not identify all Antechamber MOL2 atom lines")

    generated_charges = [
        float(lines[line_index].split()[8]) for line_index in atom_line_indices
    ]
    charge_correction = (float(target_charge) - sum(generated_charges)) / len(generated_charges)
    adjusted_charges = [charge + charge_correction for charge in generated_charges]
    # Compensate the final decimal-rounding residual on the last atom so the
    # written MOL2 sums to the requested integer charge as closely as possible.
    rounded_charges = [round(charge, 8) for charge in adjusted_charges]
    rounded_charges[-1] += float(target_charge) - sum(rounded_charges)

    for atom_index, line_index in enumerate(atom_line_indices):
        reference_atom = reference_atoms[atom_index]
        generated_atom = generated_atoms[atom_index]
        x, y, z = map(float, reference_atom["xyz_a"])
        lines[line_index] = (
            f"{atom_index + 1:7d} {str(reference_atom['name']):<8s} "
            f"{x:12.4f} {y:12.4f} {z:12.4f} "
            f"{str(generated_atom['type']):<8s} {1:4d} {ligand_resname:<8s} "
            f"{rounded_charges[atom_index]:12.8f}"
        )
    destination_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    final_atoms, _final_bonds = parse_mol2(destination_path)
    final_charge = float(sum(float(atom["charge"]) for atom in final_atoms))
    if abs(final_charge - target_charge) > 0.01:
        raise ValueError(
            f"Antechamber MOL2 charge is {final_charge:.6f}; expected {target_charge:+d}"
        )
    if not all(
        re.fullmatch(r"[a-z][a-z0-9]*", str(atom["type"])) for atom in final_atoms
    ):
        raise ValueError("Antechamber output contains non-GAFF2 atom types")
    return {
        "atom_count": len(final_atoms),
        "charge_sum_e": final_charge,
        "uniform_charge_correction_e": charge_correction,
        "maximum_coordinate_change_angstrom": maximum_coordinate_change,
        "atom_order_preserved": True,
        "connectivity_preserved": True,
    }


def connected_components(n_atoms: int, bonds: Sequence[tuple[int, int]]) -> list[list[int]]:
    adjacency = [[] for _ in range(n_atoms)]
    for left, right in bonds:
        adjacency[left].append(right)
        adjacency[right].append(left)
    components: list[list[int]] = []
    seen: set[int] = set()
    for root in range(n_atoms):
        if root in seen:
            continue
        component: list[int] = []
        stack = [root]
        seen.add(root)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in adjacency[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components


def unwrap_component(
    raw: np.ndarray,
    box: np.ndarray,
    component: Sequence[int],
    adjacency: Sequence[Sequence[int]],
) -> np.ndarray:
    result = np.array(raw, copy=True)
    root = component[0]
    visited = {root}
    queue: deque[int] = deque([root])
    while queue:
        parent = queue.popleft()
        for child in adjacency[parent]:
            if child in visited:
                continue
            delta = minimum_image(raw[child] - raw[parent], box)
            result[child] = result[parent] + delta
            visited.add(child)
            queue.append(child)
    return result


def fit_coordinates(coords: np.ndarray, reference: np.ndarray, indices: np.ndarray) -> np.ndarray:
    mobile = coords[indices]
    target = reference[indices]
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (target - target_center)
    u, _singular, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return (coords - mobile_center) @ rotation + target_center


def process_frame(
    raw_coords: np.ndarray,
    box: np.ndarray,
    components: list[list[int]],
    adjacency: list[list[int]],
    main_component_index: int,
    fit_indices: np.ndarray,
    reference: np.ndarray | None,
) -> np.ndarray:
    whole = np.array(raw_coords, copy=True)
    for component in components:
        unwrapped = unwrap_component(raw_coords, box, component, adjacency)
        whole[component] = unwrapped[component]
    main_component = components[main_component_index]
    main_center = whole[main_component].mean(axis=0)
    for component_index, component in enumerate(components):
        if component_index == main_component_index:
            continue
        center = whole[component].mean(axis=0)
        nearest_delta = minimum_image(center - main_center, box)
        whole[component] += main_center + nearest_delta - center
    if reference is not None:
        whole = fit_coordinates(whole, reference, fit_indices)
    return whole


def local_cap_frame(
    origin: np.ndarray,
    direction_atom: np.ndarray,
    plane_atom: np.ndarray,
    box: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct a right-handed local frame using minimum-image vectors."""

    x = minimum_image(direction_atom - origin, box)
    x_norm = float(np.linalg.norm(x))
    if x_norm < 1.0e-8:
        raise ValueError("Cannot orient synthetic cap: zero-length direction vector")
    x /= x_norm

    plane = minimum_image(plane_atom - origin, box)
    y = plane - np.dot(plane, x) * x
    y_norm = float(np.linalg.norm(y))
    if y_norm < 1.0e-8:
        # Deterministic fallback for a pathological collinear frame.
        trial = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(trial, x))) > 0.85:
            trial = np.array([0.0, 1.0, 0.0])
        y = trial - np.dot(trial, x) * x
        y_norm = float(np.linalg.norm(y))
    y /= y_norm
    z = np.cross(x, y)
    z /= np.linalg.norm(z)
    return x, y, z


def place_synthetic_cap(
    spec: SyntheticCapSpec,
    selected_xyz: np.ndarray,
    source_to_selected: dict[int, int],
    box: np.ndarray,
) -> dict[int, np.ndarray]:
    """Generate one ACE or NME cap in the current input frame."""

    origin = selected_xyz[source_to_selected[spec.origin_source_index]]
    direction = selected_xyz[source_to_selected[spec.direction_source_index]]
    plane = selected_xyz[source_to_selected[spec.plane_source_index]]
    x, y, z = local_cap_frame(origin, direction, plane, box)
    local_coordinates = (
        NME_LOCAL_COORDS_NM if spec.template_name == "NME" else ACE_LOCAL_COORDS_NM
    )
    placed: dict[int, np.ndarray] = {}
    for atom_name, output_index in spec.output_indices.items():
        local = local_coordinates[atom_name]
        placed[output_index] = origin + local[0] * x + local[1] * y + local[2] * z
    return placed


def unitcell_lengths_angles(box: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lengths = np.linalg.norm(box, axis=1)
    angles = []
    for left, right in ((1, 2), (0, 2), (0, 1)):
        cosine = np.dot(box[left], box[right]) / (lengths[left] * lengths[right])
        angles.append(math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0)))))
    return lengths, np.asarray(angles)


def pdb_atom_name(name: str, element: str) -> str:
    if len(name) >= 4:
        return name[:4]
    if len(element) == 1:
        return f" {name:<3s}"
    return f"{name:<4s}"


def write_pdb(path: Path, atoms: Sequence[OutputAtom], coords_nm: np.ndarray, box_nm: np.ndarray) -> None:
    lengths, angles = unitcell_lengths_angles(box_nm)
    lines = [
        "CRYST1"
        f"{lengths[0] * 10:9.3f}{lengths[1] * 10:9.3f}{lengths[2] * 10:9.3f}"
        f"{angles[0]:7.2f}{angles[1]:7.2f}{angles[2]:7.2f} P 1           1\n"
    ]
    serial = 1
    previous_chain: str | None = None
    for atom, xyz_nm in zip(atoms, coords_nm):
        if previous_chain is not None and atom.chain_id != previous_chain:
            lines.append(f"TER   {serial:5d}\n")
            serial += 1
        xyz_a = xyz_nm * 10.0
        atom_field = pdb_atom_name(atom.atom_name, atom.element)
        lines.append(
            f"{atom.record_name:<6s}{serial:5d} {atom_field:4s} "
            f"{atom.residue_name:>3s} {atom.chain_id:1s}{atom.residue_number:4d}    "
            f"{xyz_a[0]:8.3f}{xyz_a[1]:8.3f}{xyz_a[2]:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {atom.element:>2s}\n"
        )
        previous_chain = atom.chain_id
        serial += 1
    lines.append(f"TER   {serial:5d}\nEND\n")
    path.write_text("".join(lines), encoding="ascii")


def write_gro(path: Path, atoms: Sequence[OutputAtom], coords_nm: np.ndarray, box_nm: np.ndarray) -> None:
    lines = ["Amber-ordered protein-ligand complex\n", f"{len(atoms):5d}\n"]
    for output_index, (atom, xyz) in enumerate(zip(atoms, coords_nm), start=1):
        lines.append(
            f"{atom.residue_number % 100000:5d}{atom.residue_name[:5]:<5s}"
            f"{atom.atom_name[:5]:>5s}{output_index % 100000:5d}"
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}\n"
        )
    fields = [
        box_nm[0, 0], box_nm[1, 1], box_nm[2, 2],
        box_nm[0, 1], box_nm[0, 2], box_nm[1, 0],
        box_nm[1, 2], box_nm[2, 0], box_nm[2, 1],
    ]
    if np.allclose(fields[3:], 0.0, atol=1e-8):
        lines.append("".join(f"{value:10.5f}" for value in fields[:3]) + "\n")
    else:
        lines.append("".join(f"{value:10.5f}" for value in fields) + "\n")
    path.write_text("".join(lines), encoding="ascii")


def write_ndx(path: Path, n_protein: int, n_total: int, ligand_name: str) -> None:
    def group(name: str, indices: Iterable[int]) -> list[str]:
        values = list(indices)
        output = [f"[ {name} ]\n"]
        for start in range(0, len(values), 15):
            output.append(" ".join(map(str, values[start:start + 15])) + "\n")
        output.append("\n")
        return output

    lines: list[str] = []
    lines.extend(group("System", range(1, n_total + 1)))
    lines.extend(group("Protein", range(1, n_protein + 1)))
    lines.extend(group(ligand_name, range(n_protein + 1, n_total + 1)))
    lines.extend(group(f"Protein_{ligand_name}", range(1, n_total + 1)))
    path.write_text("".join(lines), encoding="ascii")


def write_atom_mapping(path: Path, atoms: Sequence[OutputAtom]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "output_index_1based", "source_index_1based", "chain", "resid",
                "resname", "atom_name", "element", "record", "is_protein",
            ]
        )
        for atom in atoms:
            writer.writerow(
                [
                    atom.output_index + 1,
                    atom.source_index + 1 if atom.source_index >= 0 else "synthetic",
                    atom.chain_id,
                    atom.residue_number, atom.residue_name, atom.atom_name,
                    atom.element, atom.record_name, int(atom.is_protein),
                ]
            )


def file_sha256(path: Path) -> str:
    """Return a stable SHA-256 digest for a prepared input file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def files_are_identical(left: Path, right: Path) -> bool:
    if not left.is_file() or not right.is_file():
        return False
    if left.stat().st_size != right.stat().st_size:
        return False
    return file_sha256(left) == file_sha256(right)


def archive_stale_preparation_state(outdir: Path) -> Path | None:
    """Archive state that must never be reused by a fresh prepare command.

    A fresh prepare into an existing directory writes new base PDB/XTC files.
    Old cap-repair reports and capfixed coordinates then refer to the previous
    preparation.  Keeping those files under their live names allowed a stale
    ``validation_passed`` report to suppress cap repair.  Move the prior state
    aside before generating new coordinates.
    """

    exact_names = {
        "preparation_report.json",
        "cap_clash_repair.json",
        "cap_rotation_angles.tsv",
        "steric_preflight.json",
        "steric_preflight_before.json",
        "steric_clashes.tsv",
        "steric_clashes_before.tsv",
        "steric_preflight.stdout.log",
        "undefined_energy_report.json",
        "FINAL_RESULTS_MMPBSA.dat",
        "FINAL_RESULTS_MMPBSA.csv",
        "binding_energy_summary.json",
        "binding_energy_summary_warning.json",
    }
    patterns = (
        "*_capfix.pdb",
        "*_capfix.xtc",
        "*_stericfix*.pdb",
        "*_stericfix*.xtc",
        "*.before_capfix",
    )
    candidates = {outdir / name for name in exact_names}
    for pattern in patterns:
        candidates.update(outdir.glob(pattern))
    existing = sorted(path for path in candidates if path.exists())
    if not existing:
        return None

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    archive = outdir / "previous_run_state" / timestamp
    counter = 1
    while archive.exists():
        archive = outdir / "previous_run_state" / f"{timestamp}_{counter}"
        counter += 1
    archive.mkdir(parents=True, exist_ok=False)
    for path in existing:
        destination = archive / path.name
        shutil.move(str(path), str(destination))
    return archive


def parse_pdb_atoms(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            padded = line.rstrip("\n").ljust(80)
            atom_name = padded[12:16].strip()
            residue_name = padded[17:20].strip()
            element = padded[76:78].strip().upper() or infer_element(atom_name, residue_name)
            atoms.append(
                {
                    "name": atom_name,
                    "resname": residue_name,
                    "chain": padded[21].strip(),
                    "resid": int(padded[22:26]),
                    "element": element,
                }
            )
    return atoms


def validate_processed_trajectory(
    xyz: np.ndarray,
    output_atoms: Sequence[OutputAtom],
    output_bonds: Sequence[tuple[int, int]],
    n_protein: int,
    max_allowed_bond_nm: float = 0.30,
) -> dict[str, object]:
    bond_array = np.asarray(output_bonds, dtype=int)
    distances = np.linalg.norm(
        xyz[:, bond_array[:, 0], :] - xyz[:, bond_array[:, 1], :], axis=2
    )
    max_by_bond = distances.max(axis=0)
    overall_max = float(max_by_bond.max())
    if overall_max > max_allowed_bond_nm:
        worst = int(np.argmax(max_by_bond))
        left, right = output_bonds[worst]
        raise ValueError(
            f"Processed trajectory has a covalent bond of {overall_max:.4f} nm "
            f"between output atoms {left + 1} and {right + 1}; PBC repair failed"
        )
    protein_heavy = np.array(
        [atom.output_index for atom in output_atoms[:n_protein] if atom.element != "H"],
        dtype=int,
    )
    ligand_heavy = np.array(
        [atom.output_index for atom in output_atoms[n_protein:] if atom.element != "H"],
        dtype=int,
    )
    minimum_distances: list[float] = []
    for frame in xyz:
        delta = frame[protein_heavy][:, None, :] - frame[ligand_heavy][None, :, :]
        minimum_distances.append(float(np.sqrt(np.sum(delta * delta, axis=2)).min()))
    worst_bonds = np.argsort(max_by_bond)[::-1][:20]
    return {
        "frames": int(xyz.shape[0]),
        "atoms": int(xyz.shape[1]),
        "covalent_bonds_checked": int(len(output_bonds)),
        "maximum_covalent_bond_nm": overall_max,
        "mean_covalent_bond_nm": float(distances.mean()),
        "protein_ligand_minimum_heavy_distance_nm": {
            "minimum": min(minimum_distances),
            "maximum": max(minimum_distances),
        },
        "largest_bonds": [
            {
                "atom1_1based": int(output_bonds[index][0] + 1),
                "atom2_1based": int(output_bonds[index][1] + 1),
                "maximum_nm": float(max_by_bond[index]),
                "mean_nm": float(distances[:, index].mean()),
            }
            for index in worst_bonds
        ],
    }


def prepare_workflow(args: argparse.Namespace) -> Path:
    gro = Path(args.gro).expanduser().resolve()
    xtc = Path(args.xtc).expanduser().resolve()
    mol2 = Path(args.mol2).expanduser().resolve()
    ffxml = Path(args.ffxml).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    if args.stride < 1 or args.chunk < 1:
        raise ValueError("--stride and --chunk must be positive integers")
    for path in (gro, xtc, mol2, ffxml):
        if not path.is_file():
            raise FileNotFoundError(path)
    outdir.mkdir(parents=True, exist_ok=True)
    archived_state = archive_stale_preparation_state(outdir)
    if archived_state is not None:
        print(f"Archived stale prior-run state in {archived_state}")

    _title, all_atoms, gro_box = parse_gro(gro)
    raw_protein_blocks = contiguous_residue_blocks(all_atoms, AA_NAMES)
    (
        protein_residues,
        hydrogen_caps,
        internal_hydrogen_corrections,
    ) = separate_noncontiguous_hydrogen_caps(raw_protein_blocks, gro_box)
    if internal_hydrogen_corrections and args.internal_duplicate_h_mode == "error":
        labels = [
            f"{item.residue_key[1]}{item.residue_key[0]} "
            f"remove {item.removed_atom.source_name} "
            f"retain {item.retained_atom.source_name if item.retained_atom else 'none'}"
            for item in internal_hydrogen_corrections
        ]
        raise ValueError(
            "Peptide-connected residues contain excess backbone-N hydrogens: "
            + "; ".join(labels)
            + ". Re-run with --internal-duplicate-h-mode planar-remove "
            "to retain the planar amide H and omit the excess atom."
        )
    if not protein_residues:
        raise ValueError("No Amber-compatible protein residues were found in the GRO")

    ligand_name = args.ligand_resname.upper()
    ligand_atoms = [atom for atom in all_atoms if atom.source_resname == ligand_name]
    if not ligand_atoms:
        raise ValueError(f"No ligand residue {ligand_name!r} found in GRO")
    ligand_residue_ids = {(atom.source_resid, atom.source_resname) for atom in ligand_atoms}
    if len(ligand_residue_ids) != 1:
        raise ValueError(f"Expected exactly one {ligand_name} residue; found {ligand_residue_ids}")

    templates, template_bonds = load_ffxml(ffxml)
    for required_cap in ("ACE", "NME"):
        if required_cap not in templates:
            raise ValueError(f"ff19SB XML lacks required cap template {required_cap}")

    segments = detect_chain_segments(protein_residues, gro_box)
    hydrogen_cap_assignments = assign_hydrogen_caps_to_segment_boundaries(
        protein_residues, segments, hydrogen_caps
    )
    if hydrogen_cap_assignments and args.hydrogen_cap_mode == "error":
        raise ValueError(
            "Neutral terminal carbonyl-H/duplicate-H caps were detected. Re-run with "
            "--hydrogen-cap-mode ace-nme to replace each detected terminal "
            "hydrogen independently with a standard ff19SB NME or ACE cap."
        )
    c_cap_by_residue = {
        assignment.residue_index: assignment.cap
        for assignment in hydrogen_cap_assignments
        if assignment.replacement_template == "NME"
    }
    n_cap_by_residue = {
        assignment.residue_index: assignment.cap
        for assignment in hydrogen_cap_assignments
        if assignment.replacement_template == "ACE"
    }

    disulfide_pairs = detect_disulfides(protein_residues, gro_box)
    cyx_indices = {
        index
        for left, right, _distance in disulfide_pairs
        for index in (left, right)
    }

    output_atoms: list[OutputAtom] = []
    output_bonds: list[tuple[int, int]] = []
    protein_report: list[dict[str, object]] = []
    template_substitutions: list[dict[str, object]] = []
    residue_atom_index: dict[tuple[int, str], int] = {}
    synthetic_caps: list[SyntheticCapSpec] = []
    segment_entries: dict[int, list[dict[str, int]]] = defaultdict(list)
    cap_output_residues: dict[tuple[str, int], int] = {}
    output_residue_count = 0

    def append_synthetic_cap(
        template_name: str,
        cap: HydrogenCap,
        chain_index: int,
    ) -> dict[str, int]:
        nonlocal output_residue_count
        output_residue_count += 1
        chain_id = CHAIN_IDS[chain_index]
        local_output: dict[str, int] = {}
        for atom_name in templates[template_name]:
            output_index = len(output_atoms)
            element = infer_element(atom_name, template_name)
            output_atoms.append(
                OutputAtom(
                    output_index=output_index,
                    source_index=-1,
                    atom_name=atom_name,
                    residue_name=template_name,
                    residue_number=output_residue_count,
                    chain_id=chain_id,
                    element=element,
                    record_name="ATOM",
                    is_protein=True,
                )
            )
            local_output[atom_name] = output_index
        for atom_left, atom_right in template_bonds[template_name]:
            output_bonds.append((local_output[atom_left], local_output[atom_right]))
        synthetic_caps.append(
            SyntheticCapSpec(
                template_name=template_name,
                output_residue=output_residue_count,
                chain_id=chain_id,
                output_indices=dict(local_output),
                origin_source_index=cap.anchor_atom.source_index,
                direction_source_index=cap.source_atom.source_index,
                plane_source_index=cap.plane_atom.source_index,
                replaced_source_index=cap.source_atom.source_index,
                source_resid=cap.residue_key[0],
                source_resname=cap.residue_key[1],
            )
        )
        cap_output_residues[(template_name, cap.source_atom.source_index)] = output_residue_count
        protein_report.append(
            {
                "output_residue": output_residue_count,
                "chain": chain_id,
                "source_resid": cap.residue_key[0],
                "source_resname": cap.residue_key[1],
                "amber_template": template_name,
                "output_resname": template_name,
                "atom_count": len(templates[template_name]),
                "source_type": "synthetic_cap",
                "replaced_source_atom_index_1based": cap.source_atom.source_index + 1,
                "replaced_source_atom_name": cap.source_atom.source_name,
            }
        )
        return local_output

    def append_source_residue(residue_index: int) -> dict[str, int]:
        nonlocal output_residue_count
        residue = protein_residues[residue_index]
        output_residue_count += 1
        residue.output_resid = output_residue_count
        template = choose_template(
            residue_index, residue, segments, cyx_indices, templates, gro_box
        )
        source_by_name = atom_lookup(residue)
        surrogate_reason = ""
        try:
            name_mapping = map_protein_names(
                list(source_by_name), templates[template], template
            )
        except ValueError as original_error:
            surrogate = find_exact_truncated_residue_surrogate(
                residue, template, templates, template_bonds, gro_box
            )
            if surrogate is None:
                raise original_error
            selected_template, name_mapping, surrogate_reason = surrogate
            template_substitutions.append(
                {
                    "source_resid": residue.source_resid,
                    "source_resname": residue.source_resname,
                    "expected_template": template,
                    "selected_template": selected_template,
                    "source_atom_names": ",".join(list(source_by_name)),
                    "reason": surrogate_reason,
                }
            )
            template = selected_template

        residue.amber_template = template
        residue.output_resname = pdb_resname_from_template(template)
        target_to_source = {target: source for source, target in name_mapping.items()}
        chain_id = CHAIN_IDS[residue.chain_index]
        local_output: dict[str, int] = {}
        for target_name in templates[template]:
            source_name = target_to_source[target_name]
            source_atom = source_by_name[source_name]
            output_index = len(output_atoms)
            output_atoms.append(
                OutputAtom(
                    output_index=output_index,
                    source_index=source_atom.source_index,
                    atom_name=target_name,
                    residue_name=residue.output_resname,
                    residue_number=residue.output_resid,
                    chain_id=chain_id,
                    element=source_atom.element,
                    record_name="ATOM",
                    is_protein=True,
                )
            )
            local_output[target_name] = output_index
            residue_atom_index[(residue_index, target_name)] = output_index
        for atom_left, atom_right in template_bonds[template]:
            output_bonds.append((local_output[atom_left], local_output[atom_right]))
        protein_report.append(
            {
                "output_residue": residue.output_resid,
                "chain": chain_id,
                "source_resid": residue.source_resid,
                "source_resname": residue.source_resname,
                "amber_template": template,
                "output_resname": residue.output_resname,
                "atom_count": len(residue.atoms),
                "source_type": (
                    "source_residue_surrogate" if surrogate_reason else "source_residue"
                ),
                "replaced_source_atom_index_1based": "",
                "replaced_source_atom_name": "",
            }
        )
        return local_output

    for chain_index, segment in enumerate(segments):
        if segment[0] in n_cap_by_residue:
            segment_entries[chain_index].append(
                append_synthetic_cap("ACE", n_cap_by_residue[segment[0]], chain_index)
            )
        for residue_index in segment:
            segment_entries[chain_index].append(append_source_residue(residue_index))
        if segment[-1] in c_cap_by_residue:
            segment_entries[chain_index].append(
                append_synthetic_cap("NME", c_cap_by_residue[segment[-1]], chain_index)
            )

    # Join all consecutive residue entries in each chain. This includes the
    # newly generated ARG-C--N-NME and ACE-C--N-PHE peptide bonds.
    for chain_index in range(len(segments)):
        entries = segment_entries[chain_index]
        for left_entry, right_entry in zip(entries, entries[1:]):
            if "C" not in left_entry or "N" not in right_entry:
                raise ValueError(
                    f"Cannot join consecutive output residues in chain {CHAIN_IDS[chain_index]}: "
                    f"left atoms={sorted(left_entry)}, right atoms={sorted(right_entry)}"
                )
            output_bonds.append((left_entry["C"], right_entry["N"]))

    for left, right, _distance in disulfide_pairs:
        output_bonds.append(
            (residue_atom_index[(left, "SG")], residue_atom_index[(right, "SG")])
        )

    n_protein = len(output_atoms)
    mol2_atoms, mol2_bonds = parse_mol2(mol2)
    ligand_parameterization = assess_mol2_parameterization(
        mol2_atoms,
        requested_charge=args.ligand_charge,
        preserve_charges=args.preserve_ligand_charges,
    )
    (
        ligand_source_order,
        ligand_mapping_report,
        _source_ligand_bonds,
        ligand_mapping_quality,
    ) = map_ligand_to_mol2(ligand_atoms, mol2_atoms, mol2_bonds, gro_box)
    ligand_resnum = output_residue_count + 1
    occupied_chain_ids = {CHAIN_IDS[index] for index in range(len(segments))}
    ligand_chain = "L" if "L" not in occupied_chain_ids else CHAIN_IDS[len(segments)]
    ligand_target_to_output: dict[int, int] = {}
    for target_index, source_index in enumerate(ligand_source_order):
        source_atom = all_atoms[source_index]
        mol_atom = mol2_atoms[target_index]
        output_index = len(output_atoms)
        output_atoms.append(
            OutputAtom(
                output_index=output_index,
                source_index=source_index,
                atom_name=str(mol_atom["name"]),
                residue_name=ligand_name[:3],
                residue_number=ligand_resnum,
                chain_id=ligand_chain,
                element=str(mol_atom["element"]),
                record_name="HETATM",
                is_protein=False,
            )
        )
        ligand_target_to_output[target_index] = output_index
    for left, right, _order in mol2_bonds:
        output_bonds.append((ligand_target_to_output[left], ligand_target_to_output[right]))

    n_total = len(output_atoms)
    real_output_atoms = [atom for atom in output_atoms if atom.source_index >= 0]
    real_source_indices = [atom.source_index for atom in real_output_atoms]
    if len(set(real_source_indices)) != len(real_source_indices):
        raise RuntimeError("Real source-to-output atom mapping is not one-to-one")

    expected_source_protein = {
        atom.source_index for residue in protein_residues for atom in residue.atoms
    }
    actual_source_protein = {
        atom.source_index
        for atom in output_atoms[:n_protein]
        if atom.source_index >= 0
    }
    if actual_source_protein != expected_source_protein:
        missing = sorted(expected_source_protein - actual_source_protein)
        extra = sorted(actual_source_protein - expected_source_protein)
        raise RuntimeError(
            f"Protein source mapping mismatch: missing={missing[:20]}, extra={extra[:20]}"
        )
    if n_total != n_protein + len(ligand_atoms):
        raise RuntimeError("Complex atom count is inconsistent")

    selected_source_set = set(real_source_indices)
    for cap in synthetic_caps:
        selected_source_set.update(
            (
                cap.origin_source_index,
                cap.direction_source_index,
                cap.plane_source_index,
            )
        )
    selected_source = sorted(selected_source_set)
    source_to_selected = {
        source: index for index, source in enumerate(selected_source)
    }
    real_output_positions = np.asarray(
        [atom.output_index for atom in real_output_atoms], dtype=int
    )
    real_selected_positions = np.asarray(
        [source_to_selected[atom.source_index] for atom in real_output_atoms], dtype=int
    )

    adjacency: list[list[int]] = [[] for _ in range(n_total)]
    for left, right in output_bonds:
        adjacency[left].append(right)
        adjacency[right].append(left)
    components = connected_components(n_total, output_bonds)
    protein_set = set(range(n_protein))
    main_component_index = max(
        range(len(components)),
        key=lambda index: sum(atom in protein_set for atom in components[index]),
    )
    fit_indices = np.array(
        [
            atom.output_index
            for atom in output_atoms
            if atom.is_protein and atom.atom_name == "CA"
        ],
        dtype=int,
    )
    if len(fit_indices) < 3:
        fit_indices = np.array(
            [
                atom.output_index
                for atom in output_atoms
                if atom.is_protein and atom.element != "H"
            ],
            dtype=int,
        )

    with md.open(str(xtc)) as trajectory_file:
        n_input_frames = len(trajectory_file)
    expected_source_frame_indices = list(range(0, n_input_frames, args.stride))

    processed_xyz: list[np.ndarray] = []
    processed_time: list[np.ndarray] = []
    processed_lengths: list[np.ndarray] = []
    processed_angles: list[np.ndarray] = []
    reference: np.ndarray | None = None
    first_box: np.ndarray | None = None
    for chunk in md.iterload(
        str(xtc),
        top=str(gro),
        chunk=args.chunk,
        stride=args.stride,
        atom_indices=np.asarray(selected_source, dtype=int),
    ):
        chunk_output = np.empty((chunk.n_frames, n_total, 3), dtype=np.float32)
        for frame_index in range(chunk.n_frames):
            selected_xyz = np.asarray(chunk.xyz[frame_index], dtype=float)
            raw = np.empty((n_total, 3), dtype=float)
            raw[real_output_positions] = selected_xyz[real_selected_positions]
            box = np.asarray(chunk.unitcell_vectors[frame_index], dtype=float)
            for cap_spec in synthetic_caps:
                for output_index, coordinate in place_synthetic_cap(
                    cap_spec, selected_xyz, source_to_selected, box
                ).items():
                    raw[output_index] = coordinate
            if first_box is None:
                first_box = box.copy()
            frame = process_frame(
                raw, box, components, adjacency, main_component_index, fit_indices, reference
            )
            if reference is None:
                reference = frame.copy()
            chunk_output[frame_index] = frame.astype(np.float32)
        processed_xyz.append(chunk_output)
        processed_time.append(np.asarray(chunk.time, dtype=np.float32))
        processed_lengths.append(np.asarray(chunk.unitcell_lengths, dtype=np.float32))
        processed_angles.append(np.asarray(chunk.unitcell_angles, dtype=np.float32))

    if not processed_xyz or reference is None or first_box is None:
        raise RuntimeError("No trajectory frames were read")
    xyz = np.concatenate(processed_xyz, axis=0)
    times = np.concatenate(processed_time, axis=0)
    lengths = np.concatenate(processed_lengths, axis=0)
    angles = np.concatenate(processed_angles, axis=0)
    n_frames = int(xyz.shape[0])
    if n_frames != len(expected_source_frame_indices):
        raise RuntimeError(
            f"Stride selection produced {n_frames} frames but expected "
            f"{len(expected_source_frame_indices)} from {n_input_frames} input frames"
        )

    complex_pdb = outdir / "complex_amber_order.pdb"
    reference_pdb = outdir / "complex_reference_noH.pdb"
    receptor_pdb = outdir / "receptor_amber_order.pdb"
    ligand_pdb = outdir / f"{ligand_name}_amber_order.pdb"
    complex_gro = outdir / "complex_amber_order.gro"
    complex_xtc = outdir / f"complex_stride{args.stride}.xtc"
    index_ndx = outdir / "index.ndx"
    normalized_mol2 = outdir / f"{ligand_name}_normalized.mol2"
    parameterized_mol2 = outdir / f"{ligand_name}_gaff2.mol2"

    write_pdb(complex_pdb, output_atoms, xyz[0], first_box)
    reference_indices = [
        index for index, atom in enumerate(output_atoms) if atom.element != "H"
    ]
    write_pdb(
        reference_pdb,
        [output_atoms[index] for index in reference_indices],
        xyz[0, reference_indices],
        first_box,
    )
    write_pdb(receptor_pdb, output_atoms[:n_protein], xyz[0, :n_protein], first_box)
    write_pdb(ligand_pdb, output_atoms[n_protein:], xyz[0, n_protein:], first_box)
    write_gro(complex_gro, output_atoms, xyz[0], first_box)
    topology = md.load_pdb(str(complex_pdb)).topology
    trajectory = md.Trajectory(
        xyz=xyz,
        topology=topology,
        time=times,
        unitcell_lengths=lengths,
        unitcell_angles=angles,
    )
    trajectory.save_xtc(str(complex_xtc))
    write_ndx(index_ndx, n_protein, n_total, ligand_name)
    charge_report = normalize_mol2_charges(
        mol2,
        normalized_mol2,
        target_charge=int(ligand_parameterization["target_charge_e"]),
        # For a SYBYL/non-GAFF input, preserve its charges only as provenance;
        # Antechamber will replace them with AM1-BCC charges during run.
        preserve=(
            args.preserve_ligand_charges
            or bool(ligand_parameterization["requires_antechamber"])
        ),
        coordinates_angstrom=xyz[0, n_protein:] * 10.0,
    )
    write_atom_mapping(outdir / "atom_mapping.tsv", output_atoms)

    with (outdir / "protein_templates.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        fieldnames = [
            "output_residue",
            "chain",
            "source_resid",
            "source_resname",
            "amber_template",
            "output_resname",
            "atom_count",
            "source_type",
            "replaced_source_atom_index_1based",
            "replaced_source_atom_name",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(protein_report)

    with (outdir / "template_substitutions.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        fieldnames = [
            "source_resid",
            "source_resname",
            "expected_template",
            "selected_template",
            "source_atom_names",
            "reason",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(template_substitutions)

    with (outdir / "ligand_mapping.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        fieldnames = [
            "output_position_1based",
            "source_index_1based",
            "source_name",
            "mol2_name",
            "element",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(ligand_mapping_report)

    internal_h_rows: list[dict[str, object]] = []
    for correction in internal_hydrogen_corrections:
        internal_h_rows.append(
            {
                "previous_source_resid": correction.previous_residue_key[0],
                "previous_source_resname": correction.previous_residue_key[1],
                "source_resid": correction.residue_key[0],
                "source_resname": correction.residue_key[1],
                "removed_source_atom_index_1based": correction.removed_atom.source_index + 1,
                "removed_source_atom_name": correction.removed_atom.source_name,
                "retained_source_atom_index_1based": (
                    correction.retained_atom.source_index + 1
                    if correction.retained_atom is not None
                    else ""
                ),
                "retained_source_atom_name": (
                    correction.retained_atom.source_name
                    if correction.retained_atom is not None
                    else ""
                ),
                "c_n_distance_nm": correction.c_n_distance_nm,
                "ca_c_n_angle_deg": correction.ca_c_n_angle_deg,
                "c_n_ca_angle_deg": correction.c_n_ca_angle_deg,
                "omega_deg": correction.omega_deg,
                "retained_o_c_n_h_dihedral_deg": (
                    correction.retained_o_c_n_h_dihedral_deg
                    if correction.retained_o_c_n_h_dihedral_deg is not None
                    else ""
                ),
                "removed_o_c_n_h_dihedral_deg": correction.removed_o_c_n_h_dihedral_deg,
                "action": "omit excess internal N-H; retain planar ff19SB amide hydrogen",
            }
        )
    with (outdir / "internal_backbone_hydrogen_corrections.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        fieldnames = [
            "previous_source_resid",
            "previous_source_resname",
            "source_resid",
            "source_resname",
            "removed_source_atom_index_1based",
            "removed_source_atom_name",
            "retained_source_atom_index_1based",
            "retained_source_atom_name",
            "c_n_distance_nm",
            "ca_c_n_angle_deg",
            "c_n_ca_angle_deg",
            "omega_deg",
            "retained_o_c_n_h_dihedral_deg",
            "removed_o_c_n_h_dihedral_deg",
            "action",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(internal_h_rows)

    cap_rows: list[dict[str, object]] = []
    for assignment in hydrogen_cap_assignments:
        cap = assignment.cap
        residue = protein_residues[assignment.residue_index]
        segment = segments[assignment.segment_index]
        if assignment.boundary == "N":
            adjacent_index = (
                segments[assignment.segment_index - 1][-1]
                if assignment.segment_index > 0
                else None
            )
        else:
            adjacent_index = (
                segments[assignment.segment_index + 1][0]
                if assignment.segment_index + 1 < len(segments)
                else None
            )
        adjacent_residue = (
            protein_residues[adjacent_index] if adjacent_index is not None else None
        )
        cap_rows.append(
            {
                "segment_1based": assignment.segment_index + 1,
                "boundary": assignment.boundary,
                "source_resid": residue.source_resid,
                "source_resname": residue.source_resname,
                "cap_kind": cap.kind,
                "source_atom_index_1based": cap.source_atom.source_index + 1,
                "source_atom_name": cap.source_atom.source_name,
                "anchor_atom_name": cap.anchor_atom.source_name,
                "distance_nm": cap.distance_nm,
                "replacement_template": assignment.replacement_template,
                "replacement_output_residue": cap_output_residues[
                    (assignment.replacement_template, cap.source_atom.source_index)
                ],
                "adjacent_source_resid": (
                    adjacent_residue.source_resid if adjacent_residue is not None else ""
                ),
                "adjacent_source_resname": (
                    adjacent_residue.source_resname if adjacent_residue is not None else ""
                ),
                "adjacent_is_standard_cap": (
                    adjacent_residue.source_resname in {"ACE", "NME"}
                    if adjacent_residue is not None
                    else False
                ),
            }
        )
    with (outdir / "hydrogen_cap_replacements.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        fieldnames = [
            "segment_1based",
            "boundary",
            "source_resid",
            "source_resname",
            "cap_kind",
            "source_atom_index_1based",
            "source_atom_name",
            "anchor_atom_name",
            "distance_nm",
            "replacement_template",
            "replacement_output_residue",
            "adjacent_source_resid",
            "adjacent_source_resname",
            "adjacent_is_standard_cap",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(cap_rows)

    disulfide_rows: list[dict[str, object]] = []
    for left, right, input_distance in disulfide_pairs:
        disulfide_rows.append(
            {
                "left_output_residue": protein_residues[left].output_resid,
                "left_source_resid": protein_residues[left].source_resid,
                "right_output_residue": protein_residues[right].output_resid,
                "right_source_resid": protein_residues[right].source_resid,
                "input_sg_distance_nm": input_distance,
            }
        )
    with (outdir / "disulfides.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        fieldnames = [
            "left_output_residue",
            "left_source_resid",
            "right_output_residue",
            "right_source_resid",
            "input_sg_distance_nm",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(disulfide_rows)

    with (outdir / "covalent_bonds.tsv").open(
        "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["atom1_1based", "atom2_1based"])
        for left, right in output_bonds:
            writer.writerow([left + 1, right + 1])

    validation = validate_processed_trajectory(xyz, output_atoms, output_bonds, n_protein)
    (outdir / "trajectory_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )

    frcmod_name = f"{ligand_name}.frcmod"
    leap_lines = [
        "source leaprc.protein.ff19SB\n",
        "source leaprc.gaff2\n",
        "set default PBRadii mbondi2\n",
        f"loadamberparams {frcmod_name}\n",
        f"LIG = loadmol2 {parameterized_mol2.name}\n",
        f"REC = loadpdb {receptor_pdb.name}\n",
    ]
    for row in disulfide_rows:
        leap_lines.append(
            f"bond REC.{row['left_output_residue']}.SG "
            f"REC.{row['right_output_residue']}.SG\n"
        )
    leap_lines.extend(
        [
            "COM = combine { REC LIG }\n",
            "check REC\n",
            "check LIG\n",
            "check COM\n",
            "saveamberparm REC REC.prmtop REC.inpcrd\n",
            "saveamberparm LIG LIG.prmtop LIG.inpcrd\n",
            "saveamberparm COM COM.prmtop COM.inpcrd\n",
            "savepdb COM COM_tleap.pdb\n",
            "quit\n",
        ]
    )
    (outdir / "tleap.in").write_text("".join(leap_lines), encoding="ascii")

    if args.solvent_model == "pb":
        solvent_block = "&pb\n  istrng=0.150,\n  radiopt=0,\n/\n"
    else:
        solvent_block = "&gb\n  igb=5,\n  saltcon=0.150,\n/\n"
    mmpbsa_text = (
        f"Protein-ligand MM/{args.solvent_model.upper()}SA on every "
        f"{args.stride}th input frame\n"
        "&general\n"
        f"  sys_name=\"{ligand_name}_ff19SB_GAFF2_{args.solvent_model.upper()}\",\n"
        "  startframe=1,\n"
        f"  endframe={n_frames},\n"
        "  interval=1,\n"
        "  PBRadii=3,\n"
        "  verbose=1,\n"
        "/\n"
        + solvent_block
    )
    (outdir / "mmpbsa.in").write_text(mmpbsa_text, encoding="ascii")

    shutil.copy2(ffxml, outdir / "protein.ff19SB.xml")
    report = {
        "miniapp_version": __version__,
        "input": {
            "gro": str(gro),
            "xtc": str(xtc),
            "mol2": str(mol2),
            "ffxml": str(ffxml),
        },
        "ligand_resname": ligand_name,
        "solvent_model": args.solvent_model,
        "hydrogen_cap_mode": args.hydrogen_cap_mode,
        "internal_duplicate_h_mode": args.internal_duplicate_h_mode,
        "protein_atoms": n_protein,
        "protein_source_atoms_retained": len(expected_source_protein),
        "protein_source_cap_hydrogens_replaced": len(synthetic_caps),
        "protein_source_internal_hydrogens_removed": len(internal_hydrogen_corrections),
        "synthetic_cap_atoms_added": sum(
            len(cap.output_indices) for cap in synthetic_caps
        ),
        "ligand_atoms": len(ligand_atoms),
        "complex_atoms": n_total,
        "source_protein_residues": len(protein_residues),
        "protein_residues": output_residue_count,
        "protein_chains": len(segments),
        "input_trajectory_frames": n_input_frames,
        "trajectory_frames_written": n_frames,
        "stride": args.stride,
        "source_frame_indices": expected_source_frame_indices,
        "hydrogen_cap_replacements": cap_rows,
        "internal_backbone_hydrogen_corrections": internal_h_rows,
        "template_substitutions": template_substitutions,
        "disulfides": disulfide_rows,
        "ligand_charge": charge_report,
        "ligand_parameterization": ligand_parameterization,
        "ligand_mapping_quality": ligand_mapping_quality,
        "connected_components": [len(component) for component in components],
        "trajectory_validation": validation,
        "groups": {
            "System": 0,
            "Protein": 1,
            ligand_name: 2,
            f"Protein_{ligand_name}": 3,
        },
        "files": {
            "complex_pdb": complex_pdb.name,
            "reference_pdb": reference_pdb.name,
            "receptor_pdb": receptor_pdb.name,
            "ligand_pdb": ligand_pdb.name,
            "complex_gro": complex_gro.name,
            "complex_xtc": complex_xtc.name,
            "index": index_ndx.name,
            "normalized_mol2": normalized_mol2.name,
            "parameterized_mol2": parameterized_mol2.name,
            "tleap_input": "tleap.in",
            "mmpbsa_input": "mmpbsa.in",
            "internal_backbone_hydrogen_corrections": (
                "internal_backbone_hydrogen_corrections.tsv"
            ),
            "template_substitutions": "template_substitutions.tsv",
        },
    }
    (outdir / "preparation_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    print(f"\nPreparation completed: {outdir}")
    return outdir

def run_logged(command: Sequence[str], cwd: Path, log_path: Path) -> None:
    printable = " ".join(subprocess.list2cmdline([part]) for part in command)
    print(f"\n$ {printable}")
    with log_path.open("wt", encoding="utf-8") as log_handle:
        process = subprocess.Popen(
            list(command),
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            log_handle.write(line)
        return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, list(command))


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(
            f"Required executable {name!r} is not in PATH. Activate the "
            "AmberTools/GROMACS/gmx_MMPBSA environment first."
        )
    return path


def command_version(command: Sequence[str]) -> str:
    try:
        completed = subprocess.run(
            list(command),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
        return completed.stdout.strip()[:4000]
    except Exception as exc:  # pragma: no cover - diagnostic only
        return f"version query failed: {exc}"


def topology_atom_signature(atom: object) -> tuple[str, int]:
    name = str(getattr(atom, "name", "")).strip()
    atomic_number = int(getattr(atom, "atomic_number", 0) or 0)
    if atomic_number == 0:
        element = infer_element(name)
        atomic_number = ATOMIC_NUMBERS.get(element, 0)
    return name, atomic_number



def _pdb_residue_blocks(path: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Return contiguous PDB atoms and residues with zero-based atom indices."""

    atoms = parse_pdb_atoms(path)
    residues: list[dict[str, object]] = []
    current_key: tuple[str, int, str] | None = None
    current: dict[str, object] | None = None
    for atom_index, atom in enumerate(atoms):
        key = (str(atom["chain"]), int(atom["resid"]), str(atom["resname"]))
        if key != current_key:
            current = {
                "chain": key[0],
                "resid": key[1],
                "resname": key[2],
                "atom_indices": [],
            }
            residues.append(current)
            current_key = key
        assert current is not None
        current["atom_indices"].append(atom_index)
    return atoms, residues


def _kabsch_reference_to_mobile(
    reference: np.ndarray, mobile: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return centers and row-vector rotation mapping reference to mobile."""

    if reference.shape != mobile.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("Kabsch inputs must have matching (N,3) shapes")
    if reference.shape[0] < 3:
        raise ValueError("At least three common atoms are required to place a missing atom")
    reference_center = reference.mean(axis=0)
    mobile_center = mobile.mean(axis=0)
    covariance = (reference - reference_center).T @ (mobile - mobile_center)
    u, _singular, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return reference_center, mobile_center, rotation


def _backup_once(path: Path, suffix: str = ".before_topology_reconcile") -> None:
    backup = path.with_name(path.name + suffix)
    if path.is_file() and not backup.exists():
        shutil.copy2(path, backup)


def reconcile_tleap_terminal_additions(
    workdir: Path,
    report: dict[str, object],
    rec: object,
    lig: object,
    com: object,
) -> dict[str, object] | None:
    """Insert terminal atoms added by tleap into the prepared PDB/XTC.

    The preparation stage normally produces an atom-exact ff19SB PDB.  Some
    Desmond exports leave a segment end in an internal-like form (for example,
    a carbonyl O without OXT).  When LEaP reads such a chain end it applies the
    standard Amber terminal patch and adds one or more terminal atoms.  Rather
    than accepting a topology/trajectory count mismatch, this routine maps the
    LEaP topology residue by residue and propagates only recognized terminal
    additions into every prepared trajectory frame.

    Accepted automatic additions are deliberately narrow: N-terminal H/H1/H2/
    H3 and C-terminal OXT.  Any side-chain or internal-residue addition remains
    a hard error because it indicates an incomplete or incorrectly mapped
    source structure.
    """

    files = report.get("files")
    if not isinstance(files, dict):
        raise ValueError("Invalid preparation report: files must be an object")
    complex_pdb = workdir / str(files["complex_pdb"])
    receptor_pdb = workdir / str(files["receptor_pdb"])
    ligand_pdb = workdir / str(files["ligand_pdb"])
    complex_gro = workdir / str(files["complex_gro"])
    complex_xtc = workdir / str(files["complex_xtc"])
    reference_pdb = workdir / str(files.get("reference_pdb", files["complex_pdb"]))
    index_path = workdir / str(files.get("index", "index.ndx"))
    atom_mapping_path = workdir / "atom_mapping.tsv"
    old_source_indices: list[int] = []
    if atom_mapping_path.is_file():
        with atom_mapping_path.open(
            "rt", encoding="utf-8", errors="replace", newline=""
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                value = str(row.get("source_index_1based", "")).strip()
                old_source_indices.append(int(value) - 1 if value.isdigit() else -1)

    source_complex_atoms, source_complex_residues = _pdb_residue_blocks(complex_pdb)
    source_rec_atoms, source_rec_residues = _pdb_residue_blocks(receptor_pdb)
    source_lig_atoms, source_lig_residues = _pdb_residue_blocks(ligand_pdb)
    if len(source_rec_atoms) + len(source_lig_atoms) != len(source_complex_atoms):
        raise RuntimeError("Prepared complex is not receptor followed by ligand")
    if len(old_source_indices) != len(source_complex_atoms):
        old_source_indices = [-1] * len(source_complex_atoms)
    if len(source_rec_residues) != len(rec.residues):
        raise RuntimeError(
            "Cannot reconcile LEaP additions: receptor residue counts differ "
            f"({len(source_rec_residues)} prepared vs {len(rec.residues)} topology)"
        )
    if len(source_lig_residues) != len(lig.residues):
        raise RuntimeError(
            "Cannot reconcile LEaP additions: ligand residue counts differ "
            f"({len(source_lig_residues)} prepared vs {len(lig.residues)} topology)"
        )

    # Determine first/last residue positions within each source PDB chain.
    first_by_chain: dict[str, int] = {}
    last_by_chain: dict[str, int] = {}
    for residue_index, residue in enumerate(source_rec_residues):
        chain = str(residue["chain"])
        first_by_chain.setdefault(chain, residue_index)
        last_by_chain[chain] = residue_index

    rec_reference_a = np.asarray(getattr(rec, "coordinates", None), dtype=float)
    if rec_reference_a.shape != (len(rec.atoms), 3):
        raise RuntimeError("REC.inpcrd coordinates are unavailable for topology reconciliation")
    rec_reference_nm = rec_reference_a * 0.1

    topology_order: list[tuple[str, int, int | None]] = []
    additions: list[dict[str, object]] = []
    residue_maps: list[dict[str, object]] = []
    rec_top_atom_offset = 0

    for residue_index, (source_residue, top_residue) in enumerate(
        zip(source_rec_residues, rec.residues)
    ):
        source_indices = list(source_residue["atom_indices"])
        source_by_name: dict[str, int] = {}
        for source_index in source_indices:
            name = str(source_rec_atoms[source_index]["name"])
            if name in source_by_name:
                raise RuntimeError(
                    f"Duplicate atom name {name!r} in prepared residue "
                    f"{source_residue['resid']} {source_residue['resname']}"
                )
            source_by_name[name] = source_index
        top_atoms = list(top_residue.atoms)
        top_names = [str(atom.name).strip() for atom in top_atoms]
        if len(set(top_names)) != len(top_names):
            raise RuntimeError(
                f"Duplicate atom names in Amber topology residue {residue_index + 1}"
            )

        # Exact mapping first, followed by terminal-amine hydrogen aliases.
        # LEaP can rename a prepared backbone H to H1 while adding H2/H3.
        top_to_source: dict[str, int] = {
            name: source_by_name[name] for name in top_names if name in source_by_name
        }
        used_source = set(top_to_source.values())
        chain = str(source_residue["chain"])
        is_first = residue_index == first_by_chain[chain]
        is_last = residue_index == last_by_chain[chain]
        if is_first:
            generic_source = [
                (name, index)
                for name, index in source_by_name.items()
                if index not in used_source
                and re.fullmatch(r"(?:H|[123]H|H[123])", name.upper())
            ]
            generic_target = [
                name
                for name in top_names
                if name not in top_to_source and name in {"H", "H1", "H2", "H3"}
            ]
            generic_source.sort(key=lambda item: natural_key(item[0]))
            generic_target.sort(key=natural_key)
            for (source_name, source_index), target_name in zip(
                generic_source, generic_target
            ):
                top_to_source[target_name] = source_index
                used_source.add(source_index)

        extra_source = sorted(
            str(source_rec_atoms[index]["name"])
            for index in source_indices
            if index not in used_source
        )
        missing_top = [name for name in top_names if name not in top_to_source]
        if extra_source:
            raise RuntimeError(
                "LEaP topology omitted prepared receptor atoms in residue "
                f"{source_residue['resid']} {source_residue['resname']}: {extra_source}"
            )
        if missing_top:
            allowed: set[str] = set()
            if is_first:
                allowed.update({"H", "H1", "H2", "H3"})
            if is_last:
                allowed.add("OXT")
            invalid = [name for name in missing_top if name not in allowed]
            if invalid:
                raise RuntimeError(
                    "LEaP added nonterminal or unsupported atoms in residue "
                    f"{source_residue['resid']} {source_residue['resname']}: {invalid}. "
                    "Only terminal H/H1/H2/H3 and OXT can be reconciled automatically."
                )

        common_names = [name for name in top_names if name in top_to_source]
        # Prefer heavy atoms for a stable local frame, then include hydrogens
        # only if fewer than three heavy atoms are available.
        heavy_common = [
            name
            for name in common_names
            if int(getattr(top_atoms[top_names.index(name)], "atomic_number", 0) or 0) != 1
        ]
        # Terminal additions are determined by the local peptide frame, not by
        # a potentially mobile side chain.  Prefer backbone heavy atoms and
        # only fall back to all common heavy atoms when necessary.
        preferred_backbone = [
            name for name in ("N", "CA", "C", "O", "CB") if name in heavy_common
        ]
        if len(preferred_backbone) >= 3:
            alignment_names = preferred_backbone
        elif len(heavy_common) >= 3:
            alignment_names = heavy_common
        else:
            alignment_names = common_names
        if missing_top and len(alignment_names) < 3:
            raise RuntimeError(
                f"Not enough common atoms to place {missing_top} in residue "
                f"{source_residue['resid']} {source_residue['resname']}"
            )

        top_index_by_name = {
            str(atom.name).strip(): rec_top_atom_offset + local_index
            for local_index, atom in enumerate(top_atoms)
        }
        residue_maps.append(
            {
                "source_residue": source_residue,
                "top_to_source": top_to_source,
                "top_names": top_names,
                "top_index_by_name": top_index_by_name,
                "alignment_names": alignment_names,
                "missing_names": missing_top,
            }
        )
        for name in top_names:
            if name in top_to_source:
                topology_order.append(("source", top_to_source[name], None))
            else:
                top_global = top_index_by_name[name]
                topology_order.append(("generated", residue_index, top_global))
                additions.append(
                    {
                        "residue_order_1based": residue_index + 1,
                        "chain": source_residue["chain"],
                        "resid": source_residue["resid"],
                        "resname": source_residue["resname"],
                        "atom_name": name,
                        "topology_atom_index_1based": top_global + 1,
                    }
                )
        rec_top_atom_offset += len(top_atoms)

    # Ligand is expected to be atom-exact.  Map by names only to guard against
    # a harmless residue renumbering while rejecting any topology alteration.
    ligand_order: list[int] = []
    lig_top_offset = 0
    for source_residue, top_residue in zip(source_lig_residues, lig.residues):
        source_indices = list(source_residue["atom_indices"])
        source_by_name = {
            str(source_lig_atoms[index]["name"]): index for index in source_indices
        }
        top_names = [str(atom.name).strip() for atom in top_residue.atoms]
        if set(source_by_name) != set(top_names) or len(source_by_name) != len(top_names):
            raise RuntimeError("Ligand topology differs from the prepared ligand atom set")
        ligand_order.extend(source_by_name[name] for name in top_names)
        lig_top_offset += len(top_names)

    if not additions:
        return None

    trajectory = md.load(str(complex_xtc), top=str(complex_pdb))
    if trajectory.n_atoms != len(source_complex_atoms):
        raise RuntimeError(
            f"Prepared XTC has {trajectory.n_atoms} atoms but PDB has {len(source_complex_atoms)}"
        )
    source_xyz = np.asarray(trajectory.xyz, dtype=np.float64)
    n_rec_top = len(rec.atoms)
    n_lig_top = len(lig.atoms)
    new_xyz = np.empty((trajectory.n_frames, n_rec_top + n_lig_top, 3), dtype=np.float64)

    # Copy/generate receptor atoms in Amber topology order.
    for frame_index in range(trajectory.n_frames):
        frame = source_xyz[frame_index]
        transforms: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for residue_index, residue_map in enumerate(residue_maps):
            missing_names = list(residue_map["missing_names"])
            if not missing_names:
                continue
            names = list(residue_map["alignment_names"])
            top_to_source = dict(residue_map["top_to_source"])
            top_index_by_name = dict(residue_map["top_index_by_name"])
            reference = np.asarray(
                [rec_reference_nm[top_index_by_name[name]] for name in names],
                dtype=float,
            )
            mobile = np.asarray(
                [frame[top_to_source[name]] for name in names],
                dtype=float,
            )
            transforms[residue_index] = _kabsch_reference_to_mobile(reference, mobile)

        for output_index, (kind, value, top_global) in enumerate(topology_order):
            if kind == "source":
                new_xyz[frame_index, output_index] = frame[value]
            else:
                residue_index = value
                assert top_global is not None
                reference_center, mobile_center, rotation = transforms[residue_index]
                new_xyz[frame_index, output_index] = (
                    (rec_reference_nm[top_global] - reference_center) @ rotation
                    + mobile_center
                )

        ligand_source_offset = len(source_rec_atoms)
        for local_output, source_lig_index in enumerate(ligand_order):
            new_xyz[frame_index, n_rec_top + local_output] = frame[
                ligand_source_offset + source_lig_index
            ]

    # Build output labels in topology order while preserving source chain and
    # residue identifiers.  Terminal Amber variant labels are intentionally
    # written in ordinary PDB form, as in the original preparation stage.
    output_atoms: list[OutputAtom] = []
    output_index = 0
    receptor_top_global = 0
    for residue_index, (source_residue, top_residue) in enumerate(
        zip(source_rec_residues, rec.residues)
    ):
        for atom in top_residue.atoms:
            atomic_number = int(getattr(atom, "atomic_number", 0) or 0)
            element = next(
                (symbol for symbol, number in ATOMIC_NUMBERS.items() if number == atomic_number),
                infer_element(str(atom.name), str(source_residue["resname"])),
            )
            kind, old_output_index, _top_index = topology_order[receptor_top_global]
            original_source_index = (
                old_source_indices[old_output_index]
                if kind == "source"
                else -1
            )
            output_atoms.append(
                OutputAtom(
                    output_index=output_index,
                    source_index=original_source_index,
                    atom_name=str(atom.name).strip(),
                    residue_name=str(source_residue["resname"]),
                    residue_number=int(source_residue["resid"]),
                    chain_id=str(source_residue["chain"]),
                    element=element,
                    record_name="ATOM",
                    is_protein=True,
                )
            )
            output_index += 1
            receptor_top_global += 1
    ligand_chain = str(source_lig_residues[0]["chain"]) if source_lig_residues else "L"
    ligand_resid = int(source_lig_residues[0]["resid"]) if source_lig_residues else len(rec.residues) + 1
    ligand_resname = str(report["ligand_resname"])
    for ligand_top_index, atom in enumerate(lig.atoms):
        atomic_number = int(getattr(atom, "atomic_number", 0) or 0)
        element = next(
            (symbol for symbol, number in ATOMIC_NUMBERS.items() if number == atomic_number),
            infer_element(str(atom.name), ligand_resname),
        )
        old_output_index = len(source_rec_atoms) + ligand_order[ligand_top_index]
        output_atoms.append(
            OutputAtom(
                output_index=output_index,
                source_index=old_source_indices[old_output_index],
                atom_name=str(atom.name).strip(),
                residue_name=ligand_resname,
                residue_number=ligand_resid,
                chain_id=ligand_chain,
                element=element,
                record_name="HETATM",
                is_protein=False,
            )
        )
        output_index += 1

    if trajectory.unitcell_vectors is not None:
        first_box = np.asarray(trajectory.unitcell_vectors[0], dtype=float)
    else:
        first_box = np.eye(3, dtype=float) * 10.0

    for path in (
        complex_pdb,
        receptor_pdb,
        ligand_pdb,
        reference_pdb,
        complex_gro,
        complex_xtc,
        index_path,
        atom_mapping_path,
    ):
        _backup_once(path)

    write_pdb(complex_pdb, output_atoms, new_xyz[0], first_box)
    write_pdb(receptor_pdb, output_atoms[:n_rec_top], new_xyz[0, :n_rec_top], first_box)
    write_pdb(ligand_pdb, output_atoms[n_rec_top:], new_xyz[0, n_rec_top:], first_box)
    heavy_indices = [
        index for index, atom in enumerate(output_atoms) if atom.element.upper() != "H"
    ]
    write_pdb(
        reference_pdb,
        [output_atoms[index] for index in heavy_indices],
        new_xyz[0, heavy_indices],
        first_box,
    )
    write_gro(complex_gro, output_atoms, new_xyz[0], first_box)
    rebuilt_topology = md.load_pdb(str(complex_pdb)).topology
    rebuilt = md.Trajectory(
        xyz=new_xyz.astype(np.float32),
        topology=rebuilt_topology,
        time=np.asarray(trajectory.time, dtype=np.float32),
        unitcell_lengths=(
            None
            if trajectory.unitcell_lengths is None
            else np.asarray(trajectory.unitcell_lengths, dtype=np.float32)
        ),
        unitcell_angles=(
            None
            if trajectory.unitcell_angles is None
            else np.asarray(trajectory.unitcell_angles, dtype=np.float32)
        ),
    )
    rebuilt.save_xtc(str(complex_xtc))
    write_ndx(index_path, n_rec_top, n_rec_top + n_lig_top, ligand_resname)
    write_atom_mapping(atom_mapping_path, output_atoms)

    report["protein_atoms"] = n_rec_top
    report["complex_atoms"] = n_rec_top + n_lig_top
    trajectory_validation = report.get("trajectory_validation")
    if isinstance(trajectory_validation, dict):
        trajectory_validation["atoms"] = n_rec_top + n_lig_top
    reconciliation = {
        "applied": True,
        "prepared_receptor_atoms_before": len(source_rec_atoms),
        "amber_receptor_atoms_after": n_rec_top,
        "prepared_complex_atoms_before": len(source_complex_atoms),
        "amber_complex_atoms_after": n_rec_top + n_lig_top,
        "additions": additions,
        "coordinate_method": (
            "Per-frame residue-local Kabsch transform of tleap-generated "
            "terminal coordinates"
        ),
        "accepted_atom_classes": ["N-terminal H/H1/H2/H3", "C-terminal OXT"],
    }
    report["topology_reconciliation"] = reconciliation
    (workdir / "topology_reconciliation.json").write_text(
        json.dumps(reconciliation, indent=2) + "\n", encoding="utf-8"
    )
    (workdir / "preparation_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return reconciliation

def validate_and_convert_topologies(workdir: Path, report: dict[str, object]) -> dict[str, object]:
    try:
        import parmed as pmd
    except ImportError as exc:
        raise RuntimeError(
            "ParmEd is required for topology validation/conversion. Install it "
            "in the active gmx_MMPBSA environment."
        ) from exc

    expected_complex = int(report["complex_atoms"])
    expected_receptor = int(report["protein_atoms"])
    expected_ligand = int(report["ligand_atoms"])
    files = report["files"]
    assert isinstance(files, dict)
    complex_pdb = workdir / str(files["complex_pdb"])
    receptor_pdb = workdir / str(files["receptor_pdb"])
    ligand_pdb = workdir / str(files["ligand_pdb"])

    pdb_complex = parse_pdb_atoms(complex_pdb)
    pdb_receptor = parse_pdb_atoms(receptor_pdb)
    pdb_ligand = parse_pdb_atoms(ligand_pdb)
    if [len(pdb_complex), len(pdb_receptor), len(pdb_ligand)] != [
        expected_complex, expected_receptor, expected_ligand
    ]:
        raise RuntimeError(
            "Prepared PDB counts changed unexpectedly: "
            f"complex/receptor/ligand={len(pdb_complex)}/{len(pdb_receptor)}/{len(pdb_ligand)}"
        )

    topologies: dict[str, object] = {}
    for prefix in ("COM", "REC", "LIG"):
        prmtop = workdir / f"{prefix}.prmtop"
        inpcrd = workdir / f"{prefix}.inpcrd"
        if not prmtop.is_file() or not inpcrd.is_file():
            raise FileNotFoundError(f"Missing {prmtop.name} or {inpcrd.name} after tleap")
        topologies[prefix] = pmd.load_file(str(prmtop), xyz=str(inpcrd))

    com = topologies["COM"]
    rec = topologies["REC"]
    lig = topologies["LIG"]
    counts = {"COM": len(com.atoms), "REC": len(rec.atoms), "LIG": len(lig.atoms)}
    expected_counts = {"COM": expected_complex, "REC": expected_receptor, "LIG": expected_ligand}
    reconciliation = None
    if counts != expected_counts:
        delta_rec = counts["REC"] - expected_counts["REC"]
        delta_com = counts["COM"] - expected_counts["COM"]
        if (
            delta_rec > 0
            and delta_com == delta_rec
            and counts["LIG"] == expected_counts["LIG"]
        ):
            reconciliation = reconcile_tleap_terminal_additions(
                workdir, report, rec, lig, com
            )
            if reconciliation is not None:
                print(
                    "WARNING: tleap added recognized terminal atoms; the prepared "
                    "PDB/XTC were reconciled to the Amber topology. See "
                    "topology_reconciliation.json."
                )
                expected_receptor = int(report["protein_atoms"])
                expected_complex = int(report["complex_atoms"])
                expected_counts = {
                    "COM": expected_complex,
                    "REC": expected_receptor,
                    "LIG": expected_ligand,
                }
                pdb_complex = parse_pdb_atoms(complex_pdb)
                pdb_receptor = parse_pdb_atoms(receptor_pdb)
                pdb_ligand = parse_pdb_atoms(ligand_pdb)
        if counts != expected_counts:
            raise RuntimeError(
                f"Amber topology atom counts {counts} do not match {expected_counts}. "
                "Automatic reconciliation is limited to terminal H/H1/H2/H3 and OXT additions."
            )

    def compare_pdb_to_topology(label: str, pdb_atoms: list[dict[str, object]], topology: object) -> None:
        mismatches: list[str] = []
        for index, (pdb_atom, top_atom) in enumerate(zip(pdb_atoms, topology.atoms), start=1):
            top_name, top_atomic_number = topology_atom_signature(top_atom)
            pdb_atomic_number = ATOMIC_NUMBERS.get(str(pdb_atom["element"]).upper(), 0)
            if str(pdb_atom["name"]) != top_name or (
                pdb_atomic_number and top_atomic_number and pdb_atomic_number != top_atomic_number
            ):
                mismatches.append(
                    f"{index}: PDB={pdb_atom['name']}/{pdb_atom['element']} "
                    f"topology={top_name}/{top_atomic_number}"
                )
                if len(mismatches) >= 20:
                    break
        if mismatches:
            raise RuntimeError(
                f"{label} PDB/topology atom order differs. First mismatches:\n" + "\n".join(mismatches)
            )

    compare_pdb_to_topology("Complex", pdb_complex, com)
    compare_pdb_to_topology("Receptor", pdb_receptor, rec)
    compare_pdb_to_topology("Ligand", pdb_ligand, lig)

    com_signatures = [topology_atom_signature(atom) for atom in com.atoms]
    rec_lig_signatures = [topology_atom_signature(atom) for atom in rec.atoms] + [
        topology_atom_signature(atom) for atom in lig.atoms
    ]
    if com_signatures != rec_lig_signatures:
        raise RuntimeError("COM.prmtop atom order is not REC.prmtop followed by LIG.prmtop")

    for index, (com_atom, separate_atom) in enumerate(
        zip(com.atoms, list(rec.atoms) + list(lig.atoms)), start=1
    ):
        if abs(float(com_atom.charge) - float(separate_atom.charge)) > 1e-7:
            raise RuntimeError(f"Charge mismatch between combined and separate topology at atom {index}")

    amber_reference = workdir / "amber_reference"
    amber_reference.mkdir(exist_ok=True)
    for prefix in ("COM", "REC", "LIG"):
        shutil.copy2(workdir / f"{prefix}.prmtop", amber_reference / f"{prefix}.prmtop")
        shutil.copy2(workdir / f"{prefix}.inpcrd", amber_reference / f"{prefix}.inpcrd")

    converted_counts: dict[str, int] = {}
    for prefix, topology in topologies.items():
        top_path = workdir / f"{prefix}.top"
        gro_path = workdir / f"{prefix}_from_amber.gro"
        topology.save(str(top_path), overwrite=True)
        topology.save(str(gro_path), overwrite=True)
        roundtrip = pmd.load_file(str(top_path), xyz=str(gro_path))
        converted_counts[prefix] = len(roundtrip.atoms)
        if converted_counts[prefix] != counts[prefix]:
            raise RuntimeError(
                f"ParmEd GROMACS round-trip changed {prefix} atom count: "
                f"{counts[prefix]} -> {converted_counts[prefix]}"
            )
        original_signatures = [topology_atom_signature(atom) for atom in topology.atoms]
        roundtrip_signatures = [topology_atom_signature(atom) for atom in roundtrip.atoms]
        if original_signatures != roundtrip_signatures:
            raise RuntimeError(f"ParmEd GROMACS round-trip changed {prefix} atom order")

    charges = {
        "COM": float(sum(atom.charge for atom in com.atoms)),
        "REC": float(sum(atom.charge for atom in rec.atoms)),
        "LIG": float(sum(atom.charge for atom in lig.atoms)),
    }
    if abs(charges["COM"] - charges["REC"] - charges["LIG"]) > 1e-6:
        raise RuntimeError("Combined topology charge is not receptor plus ligand charge")

    validation = {
        "atom_counts": counts,
        "gromacs_roundtrip_atom_counts": converted_counts,
        "charges_e": charges,
        "pdb_atom_order_matches": True,
        "combined_equals_receptor_plus_ligand": True,
        "gromacs_topology_roundtrip_matches": True,
        "parmed_version": getattr(pmd, "__version__", "unknown"),
        "topology_reconciliation": reconciliation,
    }
    (workdir / "topology_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    return validation


def parse_binding_energy(result_file: Path) -> dict[str, object]:
    """Parse the endpoint binding-energy statistics from gmx_MMPBSA output.

    gmx_MMPBSA output formats differ by release.  Older releases commonly
    printed ``DELTA TOTAL`` while newer releases print the Unicode label
    ``ΔTOTAL``.  This parser accepts both forms, optional separators, Unicode
    minus signs, and either five-column or reduced statistics tables.
    """

    text = result_file.read_text(encoding="utf-8", errors="replace")
    normalized_text = unicodedata.normalize("NFKC", text).replace("\u2212", "-")
    frame_match = re.search(
        r"Calculations performed using\s+(\d+)\s+complex frames",
        normalized_text,
        flags=re.I,
    )

    number_pattern = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    candidates: list[tuple[int, list[float], str]] = []
    for line_number, raw_line in enumerate(normalized_text.splitlines(), start=1):
        # Strip table-border characters without changing the scientific label.
        line = raw_line.strip().lstrip("|+").strip()
        label_match = re.match(
            r"^(?:(?:Δ|DELTA)\s*[_ -]*TOTAL)\b(.*)$",
            line,
            flags=re.I,
        )
        if not label_match:
            continue
        values = [
            float(token)
            for token in re.findall(number_pattern, label_match.group(1))
        ]
        if values:
            candidates.append((line_number, values, raw_line))

    if not candidates:
        # Diagnostic context makes a genuinely unsupported output format easy
        # to identify without misreporting a successful gmx_MMPBSA run.
        total_lines = [
            line for line in normalized_text.splitlines()
            if "TOTAL" in line.upper() or "Δ" in line
        ]
        preview = "\n".join(total_lines[-12:])
        raise ValueError(
            f"No DELTA TOTAL/ΔTOTAL statistics row found in {result_file}."
            + (f" Candidate TOTAL lines:\n{preview}" if preview else "")
        )

    line_number, values, source_line = candidates[-1]
    if len(values) >= 5:
        average, sd_prop, standard_deviation, sem_prop, standard_error = values[:5]
    elif len(values) == 3:
        average, standard_deviation, standard_error = values
        sd_prop = None
        sem_prop = None
    elif len(values) == 2:
        average, standard_deviation = values
        standard_error = None
        sd_prop = None
        sem_prop = None
    else:
        raise ValueError(
            f"DELTA TOTAL/ΔTOTAL row on line {line_number} of {result_file} "
            f"contains only {len(values)} numeric value: {source_line!r}"
        )

    # Determine the solvent model from the section preceding the selected row.
    prefix = "\n".join(normalized_text.splitlines()[:line_number])
    model = "unknown"
    pb_positions = [
        match.start()
        for match in re.finditer(
            r"POISSON\s+BOLTZMANN|\bPB\s+CALCULATION\b",
            prefix,
            flags=re.I,
        )
    ]
    gb_positions = [
        match.start()
        for match in re.finditer(
            r"GENERALIZED\s+BORN|\bGB\s+CALCULATION\b",
            prefix,
            flags=re.I,
        )
    ]
    if pb_positions or gb_positions:
        model = "PB" if max(pb_positions or [-1]) > max(gb_positions or [-1]) else "GB"

    stats: dict[str, float | None] = {
        "average": average,
        "standard_deviation_propagated": sd_prop,
        "standard_deviation": standard_deviation,
        "standard_error_propagated": sem_prop,
        "standard_error": standard_error,
    }
    summary: dict[str, object] = {
        "result_file": str(result_file.resolve()),
        "solvent_model": model,
        "matched_label": "ΔTOTAL" if "Δ" in source_line else "DELTA TOTAL",
        "matched_line_number": line_number,
        "delta_total_kcal_per_mol": stats,
    }
    if frame_match:
        summary["frames"] = int(frame_match.group(1))
    return summary



def run_topology_aware_steric_preflight(
    workdir: Path, report: dict[str, object]
) -> dict[str, object]:
    """Validate actual Amber LJ geometry before the expensive endpoint run.

    The preflight scans every frame with COM.prmtop exclusions and Lennard-
    Jones parameters.  It may alter only residues explicitly recorded as
    synthetic ACE/NME caps; source protein and ligand coordinates are never
    changed automatically.
    """

    script = Path(__file__).with_name("steric_preflight.py")
    if not script.is_file():
        raise FileNotFoundError(f"Bundled steric preflight script was not found: {script}")
    repair_script = Path(__file__).with_name("repair_synthetic_caps.py")
    if not repair_script.is_file():
        raise FileNotFoundError(f"Bundled cap-repair script was not found: {repair_script}")
    command = [
        sys.executable,
        str(script),
        "--workdir", str(workdir),
        "--repair-script", str(repair_script),
    ]
    run_logged(command, workdir, workdir / "steric_preflight.stdout.log")
    report_path = workdir / "steric_preflight.json"
    if not report_path.is_file():
        raise RuntimeError("Steric preflight completed without steric_preflight.json")
    result = json.loads(report_path.read_text(encoding="utf-8", errors="replace"))
    if not bool(result.get("validation_passed", False)):
        raise RuntimeError(
            "Topology-aware steric preflight failed. See steric_preflight.json "
            "and steric_clashes.tsv for the exact frame and atom pair."
        )
    report["topology_aware_steric_preflight"] = result
    (workdir / "preparation_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return result


def diagnose_undefined_energy_outputs(workdir: Path) -> dict[str, object]:
    """Record overflow/NaN fields retained in thread-specific Sander outputs."""

    candidates: list[Path] = []
    for pattern in (
        "_GMXMMPBSA_*_pb.mdout.*",
        "_GMXMMPBSA_*_gb.mdout.*",
        "_GMXMMPBSA_*mdout*",
    ):
        candidates.extend(workdir.glob(pattern))
    files = sorted(set(path for path in candidates if path.is_file()))
    findings: list[dict[str, object]] = []
    for path in files:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        energy_record = 0
        for line_number, line in enumerate(lines, start=1):
            if re.search(r"\bNSTEP\s+ENERGY\s+RMS", line):
                energy_record += 1
            undefined = "*" in line or bool(re.search(r"\bnan\b|\binf(?:inity)?\b", line, flags=re.I))
            if not undefined:
                continue
            if not re.search(
                r"BOND|ANGLE|DIHED|VDWAALS|EEL|EGB|EPB|ENPOLAR|EDISPER|1-4|ENERGY",
                line,
                flags=re.I,
            ):
                continue
            start = max(0, line_number - 7)
            stop = min(len(lines), line_number + 4)
            findings.append(
                {
                    "file": path.name,
                    "line_number": line_number,
                    "local_energy_record_1based": energy_record or None,
                    "line": line,
                    "context": lines[start:stop],
                }
            )
    report = {
        "workdir": str(workdir),
        "sander_output_files_examined": [path.name for path in files],
        "undefined_energy_fields": findings,
        "steric_preflight_report": (
            str((workdir / "steric_preflight.json").resolve())
            if (workdir / "steric_preflight.json").is_file()
            else None
        ),
    }
    (workdir / "undefined_energy_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report

def run_workflow(args: argparse.Namespace) -> Path:
    workdir = Path(args.workdir).expanduser().resolve()
    report_path = workdir / "preparation_report.json"
    if not report_path.is_file():
        raise FileNotFoundError(f"Prepared run manifest not found: {report_path}")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    ligand_name = str(report["ligand_resname"])
    files = report["files"]
    if not isinstance(files, dict):
        raise ValueError("Invalid preparation_report.json: files must be an object")

    parmchk2 = require_executable("parmchk2")
    tleap = require_executable("tleap")
    gmx = require_executable("gmx")
    gmx_mmpbsa = require_executable("gmx_MMPBSA")
    mpirun = require_executable("mpirun") if args.np > 1 else None
    parameterization = report.get("ligand_parameterization", {})
    requires_antechamber = bool(parameterization.get("requires_antechamber", False))
    antechamber = require_executable("antechamber") if requires_antechamber else None

    try:
        import parmed as pmd
    except ImportError as exc:
        raise RuntimeError("ParmEd is not importable in the active Python environment") from exc

    tool_manifest = {
        "python": sys.version,
        "parmed": getattr(pmd, "__version__", "unknown"),
        "gmx": command_version([gmx, "--version"]),
        "gmx_MMPBSA": command_version([gmx_mmpbsa, "-v"]),
        "parmchk2": parmchk2,
        "antechamber": antechamber,
        "tleap": tleap,
        "mpi_processes": args.np,
    }
    (workdir / "tool_versions.json").write_text(
        json.dumps(tool_manifest, indent=2) + "\n", encoding="utf-8"
    )

    normalized_mol2 = workdir / str(files["normalized_mol2"])
    parameterized_mol2 = workdir / str(
        files.get("parameterized_mol2", files["normalized_mol2"])
    )
    ligand_parameterization_runtime: dict[str, object]
    if requires_antechamber:
        assert antechamber is not None
        target_charge = int(parameterization["target_charge_e"])
        raw_antechamber = workdir / f"{ligand_name}_antechamber_raw.mol2"
        for stale_ligand_file in (raw_antechamber, parameterized_mol2):
            if stale_ligand_file.exists():
                stale_ligand_file.unlink()
        run_logged(
            [
                antechamber,
                "-i", normalized_mol2.name,
                "-fi", "mol2",
                "-o", raw_antechamber.name,
                "-fo", "mol2",
                "-at", "gaff2",
                "-c", "bcc",
                "-nc", str(target_charge),
                "-rn", ligand_name,
                "-s", "2",
            ],
            workdir,
            workdir / "antechamber.stdout.log",
        )
        if not raw_antechamber.is_file() or raw_antechamber.stat().st_size == 0:
            raise RuntimeError("Antechamber did not produce a nonempty MOL2 file")
        ligand_parameterization_runtime = finalize_antechamber_mol2(
            normalized_mol2,
            raw_antechamber,
            parameterized_mol2,
            ligand_name,
            target_charge,
        )
        ligand_parameterization_runtime["method"] = "GAFF2/AM1-BCC via Antechamber"
    else:
        if normalized_mol2.resolve() != parameterized_mol2.resolve():
            shutil.copy2(normalized_mol2, parameterized_mol2)
        parameterized_atoms, _parameterized_bonds = parse_mol2(parameterized_mol2)
        ligand_parameterization_runtime = {
            "method": "supplied GAFF2 atom types/charges",
            "atom_count": len(parameterized_atoms),
            "charge_sum_e": float(
                sum(float(atom["charge"]) for atom in parameterized_atoms)
            ),
            "atom_order_preserved": True,
            "connectivity_preserved": True,
        }
    (workdir / "ligand_parameterization_runtime.json").write_text(
        json.dumps(ligand_parameterization_runtime, indent=2) + "\n",
        encoding="utf-8",
    )

    frcmod = workdir / f"{ligand_name}.frcmod"
    run_logged(
        [
            parmchk2,
            "-i", parameterized_mol2.name,
            "-f", "mol2",
            "-o", frcmod.name,
            "-s", "2",
        ],
        workdir,
        workdir / "parmchk2.stdout.log",
    )
    if not frcmod.is_file() or frcmod.stat().st_size == 0:
        raise RuntimeError("parmchk2 did not produce a nonempty frcmod file")
    frcmod_text = frcmod.read_text(encoding="utf-8", errors="replace")
    if "ATTN" in frcmod_text.upper():
        print(f"WARNING: {frcmod.name} contains ATTN entries; inspect parameter analogies.")

    run_logged([tleap, "-f", "tleap.in"], workdir, workdir / "tleap.stdout.log")
    leap_text = (workdir / "tleap.stdout.log").read_text(encoding="utf-8", errors="replace")
    leap_log = workdir / "leap.log"
    if leap_log.is_file():
        leap_text += "\n" + leap_log.read_text(encoding="utf-8", errors="replace")
    if re.search(r"Fatal Error|Errors\s*=\s*[1-9]", leap_text, flags=re.I):
        raise RuntimeError("tleap reported a fatal error; inspect tleap.stdout.log and leap.log")

    topology_validation = validate_and_convert_topologies(workdir, report)
    print("\nTopology validation passed:")
    print(json.dumps(topology_validation, indent=2))

    steric_validation = run_topology_aware_steric_preflight(workdir, report)
    print("\nTopology-aware steric preflight passed:")
    print(
        json.dumps(
            {
                "frames": steric_validation.get("frames"),
                "global_minimum_interacting_pair": steric_validation.get(
                    "global_minimum_interacting_pair"
                ),
                "global_maximum_positive_pair_lj": steric_validation.get(
                    "global_maximum_positive_pair_lj"
                ),
                "repair_history": steric_validation.get("repair_history", []),
            },
            indent=2,
        )
    )

    complex_pdb = str(files["complex_pdb"])
    reference_pdb = str(files.get("reference_pdb", files["complex_pdb"]))
    complex_xtc = str(files["complex_xtc"])

    # gmx_MMPBSA intentionally retains intermediates after an error. Remove
    # only generated scratch/results from earlier attempts so a corrected run
    # cannot consume stale files. Prepared inputs and validated topologies are
    # left intact.
    stale_patterns = (
        "_GMXMMPBSA_*",
        "COM_traj_*.xtc",
        "REC_traj_*.xtc",
        "LIG_traj_*.xtc",
        "FINAL_RESULTS_MMPBSA.dat",
        "FINAL_RESULTS_MMPBSA.csv",
        "binding_energy_summary.json",
        "gmx_MMPBSA.log",
    )
    for pattern in stale_patterns:
        for stale_path in workdir.glob(pattern):
            if stale_path.is_dir():
                shutil.rmtree(stale_path)
            else:
                stale_path.unlink()

    # This is a single-trajectory (ST) calculation. Supply only the complex
    # topology and let gmx_MMPBSA derive receptor and ligand topologies from
    # the two complex index groups. Passing -rp/-lp without corresponding
    # separate receptor/ligand structures and indexes makes gmx_MMPBSA 1.6.5
    # enter its multiple-trajectory topology branch with empty REC/LIG index
    # maps and fail in cleantop().
    command = [
        gmx_mmpbsa,
        "-O",
        "-i", "mmpbsa.in",
        "-cs", complex_pdb,
        "-cr", reference_pdb,
        "-ci", "index.ndx",
        "-cg", "1", "2",
        "-ct", complex_xtc,
        "-cp", "COM.top",
        "-o", "FINAL_RESULTS_MMPBSA.dat",
        "-eo", "FINAL_RESULTS_MMPBSA.csv",
        "-nogui",
    ]
    if args.np > 1:
        assert mpirun is not None
        command = [mpirun, "-np", str(args.np)] + command
    try:
        run_logged(command, workdir, workdir / "gmx_MMPBSA.stdout.log")
    except subprocess.CalledProcessError as exc:
        diagnosis = diagnose_undefined_energy_outputs(workdir)
        if diagnosis.get("undefined_energy_fields"):
            raise RuntimeError(
                "gmx_MMPBSA failed because Sander wrote one or more undefined/"
                "overflowed energy fields. Exact mdout context was saved to "
                "undefined_energy_report.json; topology-aware geometry details "
                "are in steric_preflight.json and steric_clashes.tsv."
            ) from exc
        raise

    result_file = workdir / "FINAL_RESULTS_MMPBSA.dat"
    csv_file = workdir / "FINAL_RESULTS_MMPBSA.csv"
    if not result_file.is_file() or not csv_file.is_file():
        raise RuntimeError("gmx_MMPBSA completed without both requested result files")
    try:
        summary = parse_binding_energy(result_file)
    except Exception as exc:
        warning = {
            "result_file": str(result_file.resolve()),
            "calculation_completed": True,
            "summary_parsing_completed": False,
            "warning": str(exc),
        }
        (workdir / "binding_energy_summary_warning.json").write_text(
            json.dumps(warning, indent=2) + "\n", encoding="utf-8"
        )
        print(
            "WARNING: gmx_MMPBSA completed and produced its result files, "
            f"but automatic summary parsing failed: {exc}",
            file=sys.stderr,
        )
        return result_file

    summary["trajectory_stride"] = int(report["stride"])
    summary["prepared_frames"] = int(report["trajectory_frames_written"])
    summary["force_field"] = "Amber ff19SB + GAFF2"
    summary["entropy_included"] = False
    (workdir / "binding_energy_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print("\nBinding-energy summary:")
    print(json.dumps(summary, indent=2))
    return result_file


def apply_synthetic_cap_clash_repair(outdir: Path) -> dict[str, object]:
    """Generate and select validated synthetic-cap coordinates.

    The repaired PDB/XTC are kept as explicit ``*_capfix`` files and the run
    manifest is updated to reference them.  The base prepared coordinates are
    left immutable.  A prior report is reused only when its output files still
    exist and are demonstrably associated with the current preparation.
    """

    repair_script = Path(__file__).with_name("repair_synthetic_caps.py")
    if not repair_script.is_file():
        raise FileNotFoundError(
            f"Bundled synthetic-cap repair script was not found: {repair_script}"
        )
    preparation_report_path = outdir / "preparation_report.json"
    if not preparation_report_path.is_file():
        raise FileNotFoundError(preparation_report_path)
    preparation_report = json.loads(
        preparation_report_path.read_text(encoding="utf-8", errors="replace")
    )
    prepared_files = preparation_report.get("files", {})
    if not isinstance(prepared_files, dict):
        raise ValueError("Invalid preparation_report.json: files must be an object")

    input_pdb_name = str(prepared_files.get("complex_pdb", "complex_amber_order.pdb"))
    input_xtc_name = str(prepared_files.get("complex_xtc", "complex_stride10.xtc"))
    input_pdb_path = outdir / Path(input_pdb_name).name
    input_xtc_path = outdir / Path(input_xtc_name).name

    existing_report_path = outdir / "cap_clash_repair.json"
    if existing_report_path.is_file():
        existing_report = json.loads(
            existing_report_path.read_text(encoding="utf-8", errors="replace")
        )
        if bool(existing_report.get("validation_passed", False)):
            caps = existing_report.get("caps", [])
            if not caps:
                preparation_report["synthetic_cap_clash_repair"] = existing_report
                preparation_report["miniapp_version"] = __version__
                preparation_report_path.write_text(
                    json.dumps(preparation_report, indent=2) + "\n", encoding="utf-8"
                )
                return existing_report

            reported_pdb = outdir / Path(str(existing_report.get("output_pdb", ""))).name
            reported_xtc = outdir / Path(str(existing_report.get("output_xtc", ""))).name
            manifest_uses_reported = (
                Path(input_pdb_name).name == reported_pdb.name
                and Path(input_xtc_name).name == reported_xtc.name
            )
            if manifest_uses_reported and reported_pdb.is_file() and reported_xtc.is_file():
                preparation_report["synthetic_cap_clash_repair"] = existing_report
                preparation_report["miniapp_version"] = __version__
                preparation_report_path.write_text(
                    json.dumps(preparation_report, indent=2) + "\n", encoding="utf-8"
                )
                return existing_report

            # Compatibility recovery for folders produced by versions 4-9:
            # those versions replaced the base files and retained both a
            # pre-capfix backup and a separate validated capfix output.  A
            # later prepare could overwrite the base files while leaving the
            # old report.  Reuse is safe only when the current base files are
            # byte-identical to the recorded backups.
            backup_pdb = input_pdb_path.with_suffix(input_pdb_path.suffix + ".before_capfix")
            backup_xtc = input_xtc_path.with_suffix(input_xtc_path.suffix + ".before_capfix")
            if (
                reported_pdb.is_file()
                and reported_xtc.is_file()
                and files_are_identical(input_pdb_path, backup_pdb)
                and files_are_identical(input_xtc_path, backup_xtc)
            ):
                prepared_files["complex_pdb"] = reported_pdb.name
                prepared_files["complex_xtc"] = reported_xtc.name
                existing_report["recovered_from_overwritten_active_files"] = True
                existing_report["selected_pdb_sha256"] = file_sha256(reported_pdb)
                existing_report["selected_xtc_sha256"] = file_sha256(reported_xtc)
                preparation_report["files"] = prepared_files
                preparation_report["synthetic_cap_clash_repair"] = existing_report
                preparation_report["miniapp_version"] = __version__
                preparation_report_path.write_text(
                    json.dumps(preparation_report, indent=2) + "\n", encoding="utf-8"
                )
                existing_report_path.write_text(
                    json.dumps(existing_report, indent=2) + "\n", encoding="utf-8"
                )
                print(
                    "Recovered validated capfixed PDB/XTC after the base files "
                    "were overwritten by a later prepare run."
                )
                return existing_report

        # The report cannot be tied safely to the current input files.  Keep
        # it for audit and generate a fresh repair.
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        stale = existing_report_path.with_name(
            existing_report_path.name + f".stale.{timestamp}"
        )
        existing_report_path.rename(stale)

    # If the manifest currently points to a prior capfix name but no reusable
    # report exists, return to the immutable base preparation before creating
    # a new capfix file.
    base_pdb = outdir / "complex_amber_order.pdb"
    if base_pdb.is_file() and "_capfix" in input_pdb_path.stem:
        input_pdb_path = base_pdb
        input_pdb_name = base_pdb.name
    base_xtc_candidates = sorted(
        path for path in outdir.glob("complex_stride*.xtc")
        if "_capfix" not in path.stem and "_stericfix" not in path.stem
    )
    if base_xtc_candidates and "_capfix" in input_xtc_path.stem:
        input_xtc_path = base_xtc_candidates[0]
        input_xtc_name = input_xtc_path.name

    if not input_pdb_path.is_file():
        raise FileNotFoundError(input_pdb_path)
    if not input_xtc_path.is_file():
        raise FileNotFoundError(input_xtc_path)

    output_pdb_name = f"{input_pdb_path.stem}_capfix{input_pdb_path.suffix}"
    output_xtc_name = f"{input_xtc_path.stem}_capfix{input_xtc_path.suffix}"
    command = [
        sys.executable,
        str(repair_script),
        "--workdir", str(outdir),
        "--pdb", input_pdb_name,
        "--xtc", input_xtc_name,
        "--output-pdb", output_pdb_name,
        "--output-xtc", output_xtc_name,
        "--allow-no-caps",
    ]
    subprocess.run(command, check=True)
    repair_report_path = outdir / "cap_clash_repair.json"
    if not repair_report_path.is_file():
        raise RuntimeError("Synthetic-cap repair did not produce cap_clash_repair.json")
    repair_report = json.loads(repair_report_path.read_text(encoding="utf-8"))
    if not bool(repair_report.get("validation_passed", False)):
        raise RuntimeError("Synthetic-cap clash validation failed")

    if repair_report.get("caps"):
        output_pdb = outdir / output_pdb_name
        output_xtc = outdir / output_xtc_name
        if not output_pdb.is_file() or not output_xtc.is_file():
            raise RuntimeError("Validated cap-repair outputs are missing")
        prepared_files["complex_pdb"] = output_pdb.name
        prepared_files["complex_xtc"] = output_xtc.name
        repair_report["selected_pdb_sha256"] = file_sha256(output_pdb)
        repair_report["selected_xtc_sha256"] = file_sha256(output_xtc)
    else:
        prepared_files["complex_pdb"] = input_pdb_path.name
        prepared_files["complex_xtc"] = input_xtc_path.name

    preparation_report["files"] = prepared_files
    preparation_report["synthetic_cap_clash_repair"] = repair_report
    preparation_report["miniapp_version"] = __version__
    preparation_report_path.write_text(
        json.dumps(preparation_report, indent=2) + "\n", encoding="utf-8"
    )
    repair_report_path.write_text(
        json.dumps(repair_report, indent=2) + "\n", encoding="utf-8"
    )
    return repair_report


def add_prepare_arguments(parser: argparse.ArgumentParser) -> None:
    default_ffxml = Path(__file__).with_name("protein.ff19SB.xml")
    parser.add_argument("--gro", required=True, help="Full-system GRO matching the XTC atom order")
    parser.add_argument("--xtc", required=True, help="Full-system XTC converted from Desmond")
    parser.add_argument("--mol2", required=True, help="Explicit-H Antechamber-compatible ligand MOL2")
    parser.add_argument("--ligand-resname", required=True, help="Ligand residue name in the GRO, e.g. QNB")
    parser.add_argument("--outdir", required=True, help="Output/work directory")
    parser.add_argument("--ffxml", default=str(default_ffxml), help="Bundled ff19SB OpenMM XML")
    parser.add_argument("--stride", type=int, default=10, help="Retain every Nth input frame (default: 10)")
    parser.add_argument("--chunk", type=int, default=20, help="MDTraj input chunk size (default: 20)")
    parser.add_argument(
        "--ligand-charge",
        type=int,
        default=None,
        help="Integer ligand net charge. Default: nearest integer to the MOL2 sum.",
    )
    parser.add_argument(
        "--preserve-ligand-charges",
        action="store_true",
        help="Do not apply the tiny uniform correction that makes MOL2 charges sum to an integer.",
    )
    parser.add_argument(
        "--solvent-model",
        choices=("pb", "gb"),
        default="pb",
        help="Implicit solvent model used by gmx_MMPBSA (default: pb).",
    )
    parser.add_argument(
        "--hydrogen-cap-mode",
        choices=("ace-nme", "error"),
        default="ace-nme",
        help=(
            "How to handle non-contiguous Desmond terminal carbonyl-H/duplicate-H caps: "
            "replace each independently with an ff19SB NME/ACE cap (default) "
            "or stop with an error."
        ),
    )
    parser.add_argument(
        "--internal-duplicate-h-mode",
        choices=("planar-remove", "error"),
        default="planar-remove",
        help=(
            "How to handle a peptide-connected internal residue carrying one "
            "excess backbone-N hydrogen: retain the planar amide H and omit "
            "the excess source H from ff19SB rescoring (default), or stop."
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare and run ff19SB/GAFF2 gmx_MMPBSA rescoring of Desmond GRO/XTC trajectories."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="Prepare dry, ordered, stride-selected inputs")
    add_prepare_arguments(prepare_parser)

    run_parser = subparsers.add_parser("run", help="Build topologies and run gmx_MMPBSA")
    run_parser.add_argument("--workdir", required=True, help="Folder produced by the prepare command")
    run_parser.add_argument("--np", type=int, default=1, help="MPI processes for gmx_MMPBSA (default: 1)")

    all_parser = subparsers.add_parser("all", help="Prepare, build, validate, and run in one command")
    add_prepare_arguments(all_parser)
    all_parser.add_argument("--np", type=int, default=1, help="MPI processes for gmx_MMPBSA (default: 1)")

    summarize_parser = subparsers.add_parser("summarize", help="Parse DELTA TOTAL from a completed result")
    summarize_parser.add_argument("--result", required=True, help="FINAL_RESULTS_MMPBSA.dat")
    summarize_parser.add_argument("--output", default=None, help="Optional summary JSON path")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "prepare":
        outdir = prepare_workflow(args)
        apply_synthetic_cap_clash_repair(outdir)
    elif args.command == "run":
        if args.np < 1:
            raise ValueError("--np must be at least 1")
        workdir = Path(args.workdir).expanduser().resolve()
        apply_synthetic_cap_clash_repair(workdir)
        run_workflow(args)
    elif args.command == "all":
        if args.np < 1:
            raise ValueError("--np must be at least 1")
        outdir = prepare_workflow(args)
        apply_synthetic_cap_clash_repair(outdir)
        run_args = argparse.Namespace(workdir=str(outdir), np=args.np)
        run_workflow(run_args)
    elif args.command == "summarize":
        result = Path(args.result).expanduser().resolve()
        summary = parse_binding_energy(result)
        output = Path(args.output).expanduser().resolve() if args.output else result.with_name("binding_energy_summary.json")
        output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(summary, indent=2))
    else:  # pragma: no cover
        parser.error(f"Unknown command {args.command}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
