#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "app" / "desmond_gmxmmpbsa.py"
FFXML = ROOT / "app" / "protein.ff19SB.xml"
spec = importlib.util.spec_from_file_location("workflow_engine_v9", ENGINE)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

templates, template_bonds = module.load_ffxml(FFXML)

# Explicit ACE/NME naming variants observed in Desmond/VMD GRO exports.
assert module.map_protein_names(
    ["CH3", "C", "O", "HH31", "HH32", "HH33"],
    templates["ACE"],
    "ACE",
) == {
    "CH3": "CH3", "C": "C", "O": "O",
    "HH31": "H1", "HH32": "H2", "HH33": "H3",
}
assert module.map_protein_names(
    ["N", "CH3", "H", "HH31", "HH32", "HH33"],
    templates["NME"],
    "NME",
)["HH33"] == "H3"
assert module.map_protein_names(
    ["N", "CA", "H", "HA1", "HA2", "HA3"],
    templates["NME"],
    "NME",
) == {
    "N": "N", "CA": "C", "H": "H",
    "HA1": "H1", "HA2": "H2", "HA3": "H3",
}
assert module.map_protein_names(
    ["N", "CA", "H", "1HA", "2HA", "3HA"],
    templates["NME"],
    "NME",
)["3HA"] == "H3"

# NMA is a source-residue alias for the ff19SB NME cap.
nma = module.ResidueBlock(457, "NMA", [])
assert module.base_variant(nma, set(), 0) == "NME"

# Synthetic reproduction of the supplied terminal ARG427: after removal of
# its carbonyl-H cap, the complete sampled graph is exactly ALA-like.
raw = {
    "N": ("N", [1.457, 1.893, -2.151]),
    "CA": ("C", [1.461, 2.026, -2.092]),
    "C": ("C", [1.438, 2.134, -2.197]),
    "O": ("O", [1.325, 2.181, -2.216]),
    "CB": ("C", [1.592, 2.049, -2.015]),
    "H": ("H", [1.544, 1.849, -2.176]),
    "HA": ("H", [1.379, 2.032, -2.020]),
    "HB1": ("H", [1.592, 2.149, -1.972]),
    "HB2": ("H", [1.601, 1.976, -1.935]),
    "HB3": ("H", [1.677, 2.039, -2.083]),
}
atoms = [
    module.GroAtom(index, 427, "ARG", name, index + 1, np.array(xyz), element)
    for index, (name, (element, xyz)) in enumerate(raw.items())
]
residue = module.ResidueBlock(427, "ARG", atoms)
box = np.diag([10.0, 10.0, 10.0])
surrogate = module.find_exact_truncated_residue_surrogate(
    residue, "ARG", templates, template_bonds, box
)
assert surrogate is not None
selected, mapping, reason = surrogate
assert selected == "ALA"
assert set(mapping.values()) == set(templates["ALA"])
assert "invent missing side-chain atoms" in reason

# A single coordinate perturbation that breaks the exact covalent graph must
# prevent the automatic surrogate.
broken_atoms = [
    module.GroAtom(
        atom.source_index,
        atom.source_resid,
        atom.source_resname,
        atom.source_name,
        atom.atom_number,
        atom.xyz_nm.copy(),
        atom.element,
    )
    for atom in atoms
]
next(atom for atom in broken_atoms if atom.source_name == "HB1").xyz_nm += np.array([0.5, 0.0, 0.0])
broken = module.ResidueBlock(427, "ARG", broken_atoms)
assert module.find_exact_truncated_residue_surrogate(
    broken, "ARG", templates, template_bonds, box
) is None

print("Engine-9 atom-mapping regression tests passed")
