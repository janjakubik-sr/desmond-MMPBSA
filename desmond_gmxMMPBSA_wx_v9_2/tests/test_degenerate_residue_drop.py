#!/usr/bin/env python3
"""Regression test: degenerate chain-break debris residues.

A Desmond export can contain a residue key whose backbone carbonyl carbon
carries two C-bound capping hydrogens and whose CA atom is missing entirely
(upstream editing deleted the backbone and capped the dangling bonds).  The
engine must exclude such a residue from the ff19SB model in the default
``drop`` mode and still stop in the strict ``error`` mode.
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

spec = importlib.util.spec_from_file_location("workflow_engine_degenerate_v912", ENGINE)
assert spec and spec.loader
engine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = engine
spec.loader.exec_module(engine)
assert engine.__version__ == "9.1.4"


def atom_line(resid: int, resname: str, name: str, number: int, x: float, y: float, z: float) -> str:
    return f"{resid:5d}{resname:<5.5s}{name:<5.5s}{number:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def write_gro(path: Path, lines: list[str]) -> None:
    path.write_text(
        "regression debris test\n"
        f"{len(lines)}\n" + "".join(line + "\n" for line in lines)
        + "   5.00000   5.00000   5.00000\n",
        encoding="utf-8",
    )


# Residues 1-2 form a peptide with a neutral carbonyl-H (HXT) cap on residue 2.
# Residue 5 LEU is degenerate chain-break debris: its carbonyl C carries both
# 1HXT and 2HXT and the residue has no CA atom.  Residue 6 starts a new segment
# as a neutral NH2 terminus.
ATOMS = [
    atom_line(1, "VAL", "N", 1, 1.000, 1.000, 1.000),
    atom_line(1, "VAL", "CA", 2, 1.149, 1.000, 1.000),
    atom_line(1, "VAL", "C", 3, 1.202, 1.138, 1.000),
    atom_line(1, "VAL", "O", 4, 1.140, 1.240, 1.000),
    atom_line(1, "VAL", "CB", 5, 1.200, 0.950, 1.120),
    atom_line(2, "ALA", "N", 6, 1.150, 1.260, 1.000),
    atom_line(2, "ALA", "CA", 7, 1.240, 1.380, 1.000),
    atom_line(2, "ALA", "C", 8, 1.200, 1.510, 1.000),
    atom_line(2, "ALA", "O", 9, 1.250, 1.610, 1.000),
    atom_line(2, "ALA", "CB", 10, 1.330, 1.350, 1.090),
    atom_line(5, "LEU", "N", 11, 2.500, 2.500, 2.500),
    atom_line(5, "LEU", "C", 12, 2.600, 2.550, 2.500),
    atom_line(5, "LEU", "O", 13, 2.640, 2.580, 2.470),
    atom_line(5, "LEU", "CB", 14, 2.450, 2.400, 2.560),
    atom_line(5, "LEU", "CG", 15, 2.360, 2.330, 2.520),
    atom_line(5, "LEU", "CD1", 16, 2.420, 2.250, 2.470),
    atom_line(5, "LEU", "CD2", 17, 2.270, 2.380, 2.560),
    atom_line(1, "VAL", "H1", 18, 0.900, 1.000, 1.000),
    atom_line(1, "VAL", "HA", 19, 1.149, 0.900, 1.000),
    atom_line(2, "ALA", "H", 20, 1.150, 1.160, 1.000),
    atom_line(2, "ALA", "HA", 21, 1.240, 1.380, 1.100),
    atom_line(2, "ALA", "HXT", 22, 1.200, 1.401, 1.000),
    atom_line(5, "LEU", "H1", 23, 2.400, 2.500, 2.500),
    atom_line(5, "LEU", "1HXT", 24, 2.550, 2.450, 2.530),
    atom_line(5, "LEU", "2HXT", 25, 2.650, 2.500, 2.570),
    atom_line(5, "LEU", "HB1", 26, 2.470, 2.340, 2.640),
    atom_line(5, "LEU", "HB2", 27, 2.530, 2.350, 2.510),
    atom_line(5, "LEU", "HB3", 28, 2.420, 2.460, 2.480),
    atom_line(6, "ALA", "N", 29, 3.200, 3.100, 3.000),
    atom_line(6, "ALA", "CA", 30, 3.300, 3.180, 3.040),
    atom_line(6, "ALA", "C", 31, 3.420, 3.150, 3.000),
    atom_line(6, "ALA", "O", 32, 3.470, 3.240, 2.980),
    atom_line(6, "ALA", "CB", 33, 3.330, 3.150, 3.190),
    atom_line(6, "ALA", "H1", 34, 3.150, 3.030, 2.950),
    atom_line(6, "ALA", "H2", 35, 3.230, 3.180, 2.930),
    atom_line(6, "ALA", "HA", 36, 3.320, 3.270, 3.090),
]

with tempfile.TemporaryDirectory() as directory:
    gro = Path(directory) / "debris.gro"
    write_gro(gro, ATOMS)
    _title, all_atoms, box = engine.parse_gro(gro)
    raw_blocks = engine.contiguous_residue_blocks(all_atoms, engine.AA_NAMES)

    # Default mode drops only the degenerate residue key.
    residues, caps, corrections, dropped = engine.separate_noncontiguous_hydrogen_caps(raw_blocks, box)
    assert len(dropped) == 1, dropped
    assert dropped[0].residue_key == (5, "LEU")
    assert dropped[0].carbonyl_hydrogen_names == ["1HXT", "2HXT"]
    assert dropped[0].atom_count == 13
    assert dropped[0].reason
    kept_keys = [(r.source_resid, r.source_resname) for r in residues]
    assert (5, "LEU") not in kept_keys
    assert kept_keys == [(1, "VAL"), (2, "ALA"), (6, "ALA")]
    # The clean carbonyl-H cap on residue 2 is still detected.
    assert [cap.kind for cap in caps[(2, "ALA")]] == ["C_HXT"]
    # The neutral NH2 start of the next segment is still detected.
    assert [cap.kind for cap in caps[(6, "ALA")]] == ["N_H"]
    assert corrections == []

    # Strict mode still refuses the debris residue.
    try:
        engine.separate_noncontiguous_hydrogen_caps(
            raw_blocks, box, degenerate_residue_mode="error"
        )
    except ValueError as exc:
        assert "Multiple hydrogens bonded to backbone carbonyl C" in str(exc)
        assert "(5, 'LEU')" in str(exc)
    else:
        raise AssertionError("error mode did not raise")

    # A residue with two carbonyl hydrogens but an intact CA is never dropped,
    # even in the default drop mode: that atom set is chemically impossible.
    gro2 = Path(directory) / "intact_ca.gro"
    write_gro(gro2, ATOMS + [atom_line(5, "LEU", "CA", 40, 2.600, 2.450, 2.600)])
    _title, all_atoms2, box2 = engine.parse_gro(gro2)
    raw_blocks2 = engine.contiguous_residue_blocks(all_atoms2, engine.AA_NAMES)
    try:
        engine.separate_noncontiguous_hydrogen_caps(raw_blocks2, box2)
    except ValueError as exc:
        assert "Multiple hydrogens bonded to backbone carbonyl C" in str(exc)
        assert "(5, 'LEU')" in str(exc)
    else:
        raise AssertionError("debris with an intact CA was silently accepted")

print("test_degenerate_residue_drop: OK")
