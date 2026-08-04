#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "app" / "desmond_gmxmmpbsa.py"
PREFLIGHT = ROOT / "app" / "steric_preflight.py"

# Engine diagnostics do not need ParmEd at import time.
spec = importlib.util.spec_from_file_location("workflow_engine_v9", ENGINE)
assert spec and spec.loader
engine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = engine
spec.loader.exec_module(engine)
assert engine.__version__ == "9.1.0"

with tempfile.TemporaryDirectory() as directory:
    workdir = Path(directory)
    mdout = workdir / "_GMXMMPBSA_complex_pb.mdout.3"
    mdout.write_text(
        """ NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER\n"
        "     1       3.4933E+17     2.5182E+17     3.4933E+19     CB       5614\n"
        " BOND    =  1010877.6344  ANGLE   =     2623.4665  DIHED      =     9028.8882\n"
        " VDWAALS = *************  EEL     =   -17036.0135  EGB        =    -2042.6948\n"
        " 1-4 VDW =     5128.4701  1-4 EEL =    11580.9554\n"""
    )
    report = engine.diagnose_undefined_energy_outputs(workdir)
    assert report["undefined_energy_fields"]
    assert report["undefined_energy_fields"][0]["file"] == mdout.name
    assert "VDWAALS" in report["undefined_energy_fields"][0]["line"]
    assert (workdir / "undefined_energy_report.json").is_file()
    reloaded = json.loads((workdir / "undefined_energy_report.json").read_text())
    assert reloaded["undefined_energy_fields"]

# Import the topology-aware helper with a minimal ParmEd placeholder.  The
# numerical LJ function and topology graph construction are pure Python.
fake_parmed = types.ModuleType("parmed")
fake_parmed.load_file = lambda _path: None
sys.modules.setdefault("parmed", fake_parmed)
spec2 = importlib.util.spec_from_file_location("steric_preflight_v9", PREFLIGHT)
assert spec2 and spec2.loader
preflight = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = preflight
spec2.loader.exec_module(preflight)

class Atom:
    def __init__(self, idx: int, epsilon: float = 0.1, rmin: float = 1.5):
        self.idx = idx
        self.epsilon = epsilon
        self.rmin = rmin
        self.exclusion_partners = []

class Term:
    def __init__(self, *atoms):
        names = ("atom1", "atom2", "atom3", "atom4")
        for name, atom in zip(names, atoms):
            setattr(self, name, atom)

class DType:
    scnb = 2.0

atoms = [Atom(i) for i in range(4)]
bond = Term(atoms[0], atoms[1])
angle = Term(atoms[0], atoms[1], atoms[2])
dihedral = Term(atoms[0], atoms[1], atoms[2], atoms[3])
dihedral.improper = False
dihedral.type = DType()
atoms[0].exclusion_partners = [atoms[1], atoms[2], atoms[3]]
structure = types.SimpleNamespace(
    atoms=atoms,
    bonds=[bond],
    angles=[angle],
    dihedrals=[dihedral],
)
excluded, one_four = preflight.topology_sets(structure)
assert (0, 1) in excluded
assert (0, 2) in excluded
assert (0, 3) not in excluded
assert one_four[(0, 3)] == 2.0

# A 0.05-nm contact between ordinary Amber LJ atoms must be unmistakably unsafe.
energy = preflight.lj_energy_kcal(atoms[0], atoms[3], 0.05)
assert energy > 1.0e4

# Exercise the full frame scanner with a synthetic ACE/source-atom clash.
import mdtraj as md
import numpy as np

mdtop = md.Topology()
chain = mdtop.add_chain()
ace = mdtop.add_residue("ACE", chain, resSeq=1)
ala = mdtop.add_residue("ALA", chain, resSeq=2)
mdtop.add_atom("CH3", md.element.carbon, ace)
mdtop.add_atom("CA", md.element.carbon, ala)
trajectory = md.Trajectory(
    xyz=np.asarray([[[0.0, 0.0, 0.0], [0.02, 0.0, 0.0]]], dtype=np.float32),
    topology=mdtop,
)
scan_structure = types.SimpleNamespace(
    atoms=[Atom(0), Atom(1)], bonds=[], angles=[], dihedrals=[]
)
with tempfile.TemporaryDirectory() as directory:
    workdir = Path(directory)
    (workdir / "protein_templates.tsv").write_text(
        "output_residue\tsource_type\n1\tsynthetic_cap\n"
    )
    old_load_file = preflight.pmd.load_file
    old_load_xtc = preflight.md.load_xtc
    try:
        preflight.pmd.load_file = lambda _path: scan_structure
        preflight.md.load_xtc = lambda _path, top=None: trajectory
        result = preflight.scan(
            workdir=workdir,
            pdb=workdir / "mock.pdb",
            xtc=workdir / "mock.xtc",
            prmtop=workdir / "COM.prmtop",
            report={"source_frame_indices": [100], "protein_atoms": 2},
            cutoff_nm=0.25,
            pair_limit_kcal=1.0e4,
            frame_limit_kcal=1.0e5,
            hard_distance_nm=0.055,
            top_pairs_n=10,
        )
    finally:
        preflight.pmd.load_file = old_load_file
        preflight.md.load_xtc = old_load_xtc
    assert not result["validation_passed"]
    assert result["unsafe_synthetic_cap_pair_present"]
    assert not result["unsafe_source_pair_present"]
    assert result["unsafe_source_frames_0based"] == [100]
    assert result["top_pairs"][0]["synthetic_cap_involved"]

source = ENGINE.read_text(encoding="utf-8")
validation_pos = source.index("topology_validation = validate_and_convert_topologies")
preflight_pos = source.index("steric_validation = run_topology_aware_steric_preflight")
command_pos = source.index('"-cp", "COM.top"')
assert validation_pos < preflight_pos < command_pos
assert '"-rp", "REC.top"' not in source
assert '"-lp", "LIG.top"' not in source

print("Version-9.1 preflight and undefined-energy diagnostic tests passed")
