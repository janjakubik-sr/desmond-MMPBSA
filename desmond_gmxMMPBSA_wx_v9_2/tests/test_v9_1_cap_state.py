from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import types
from pathlib import Path

import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "app"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def write_tiny_fixture(workdir: Path) -> None:
    pdb = workdir / "complex_amber_order.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       1.900  -0.800   1.200  1.00  0.00           C\n"
        "TER\nEND\n",
        encoding="utf-8",
    )
    traj = md.load(str(pdb))
    traj.save_xtc(str(workdir / "complex_stride10.xtc"))
    (workdir / "covalent_bonds.tsv").write_text(
        "atom1_1based\tatom2_1based\n1\t2\n2\t3\n3\t4\n2\t5\n",
        encoding="utf-8",
    )
    (workdir / "preparation_report.json").write_text(
        json.dumps(
            {
                "complex_atoms": 5,
                "trajectory_frames_written": 1,
                "source_frame_indices": [0],
                "files": {
                    "complex_pdb": "complex_amber_order.pdb",
                    "complex_xtc": "complex_stride10.xtc",
                },
            }
        ),
        encoding="utf-8",
    )


def test_frame_repulsion_classification_prefers_dominant_synthetic_component(monkeypatch):
    monkeypatch.setitem(sys.modules, "parmed", types.ModuleType("parmed"))
    module = load_module("steric_preflight_test", APP / "steric_preflight.py")
    source, synthetic = module.classify_frame_repulsion(
        total_positive=31_466_300_908.0,
        source_positive=83_000.0,
        synthetic_positive=31_466_217_908.0,
        frame_limit=100_000.0,
    )
    assert source is False
    assert synthetic is True


def test_repair_report_records_active_hashes(tmp_path: Path):
    write_tiny_fixture(tmp_path)
    command = [
        sys.executable,
        str(APP / "repair_synthetic_caps.py"),
        "--workdir",
        str(tmp_path),
        "--allow-no-caps",
        "--backup-and-replace",
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    report = json.loads((tmp_path / "cap_clash_repair.json").read_text())
    assert report["validation_passed"] is True
    assert report["input_pdb_sha256_before"] == report["active_pdb_sha256_after"]
    assert report["input_xtc_sha256_before"] == report["active_xtc_sha256_after"]


def test_engine_does_not_trust_stale_success_report(tmp_path: Path):
    write_tiny_fixture(tmp_path)
    (tmp_path / "cap_clash_repair.json").write_text(
        json.dumps({"validation_passed": True}), encoding="utf-8"
    )
    engine = load_module("engine_v9_1_test", APP / "desmond_gmxmmpbsa.py")
    report = engine.apply_synthetic_cap_clash_repair(tmp_path)
    assert report["validation_passed"] is True
    assert "active_pdb_sha256_after" in report
    assert "active_xtc_sha256_after" in report
