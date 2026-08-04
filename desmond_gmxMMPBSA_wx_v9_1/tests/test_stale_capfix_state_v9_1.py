#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "app" / "desmond_gmxmmpbsa.py"
spec = importlib.util.spec_from_file_location("workflow_engine_v9_1", ENGINE)
assert spec and spec.loader
engine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = engine
spec.loader.exec_module(engine)

with tempfile.TemporaryDirectory() as directory:
    workdir = Path(directory)
    base_pdb = workdir / "complex_amber_order.pdb"
    base_xtc = workdir / "complex_stride10.xtc"
    backup_pdb = workdir / "complex_amber_order.pdb.before_capfix"
    backup_xtc = workdir / "complex_stride10.xtc.before_capfix"
    fixed_pdb = workdir / "complex_amber_order_capfix.pdb"
    fixed_xtc = workdir / "complex_stride10_capfix.xtc"
    base_pdb.write_bytes(b"BASE-PDB")
    base_xtc.write_bytes(b"BASE-XTC")
    backup_pdb.write_bytes(b"BASE-PDB")
    backup_xtc.write_bytes(b"BASE-XTC")
    fixed_pdb.write_bytes(b"FIXED-PDB")
    fixed_xtc.write_bytes(b"FIXED-XTC")
    (workdir / "preparation_report.json").write_text(
        json.dumps(
            {
                "miniapp_version": "9.0.0",
                "files": {
                    "complex_pdb": base_pdb.name,
                    "complex_xtc": base_xtc.name,
                },
            }
        )
    )
    (workdir / "cap_clash_repair.json").write_text(
        json.dumps(
            {
                "validation_passed": True,
                "caps": [{"residue": "ACE1"}],
                "input_pdb": str(base_pdb),
                "input_xtc": str(base_xtc),
                "output_pdb": str(fixed_pdb),
                "output_xtc": str(fixed_xtc),
            }
        )
    )
    result = engine.apply_synthetic_cap_clash_repair(workdir)
    assert result["validation_passed"]
    assert result["recovered_from_overwritten_active_files"]
    report = json.loads((workdir / "preparation_report.json").read_text())
    assert report["files"]["complex_pdb"] == fixed_pdb.name
    assert report["files"]["complex_xtc"] == fixed_xtc.name
    assert report["miniapp_version"] == "9.1.0"
    # Base files remain immutable; only the manifest selects repaired data.
    assert base_pdb.read_bytes() == b"BASE-PDB"
    assert base_xtc.read_bytes() == b"BASE-XTC"

with tempfile.TemporaryDirectory() as directory:
    workdir = Path(directory)
    for name in (
        "cap_clash_repair.json",
        "complex_stride10_capfix.xtc",
        "FINAL_RESULTS_MMPBSA.dat",
    ):
        (workdir / name).write_text(name)
    archive = engine.archive_stale_preparation_state(workdir)
    assert archive is not None and archive.is_dir()
    assert (archive / "cap_clash_repair.json").is_file()
    assert (archive / "complex_stride10_capfix.xtc").is_file()
    assert (archive / "FINAL_RESULTS_MMPBSA.dat").is_file()
    assert not (workdir / "cap_clash_repair.json").exists()

preflight_source = (ROOT / "app" / "steric_preflight.py").read_text()
synthetic_check = 'and bool(result["unsafe_synthetic_cap_pair_present"])'
source_stop = 'if bool(result["unsafe_source_pair_present"]):'
assert synthetic_check in preflight_source
assert preflight_source.index(synthetic_check) < preflight_source.rindex(source_stop)

print("Version-9.1 stale cap-repair state regression tests passed")
