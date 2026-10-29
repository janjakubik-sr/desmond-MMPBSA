#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import tempfile
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "app" / "desmond_gmxmmpbsa.py"
spec = importlib.util.spec_from_file_location("workflow_engine", ENGINE)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

CASES = {
    "unicode": """
Calculations performed using 101 complex frames
POISSON BOLTZMANN:
Delta (Complex - Receptor - Ligand):
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
ΔTOTAL                  -42.125        1.250      3.500        0.125      0.350
""",
    "ascii": """
Calculations performed using 51 complex frames
GENERALIZED BORN:
DELTA TOTAL             -17.25         2.00       2.50         0.20       0.25
""",
    "bordered": """
PB CALCULATION
| Δ TOTAL               −31.50         1.10       2.20         0.11       0.22 |
""",
}

with tempfile.TemporaryDirectory() as tmp:
    tmp_path = Path(tmp)
    summaries = {}
    for name, text in CASES.items():
        result = tmp_path / f"{name}.dat"
        result.write_text(text, encoding="utf-8")
        summaries[name] = module.parse_binding_energy(result)

assert summaries["unicode"]["delta_total_kcal_per_mol"]["average"] == -42.125
assert summaries["unicode"]["matched_label"] == "ΔTOTAL"
assert summaries["unicode"]["solvent_model"] == "PB"
assert summaries["unicode"]["frames"] == 101
assert summaries["ascii"]["matched_label"] == "DELTA TOTAL"
assert summaries["ascii"]["solvent_model"] == "GB"
assert summaries["bordered"]["delta_total_kcal_per_mol"]["average"] == -31.5
print(json.dumps(summaries, indent=2))
print("Parser tests passed.")
