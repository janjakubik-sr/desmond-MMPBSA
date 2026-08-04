#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_capfixed_recovery.sh WORKDIR [MPI_PROCESSES] [RERUN_DIR]

Example:
  ./run_capfixed_recovery.sh \
    /home/roshi/Data/.../LIG_mmpbsa \
    8

The script does not alter WORKDIR. It validates the previously generated
*_capfix.pdb and *_capfix.xtc files, copies the required inputs into a clean
subdirectory, and runs the single-trajectory gmx_MMPBSA calculation there.
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" || $# -lt 1 || $# -gt 3 ]]; then
  usage
  [[ $# -ge 1 ]] && exit 0 || exit 2
fi

WORKDIR=$(readlink -f "$1")
NP=${2:-8}
RERUN_DIR=${3:-"$WORKDIR/capfixed_recovery"}
RERUN_DIR=$(readlink -m "$RERUN_DIR")

if [[ ! -d "$WORKDIR" ]]; then
  echo "ERROR: work directory not found: $WORKDIR" >&2
  exit 2
fi
if ! [[ "$NP" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: MPI_PROCESSES must be a positive integer" >&2
  exit 2
fi

# Prefer the active environment. Fall back to the user's known Miniforge env.
find_program() {
  local name=$1
  if command -v "$name" >/dev/null 2>&1; then
    command -v "$name"
    return 0
  fi
  local candidate
  for candidate in \
    "$HOME/miniforge3/envs/desmond-mmpbsa/bin/$name" \
    "$HOME/miniconda3/envs/desmond-mmpbsa/bin/$name"; do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

PYTHON=$(find_program python) || {
  echo "ERROR: Python was not found in the active shell or desmond-mmpbsa environment" >&2
  exit 3
}
GMXMMPBSA=$(find_program gmx_MMPBSA) || {
  echo "ERROR: gmx_MMPBSA was not found" >&2
  exit 3
}
# Make an absolute Miniforge/Conda fallback behave like an activated environment.
ENV_BIN=$(dirname "$GMXMMPBSA")
ENV_ROOT=$(dirname "$ENV_BIN")
export AMBERHOME=${AMBERHOME:-$ENV_ROOT}
export PATH="$ENV_BIN:$ENV_ROOT/bin.AVX2_256:$PATH"
if (( NP > 1 )); then
  MPIRUN=$(find_program mpirun) || {
    echo "ERROR: mpirun was not found" >&2
    exit 3
  }
fi

mkdir -p "$RERUN_DIR"

# Validate the cap-repair provenance and create a clean, immutable rerun input set.
"$PYTHON" - "$WORKDIR" "$RERUN_DIR" <<'PY'
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

workdir = Path(sys.argv[1]).resolve()
outdir = Path(sys.argv[2]).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_required(source: Path, destination_name: str) -> Path:
    if not source.is_file():
        raise FileNotFoundError(source)
    destination = outdir / destination_name
    shutil.copy2(source, destination)
    return destination

manifest_path = workdir / "preparation_report.json"
repair_path = workdir / "cap_clash_repair.json"
if not manifest_path.is_file():
    raise FileNotFoundError(manifest_path)
if not repair_path.is_file():
    raise FileNotFoundError(repair_path)

manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
repair = json.loads(repair_path.read_text(encoding="utf-8"))
if not repair.get("validation_passed", False):
    raise RuntimeError("cap_clash_repair.json does not report validation_passed=true")

files = manifest.get("files", {})
if not isinstance(files, dict):
    raise RuntimeError("preparation_report.json has no valid files object")

active_pdb = workdir / Path(str(files.get("complex_pdb", "complex_amber_order.pdb"))).name
active_xtc = workdir / Path(str(files.get("complex_xtc", "complex_stride10.xtc"))).name
repaired_pdb = workdir / Path(str(repair["output_pdb"])).name
repaired_xtc = workdir / Path(str(repair["output_xtc"])).name
backup_pdb = active_pdb.with_suffix(active_pdb.suffix + ".before_capfix")
backup_xtc = active_xtc.with_suffix(active_xtc.suffix + ".before_capfix")

for path in (active_pdb, active_xtc, repaired_pdb, repaired_xtc):
    if not path.is_file():
        raise FileNotFoundError(path)

active_pdb_hash = sha256(active_pdb)
active_xtc_hash = sha256(active_xtc)
repaired_pdb_hash = sha256(repaired_pdb)
repaired_xtc_hash = sha256(repaired_xtc)

state = "unknown"
if active_pdb_hash == repaired_pdb_hash and active_xtc_hash == repaired_xtc_hash:
    state = "active files are already repaired"
elif backup_pdb.is_file() and backup_xtc.is_file():
    if active_pdb_hash == sha256(backup_pdb) and active_xtc_hash == sha256(backup_xtc):
        state = "active files were overwritten by the unrepaired preparation"

if state == "unknown":
    raise RuntimeError(
        "The active files match neither the validated capfixed files nor the recorded "
        "pre-capfix backups. Refusing to combine inputs from different preparations."
    )

# Basic PDB atom-count consistency.
def pdb_atoms(path: Path) -> int:
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        return sum(line.startswith(("ATOM  ", "HETATM")) for line in handle)

active_atoms = pdb_atoms(active_pdb)
repaired_atoms = pdb_atoms(repaired_pdb)
expected_atoms = int(manifest.get("complex_atoms", repaired_atoms))
if active_atoms != expected_atoms or repaired_atoms != expected_atoms:
    raise RuntimeError(
        f"PDB atom-count mismatch: active={active_atoms}, repaired={repaired_atoms}, "
        f"manifest={expected_atoms}"
    )

outdir.mkdir(parents=True, exist_ok=True)
# Remove only scratch/results from a prior recovery attempt.
for pattern in (
    "_GMXMMPBSA_*", "COM_traj_*.xtc", "REC_traj_*.xtc", "LIG_traj_*.xtc",
    "FINAL_RESULTS_MMPBSA.dat", "FINAL_RESULTS_MMPBSA.csv",
    "binding_energy_summary.json", "gmx_MMPBSA.log",
):
    for path in outdir.glob(pattern):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

copy_required(repaired_pdb, "complex_capfixed.pdb")
copy_required(repaired_xtc, "complex_capfixed.xtc")
copy_required(workdir / Path(str(files.get("reference_pdb", "complex_reference_noH.pdb"))).name,
              "reference.pdb")
copy_required(workdir / Path(str(files.get("index", "index.ndx"))).name, "index.ndx")
copy_required(workdir / "COM.top", "COM.top")
copy_required(workdir / Path(str(files.get("mmpbsa_input", "mmpbsa.in"))).name,
              "mmpbsa.in")

recovery = {
    "source_workdir": str(workdir),
    "recovery_workdir": str(outdir),
    "detected_state": state,
    "cap_repair_validation_passed": True,
    "complex_atoms": expected_atoms,
    "active_pdb_sha256": active_pdb_hash,
    "active_xtc_sha256": active_xtc_hash,
    "repaired_pdb_sha256": repaired_pdb_hash,
    "repaired_xtc_sha256": repaired_xtc_hash,
    "minimum_cap_contact_before_nm": repair.get("minimum_cap_nonbonded_contact_before", {}).get("distance_nm"),
    "minimum_cap_contact_after_nm": repair.get("minimum_cap_nonbonded_contact_after", {}).get("distance_nm"),
    "prepared_frames": manifest.get("trajectory_frames_written"),
    "stride": manifest.get("stride"),
}
(outdir / "recovery_manifest.json").write_text(
    json.dumps(recovery, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(recovery, indent=2))
PY

cd "$RERUN_DIR"

CMD=(
  "$GMXMMPBSA" -O
  -i mmpbsa.in
  -cs complex_capfixed.pdb
  -cr reference.pdb
  -ci index.ndx
  -cg 1 2
  -ct complex_capfixed.xtc
  -cp COM.top
  -o FINAL_RESULTS_MMPBSA.dat
  -eo FINAL_RESULTS_MMPBSA.csv
  -nogui
)
if (( NP > 1 )); then
  CMD=("$MPIRUN" -np "$NP" "${CMD[@]}")
fi

printf '\nRunning in %s\n' "$RERUN_DIR"
printf 'Command:'
printf ' %q' "${CMD[@]}"
printf '\n\n'

set +e
"${CMD[@]}" 2>&1 | tee gmx_MMPBSA_capfixed.log
status=${PIPESTATUS[0]}
set -e
if (( status != 0 )); then
  echo "ERROR: gmx_MMPBSA exited with status $status" >&2
  exit "$status"
fi

"$PYTHON" - FINAL_RESULTS_MMPBSA.dat binding_energy_summary.json <<'PY'
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

result = Path(sys.argv[1])
output = Path(sys.argv[2])
text = result.read_text(encoding="utf-8", errors="replace")
number = r"[-+−]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
matched = None
for line_number, line in enumerate(text.splitlines(), start=1):
    normalized = line.replace("−", "-").replace("Δ", "DELTA ")
    normalized = re.sub(r"DELTA\s*_?\s*TOTAL", "DELTA TOTAL", normalized, flags=re.I)
    if "DELTA TOTAL" not in normalized.upper():
        continue
    values = [float(token.replace("−", "-")) for token in re.findall(number, line)]
    if values:
        matched = (line_number, line.strip(), values)
if matched is None:
    raise SystemExit(f"No DELTA TOTAL/ΔTOTAL row found in {result}")
line_number, source_line, values = matched
keys = ["average", "standard_deviation_propagated", "standard_deviation",
        "standard_error_propagated", "standard_error"]
stats = {key: (values[i] if i < len(values) else None) for i, key in enumerate(keys)}
summary = {
    "result_file": str(result.resolve()),
    "matched_line_number": line_number,
    "matched_line": source_line,
    "delta_total_kcal_per_mol": stats,
}
output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print("\nBinding-energy summary")
print(json.dumps(summary, indent=2))
PY

echo ""
echo "Completed successfully."
echo "Results: $RERUN_DIR/FINAL_RESULTS_MMPBSA.dat"
echo "Summary: $RERUN_DIR/binding_energy_summary.json"
