#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_gmx_structure_route.sh WORKDIR [MPI_PROCESSES]

Runs gmx_MMPBSA by building Amber topologies directly from the prepared
PDB plus Antechamber MOL2. It deliberately omits -cp COM.top and bypasses
the ParmEd Amber -> GROMACS -> Amber topology round-trip.

WORKDIR must contain preparation_report.json and the prepared files created
by desmond_gmxmmpbsa.py. Results are written below:
  WORKDIR/gmx_structure_route/
USAGE
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage >&2
  exit 2
fi

WORKDIR=$(readlink -f "$1")
NP=${2:-1}

if [[ ! -d "$WORKDIR" ]]; then
  echo "ERROR: not a directory: $WORKDIR" >&2
  exit 2
fi
if [[ ! "$NP" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: MPI_PROCESSES must be a positive integer" >&2
  exit 2
fi

for exe in python gmx_MMPBSA parmchk2 tleap sander cpptraj; do
  if ! command -v "$exe" >/dev/null 2>&1; then
    echo "ERROR: required executable not found in PATH: $exe" >&2
    echo "Activate the desmond-mmpbsa Miniforge environment first." >&2
    exit 127
  fi
done
if (( NP > 1 )) && ! command -v mpirun >/dev/null 2>&1; then
  echo "ERROR: mpirun is required for NP=$NP" >&2
  exit 127
fi

REPORT="$WORKDIR/preparation_report.json"
if [[ ! -s "$REPORT" ]]; then
  echo "ERROR: missing $REPORT" >&2
  exit 2
fi

mapfile -t META < <(python - "$REPORT" <<'PY'
import json
import sys
from pathlib import Path

report_path = Path(sys.argv[1])
report = json.loads(report_path.read_text())
files = report.get("files", {})
required = {
    "complex_pdb": files.get("complex_pdb", "complex_amber_order.pdb"),
    "reference_pdb": files.get(
        "reference_pdb", files.get("complex_pdb", "complex_amber_order.pdb")
    ),
    "complex_xtc": files.get("complex_xtc", "complex_stride10.xtc"),
    "index": files.get("index", "index.ndx"),
    "mol2": files.get("parameterized_mol2", files.get("normalized_mol2")),
    "ligand": report.get("ligand_resname"),
    "frames": report.get("trajectory_frames_written"),
}
for key in ("mol2", "ligand", "frames"):
    if required[key] in (None, ""):
        raise SystemExit(f"preparation_report.json does not define {key}")
for key in (
    "complex_pdb", "reference_pdb", "complex_xtc", "index",
    "mol2", "ligand", "frames"
):
    print(required[key])
PY
)

COMPLEX_PDB=${META[0]}
REFERENCE_PDB=${META[1]}
COMPLEX_XTC=${META[2]}
INDEX_FILE=${META[3]}
LIGAND_MOL2=${META[4]}
LIGAND_NAME=${META[5]}
NFRAMES=${META[6]}

for rel in "$COMPLEX_PDB" "$REFERENCE_PDB" "$COMPLEX_XTC" "$INDEX_FILE" "$LIGAND_MOL2"; do
  if [[ ! -s "$WORKDIR/$rel" ]]; then
    echo "ERROR: required prepared file is missing or empty: $WORKDIR/$rel" >&2
    exit 2
  fi
done

RUN_DIR="$WORKDIR/gmx_structure_route"
mkdir -p "$RUN_DIR"

cp -f "$WORKDIR/$COMPLEX_PDB" "$RUN_DIR/complex.pdb"
cp -f "$WORKDIR/$REFERENCE_PDB" "$RUN_DIR/reference.pdb"
cp -f "$WORKDIR/$COMPLEX_XTC" "$RUN_DIR/complex.xtc"
cp -f "$WORKDIR/$INDEX_FILE" "$RUN_DIR/index.ndx"
cp -f "$WORKDIR/$LIGAND_MOL2" "$RUN_DIR/ligand.mol2"

python - "$WORKDIR/mmpbsa.in" "$RUN_DIR/mmpbsa_full.in" "$RUN_DIR/mmpbsa_test.in" "$NFRAMES" <<'PY'
import re
import sys
from pathlib import Path

source = Path(sys.argv[1])
full_out = Path(sys.argv[2])
test_out = Path(sys.argv[3])
nframes = int(sys.argv[4])

if source.is_file():
    text = source.read_text(encoding="utf-8", errors="replace")
else:
    text = f'''Protein-ligand MM/PBSA
&general
  sys_name="Desmond_rescoring",
  startframe=1,
  endframe={nframes},
  interval=1,
  PBRadii=3,
  verbose=1,
/
&pb
  istrng=0.150,
  radiopt=0,
/
'''

match = re.search(r"(?ims)(&general\b)(.*?)(^\s*/\s*$)", text)
if not match:
    raise SystemExit("Could not find an &general namelist in mmpbsa.in")

body = match.group(2)
for key in ("forcefields", "startframe", "endframe", "interval", "keep_files"):
    body = re.sub(
        rf"(?im)^\s*{key}\s*=.*?(?:,\s*)?$",
        "",
        body,
    )

insert = (
    "\n  startframe=1,\n"
    f"  endframe={nframes},\n"
    "  interval=1,\n"
    "  keep_files=2,\n"
    "  forcefields=\"leaprc.protein.ff19SB,leaprc.gaff2\",\n"
)
full = (
    text[:match.start()]
    + match.group(1)
    + insert
    + body
    + match.group(3)
    + text[match.end():]
)
full_out.write_text(full, encoding="utf-8")
test = re.sub(
    r"(?im)^\s*endframe\s*=\s*\d+\s*,?",
    "  endframe=1,",
    full,
)
test_out.write_text(test, encoding="utf-8")
PY

cd "$RUN_DIR"

clean_scratch() {
  rm -rf -- \
    _GMXMMPBSA_* COM_traj_*.xtc REC_traj_*.xtc LIG_traj_*.xtc \
    COM.prmtop REC.prmtop LIG.prmtop COM.inpcrd REC.inpcrd LIG.inpcrd \
    FINAL_RESULTS_MMPBSA_TEST.dat FINAL_RESULTS_MMPBSA_TEST.csv \
    FINAL_RESULTS_MMPBSA.dat FINAL_RESULTS_MMPBSA.csv \
    gmx_MMPBSA.log _MMPBSA_* 2>/dev/null || true
}

run_one() {
  local input_file=$1
  local out_dat=$2
  local out_csv=$3
  local log_file=$4
  local processes=${5:-$NP}
  local -a cmd=(
    gmx_MMPBSA -O
    -i "$input_file"
    -cs complex.pdb
    -cr reference.pdb
    -ci index.ndx
    -cg 1 2
    -ct complex.xtc
    -lm ligand.mol2
    -o "$out_dat"
    -eo "$out_csv"
    -nogui
  )
  if (( processes > 1 )); then
    cmd=(mpirun -np "$processes" "${cmd[@]}")
  fi
  printf '[RUN]'
  printf ' %q' "${cmd[@]}"
  printf '\n'
  "${cmd[@]}" 2>&1 | tee "$log_file"
}

clean_scratch

echo "=== One-frame topology/energy test using structure-built topology ==="
run_one \
  mmpbsa_test.in \
  FINAL_RESULTS_MMPBSA_TEST.dat \
  FINAL_RESULTS_MMPBSA_TEST.csv \
  gmx_MMPBSA_structure_test.log \
  1

if [[ ! -s FINAL_RESULTS_MMPBSA_TEST.dat ]]; then
  echo "ERROR: one-frame structure-route test did not produce a result file" >&2
  exit 1
fi

mkdir -p topology_from_structure_test
for f in COM.prmtop REC.prmtop LIG.prmtop COM.inpcrd REC.inpcrd LIG.inpcrd; do
  [[ -s "$f" ]] && cp -f "$f" topology_from_structure_test/
done

python - "$WORKDIR" "$RUN_DIR" <<'PY'
import json
import sys
from pathlib import Path

import parmed as pmd

workdir = Path(sys.argv[1])
run_dir = Path(sys.argv[2])
source_dir = workdir / "amber_reference"
generated_dir = run_dir / "topology_from_structure_test"
prep_report = json.loads((workdir / "preparation_report.json").read_text())
expected = {
    "COM": int(prep_report["complex_atoms"]),
    "REC": int(prep_report["protein_atoms"]),
    "LIG": int(prep_report["ligand_atoms"]),
}
report = {
    "reference_available": source_dir.is_dir(),
    "systems": {},
    "validation_passed": True,
}
errors = []

for prefix in ("COM", "REC", "LIG"):
    generated_path = generated_dir / f"{prefix}.prmtop"
    reference_path = source_dir / f"{prefix}.prmtop"
    if not generated_path.is_file():
        errors.append(f"Missing topology generated by gmx_MMPBSA: {generated_path}")
        continue

    generated = pmd.load_file(str(generated_path))
    item = {
        "generated_atoms": len(generated.atoms),
        "expected_atoms": expected[prefix],
        "generated_charge": float(sum(atom.charge for atom in generated.atoms)),
    }
    if len(generated.atoms) != expected[prefix]:
        errors.append(
            f"{prefix}: generated {len(generated.atoms)} atoms; "
            f"expected {expected[prefix]}"
        )

    if reference_path.is_file():
        reference = pmd.load_file(str(reference_path))
        generated_signature = [
            (atom.name, atom.atomic_number) for atom in generated.atoms
        ]
        reference_signature = [
            (atom.name, atom.atomic_number) for atom in reference.atoms
        ]
        order_matches = generated_signature == reference_signature
        max_charge_difference = (
            max(
                abs(float(a.charge) - float(b.charge))
                for a, b in zip(generated.atoms, reference.atoms)
            )
            if len(generated.atoms) == len(reference.atoms)
            else None
        )
        item.update(
            {
                "reference_atoms": len(reference.atoms),
                "reference_charge": float(
                    sum(atom.charge for atom in reference.atoms)
                ),
                "atom_name_element_order_matches": order_matches,
                "maximum_charge_difference": max_charge_difference,
            }
        )
        if not order_matches:
            errors.append(
                f"{prefix}: atom name/element order differs from amber_reference"
            )
        if (
            max_charge_difference is not None
            and max_charge_difference > 1.0e-5
        ):
            errors.append(
                f"{prefix}: maximum per-atom charge difference "
                f"{max_charge_difference:.6g} exceeds 1e-5"
            )

    report["systems"][prefix] = item

report["errors"] = errors
report["validation_passed"] = not errors
(run_dir / "structure_route_topology_audit.json").write_text(
    json.dumps(report, indent=2) + "\n",
    encoding="utf-8",
)
print(json.dumps(report, indent=2))
if errors:
    raise SystemExit(
        "Structure-built topology audit failed: " + "; ".join(errors)
    )
PY

clean_scratch

echo "=== Full calculation using structure-built topology ==="
run_one \
  mmpbsa_full.in \
  FINAL_RESULTS_MMPBSA.dat \
  FINAL_RESULTS_MMPBSA.csv \
  gmx_MMPBSA_structure_full.log \
  "$NP"

echo
echo "Completed. Results:"
echo "  $RUN_DIR/FINAL_RESULTS_MMPBSA.dat"
echo "  $RUN_DIR/FINAL_RESULTS_MMPBSA.csv"
echo "  $RUN_DIR/structure_route_topology_audit.json"
