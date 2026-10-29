# Diagnosis of the uploaded `LIG_mmpbsa` work directory

## Executive finding

The failure is not caused by ligand mapping, the Amber-to-GROMACS topology conversion, or the structure-built topology route. The work directory contains files from three successive runs. The first run used repaired synthetic-cap coordinates and completed successfully. A later `all` run regenerated the base PDB/XTC in the same directory, overwriting the active repaired coordinates, while leaving an older `cap_clash_repair.json` with `validation_passed=true`. The engine trusted that stale report and skipped cap repair. Both the normal and structure-built routes then used the unrepaired trajectory.

## File-state proof

| File | SHA-256 | Interpretation |
|---|---|---|
| `complex_stride10.xtc` | `28ecbfdb218447afc84d4da44694956af44156e06791b9299f0f77e3e4602120` | active, unrepaired |
| `complex_stride10.xtc.before_capfix` | `28ecbfdb218447afc84d4da44694956af44156e06791b9299f0f77e3e4602120` | identical pre-repair backup |
| `complex_stride10_capfix.xtc` | `a6dd6cda5980815e5c28cc06e0defd689c533026fbd4c57309ae42f09c292d19` | validated repaired trajectory |
| `complex_amber_order.pdb` | `75643e87c1a5df7d9772cd1bbff88dea4690cc33682f6f3cc376bbe02ecb2326` | active, unrepaired |
| `complex_amber_order.pdb.before_capfix` | `75643e87c1a5df7d9772cd1bbff88dea4690cc33682f6f3cc376bbe02ecb2326` | identical pre-repair backup |
| `complex_amber_order_capfix.pdb` | `cd9cb65063711519a968bd0ec0591d807fe05630b02b88434491cf8b6cb38ecb` | validated repaired structure |

The structure-route copy `gmx_structure_route/complex.xtc` has the same SHA-256 as the unrepaired active XTC, so the route change could not solve the failure.

## Exact Sander overflows

The overflows occur in both complex and receptor outputs, excluding receptor-ligand cross interactions as the primary cause.

| Prepared frame | Original source frame | Pair | Unrepaired distance | Repaired distance |
|---:|---:|---|---:|---:|
| 17 | 170 | TYR192 HE2 (3141) — synthetic ACE199 H3 (3243) | 0.027749 nm | 0.278614 nm |
| 28 | 280 | ARG103 HH21 (1653) — synthetic ACE199 H3 (3243) | 0.020273 nm | 0.236047 nm |
| 76 | 760 | ARG103 HG2 (1641) — synthetic ACE199 H2 (3242) | 0.038223 nm | 0.411147 nm |

These pairs produce the `VDWAALS = *************` records in the retained Sander mdout files. The cap-repair report states a global minimum cap contact of 0.020273 nm before repair and 0.112883 nm after repair, with 150 cap/frame combinations modified and no non-cap coordinates changed.

## Run chronology

1. **2026-07-30 16:03:24:** preparation, cap repair, topology generation, and 101-frame PB calculation completed successfully. The persistent GUI log records `ΔTOTAL = -24.02 kcal/mol`, SD 4.41 kcal/mol, SEM 0.44 kcal/mol.
2. **2026-07-31 14:02:45:** another `all` command reused the same work directory. Preparation overwrote the active repaired PDB/XTC with new unrepaired base files. The stale successful repair report caused the repair wrapper to return without reapplying the repaired coordinates. Sander overflowed.
3. **2026-07-31 16:02:12:** the topology-aware preflight correctly detected severe synthetic-cap overlaps and stopped. Its JSON says `validation_passed=false`.
4. **Structure route:** the one-frame test passed because the selected first frame was finite. The 101-frame run copied and used the same unrepaired PDB/XTC and failed at the same later frames.

## Correct recovery

Run gmx_MMPBSA against `complex_amber_order_capfix.pdb` and `complex_stride10_capfix.xtc`, preferably in a clean subdirectory. Do not use the current active `complex_amber_order.pdb` and `complex_stride10.xtc` for this work directory.

The accompanying `run_capfixed_recovery.sh` script validates the provenance and checksums, creates a clean rerun directory, uses only the validated capfixed files, and leaves the original work directory unchanged.

## Permanent application correction

Engine 9.1.0 changes state handling as follows:

- a fresh `prepare` archives prior repair reports, capfixed coordinates, and final results instead of silently reusing or deleting them;
- repaired PDB/XTC files remain explicit immutable inputs and are selected through `preparation_report.json`;
- a successful repair report is reused only when its files can be tied to the current preparation by content checks;
- old folders where active files equal the pre-capfix backups are recovered by selecting the validated capfix outputs;
- the topology-aware preflight attempts synthetic-cap-only repair before treating a frame-wide source-pair heuristic as fatal.
