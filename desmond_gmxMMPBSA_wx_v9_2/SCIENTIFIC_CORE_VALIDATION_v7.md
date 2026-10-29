# Validation report: desmond_gmxMMPBSA v7.0.0

## Target M4/5DSG GRO

The uploaded M4 GRO was parsed directly.

| Check | Result |
|---|---:|
| Source protein residues | 290 |
| Source protein atoms | 4,681 |
| Internal excess N-H atoms omitted | 1 |
| Terminal carbonyl-H atoms replaced | 1 |
| Synthetic NME atoms added | 6 |
| Expected prepared protein atoms | 4,685 |
| Protein-template mapping failures | 0 |

### MET345 correction

| Quantity | Value |
|---|---:|
| GLN344 C--MET345 N | 0.137996 nm |
| GLN344 CA--C--N angle | 115.396 degrees |
| C--N--MET345 CA angle | 125.003 degrees |
| CA--C--N--CA omega | 167.290 degrees |
| Retained H | `2H` (source atom 3312) |
| Retained O--C--N--H | -179.406 degrees |
| Omitted H | `1H` (source atom 3311) |
| Omitted O--C--N--H | -70.073 degrees |

### ARG427 terminal cap

| Quantity | Value |
|---|---:|
| Source atom | `HC` (source atom 4659) |
| C--HC distance | 0.100215 nm |
| Replacement | synthetic ff19SB NME |

## Four-receptor GRO regression

| Input | Residues | Segments | Terminal cap replacements | Internal H corrections | Mapping failures |
|---|---:|---:|---:|---:|---:|
| M1 input | 284 | 2 | 4 | 0 | 0 |
| M2 input | 277 | 1 | 0 | 0 | 0 |
| M3 input | 264 | 2 | 2 | 0 | 0 |
| M4 input | 290 | 1 | 1 | 1 | 0 |

## Full preparation regressions

Two complete GRO/XTC/MOL2 preparation runs were executed with stride 100.

| System | Protein | Ligand | Complex | Frames | Maximum covalent bond |
|---|---:|---:|---:|---:|---:|
| M2/QNB | 4,513 | 49 | 4,562 | 11 | 0.21114 nm |
| M1/0HK | 4,650 | 48 | 4,698 | 11 | 0.21257 nm |

Both produced PDB, GRO, XTC, index, mapping, trajectory-validation, and
preparation-report files without regression errors.

## Limits of validation

The M4 KH-5 XTC and UNK MOL2 were not uploaded with the failing GRO, so the
complete M4 Antechamber/LEaP/gmx_MMPBSA energy stage could not be executed in
this runtime. The exact failed protein-mapping stage was reproduced and now
passes with zero unresolved or unused protein atoms.
