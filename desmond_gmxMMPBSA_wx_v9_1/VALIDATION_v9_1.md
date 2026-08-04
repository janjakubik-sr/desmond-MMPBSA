# Validation — engine 9.1.0 / GUI 1.0.6

## Audited failure reproduced from the complete uploaded work directory

The uploaded M3/4DAJ–KH-5 work directory contained three successive runs. The first run completed; later runs reused the same directory and overwrote the active repaired coordinates.

Observed hashes:

```text
complex_stride10.xtc
28ecbfdb218447afc84d4da44694956af44156e06791b9299f0f77e3e4602120

complex_stride10.xtc.before_capfix
28ecbfdb218447afc84d4da44694956af44156e06791b9299f0f77e3e4602120

complex_stride10_capfix.xtc
a6dd6cda5980815e5c28cc06e0defd689c533026fbd4c57309ae42f09c292d19
```

Thus, the active XTC was byte-identical to its unrepaired backup and not to the validated capfixed XTC. The structure-route copy had the same unrepaired hash.

## Sander overflow localization

The retained receptor and complex mdout files identified the same synthetic-cap clashes:

| Prepared frame | Source frame | Pair | Before | After cap repair |
|---:|---:|---|---:|---:|
| 17 | 170 | TYR192 HE2 — ACE199 H3 | 0.027749 nm | 0.278614 nm |
| 28 | 280 | ARG103 HH21 — ACE199 H3 | 0.020273 nm | 0.236047 nm |
| 76 | 760 | ARG103 HG2 — ACE199 H2 | 0.038223 nm | 0.411147 nm |

The overflow was present in receptor output as well as complex output, excluding ligand cross-interactions as its source.

## Recovery test

On a copy of the uploaded directory, `apply_synthetic_cap_clash_repair()` detected that:

- active PDB/XTC matched the `.before_capfix` files;
- validated capfixed outputs existed;
- the old report could therefore be safely associated with the current preparation.

The manifest was updated to:

```text
complex_pdb = complex_amber_order_capfix.pdb
complex_xtc = complex_stride10_capfix.xtc
```

The base files remained unchanged.

## Automated tests

All packaged regression scripts completed successfully:

- atom mapping;
- Miniforge/Conda discovery;
- GUI version and SciPy dependency checks;
- Unicode/ASCII `ΔTOTAL` parsing;
- stale capfix state recovery and prior-run archiving;
- topology-aware preflight and undefined-energy diagnostics;
- wx parent/sizer checks.

The final AmberTools/gmx_MMPBSA numerical run cannot be executed in the artifact build container because it lacks AmberTools and GROMACS. The uploaded persistent GUI log records that the same capfixed trajectory completed a 101-frame PB run before the later overwrite.
