# Release notes — engine 9.1.0 / wxPython GUI 1.0.6

## Corrected stale cap-repair state

A fresh `all` run could regenerate `complex_amber_order.pdb` and the stride-selected XTC in an existing work directory while leaving an older successful `cap_clash_repair.json`. The repair wrapper returned immediately on `validation_passed=true`, even though the active files no longer contained the repaired coordinates. Later Sander calculations could then overflow `VDWAALS` at synthetic ACE/NME contacts.

Engine 9.1.0:

- archives previous cap-repair state and final result files before a fresh preparation;
- keeps repaired PDB/XTC files explicit and immutable rather than replacing the base prepared files;
- updates `preparation_report.json` to select the exact capfixed inputs;
- validates report/file provenance before reusing an existing repair;
- recovers compatible older folders by comparing active files with `.before_capfix` backups;
- attempts synthetic-cap repair before treating a frame-wide source-pair heuristic as fatal.

## Recovery utility

`run_capfixed_recovery.sh` validates an existing `cap_clash_repair.json`, confirms that the current active files were overwritten by the unrepaired preparation, creates a clean rerun directory, and executes single-trajectory gmx_MMPBSA using the validated `*_capfix` PDB/XTC.

## Preserved behavior

No atom-mapping, force-field, ligand-parameterization, PB/GB, or result-parser algorithms were changed. GUI 1.0.6 changes only the displayed engine/GUI version labels.
