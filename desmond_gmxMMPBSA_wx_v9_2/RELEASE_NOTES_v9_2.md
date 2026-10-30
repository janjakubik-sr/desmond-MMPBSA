Release notes — v9.2 (engine 9.1.4 / GUI 1.3.0), 2026-10-29–30
============================================================

Focus: remote execution on servers without the shared Data mount
-----------------------------------------------------------------

The remote Slurm submission module `remote_alkyone.py` was renamed to the
neutral `remote_run.py`. It now supports two staging modes:

1. **mount** (unchanged): the server reads the inputs from the same
   sshfs-mounted data directory as the local machine and copies the
   results back onto that mount itself.
2. **rsync** (new): for servers that cannot mount the data directory
   (e.g. `fgusrvvs03`). The local driver pushes the input files — or the
   whole prepared work directory for *Run an existing prepared directory* —
   to remote scratch with rsync before the job, and pulls the results and
   Slurm logs back afterwards. Scratch is removed on success and kept for
   inspection on failure.

New GUI fields (Advanced tab, persisted in settings/profiles):
**Remote SSH user** (optional `user@host` prefix for all ssh/rsync calls)
and **Remote staging** (rsync-mode checkbox). The submission summary and
the remote *Check environment* honour the SSH user.

An example profile for the `fgusrvvs03` server (user `jakubik`, partition
`CPU_job_Q`, scratch `/scratch/jakubik`, conda `~/miniforge3/bin/conda`,
rsync staging) is provided as `fgusrvvs03_profile.json`.

Engine 9.1.1 fixes Antechamber aborting on normalized MOL2 files with an
empty molecule-name line (zero bonds / open valences); the engine now
substitutes the destination file stem as the molecule name.

Verification: engine/GUI unit tests re-run for the new version strings
(`tests/test_v9_1_cap_state.py` failures are pre-existing from v9_1 and
unchanged); remote rsync staging was exercised end-to-end on `fgusrvvs03`
through its Slurm queue (`CPU_job_Q`), including input push, live log
streaming, result pull-back and scratch cleanup.

Post-release engine updates within v9.2 (2026-10-29–30)
--------------------------------------------------------

**Engine 9.1.2** adds degenerate chain-break residue handling: a residue
key covering disconnected fragments (carbonyl C with two capping
hydrogens, CA deleted) is dropped from the ff19SB receptor model in the
default `drop` mode or aborts in `error` mode, driven by the new
*Degenerate chain-break residue* GUI choice and
`--degenerate-residue-mode` flag.

**Engine 9.1.3/9.1.4** fix duplicate-H caps on internal peptide residues
(M4 5DSG MET345: an sp3-like backbone N bonded to the preceding carbonyl
yet carrying two N-H hydrogens). Bondedness now uses the segment
detector's C-N criterion (0.20 nm), amide geometry is advisory and
recorded, and the retained hydrogen is chosen by O-C-N-H planarity or
conventional naming. Regression tests use the real frames of the three
failing M4 jobs.

**remote_run.py** derives the driver/job script names from the SSH
target host (`fgusrvvs03_remote.sh` / `fgusrvvs03_job.sh`) instead of the
fixed `alkyone_*` names.

Verification: all engine/GUI unit tests pass for the new version strings;
`prepare` was exercised end-to-end on all three M4 replicates
(`desmond_md_job_1/2/3`), confirming the MET345 internal correction
(removed/retained hydrogens recorded in
`internal_backbone_hydrogen_corrections.tsv`) and the HIE24→ACE /
ARG427→NME terminal caps.
