# Validation: GUI 1.0.3 / engine 7.1.0

## Static validation

- `desmond_gmxmmpbsa_wx.py` passes `python3 -m py_compile`.
- Previous wxPython 4.2.3 font fixes are retained.
- Previous Results-panel parent/sizer fixes are retained.
- The scientific engine and ff19SB resources are byte-identical to GUI 1.0.2.

## Conda/Miniforge discovery test

A temporary test installation was constructed with:

```text
$HOME/bin/conda                 activation-only failing helper
$HOME/miniforge3/bin/conda      functional mock Conda CLI
CONDA_EXE=$HOME/bin/conda
PATH=$HOME/bin:...
```

The version-1.0.3 resolver rejected both the selected helper and the bad
`CONDA_EXE`, then returned the functional Miniforge path. It also confirmed
that the candidate supported both:

```text
conda --version
conda run --help
```

## Installer hotfix test

`apply_miniforge_conda_hotfix.sh` was run on an unmodified GUI 1.0.2 source.
The patched file:

- compiled successfully;
- was byte-identical to the clean GUI 1.0.3 source;
- contained no remaining unvalidated `self.conda_file.get() or
  discover_conda()` invocation.

Interactive rendering was not exercised in this headless build environment.
