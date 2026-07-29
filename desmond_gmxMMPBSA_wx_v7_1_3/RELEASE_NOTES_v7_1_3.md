# Release notes: GUI 1.0.3 / engine 7.1.0

## Fixed: Miniforge environment check selected an activation shim

GUI 1.0.2 accepted the first file named `conda` found through `PATH` without
checking whether it was the actual Conda command-line executable. On systems
where `~/bin/conda` is an activation helper, the GUI constructed:

```text
~/bin/conda run --no-capture-output ...
```

and the helper interpreted every argument as an argument to `conda activate`.
The resulting error was:

```text
ArgumentError: activate does not accept more than one argument
```

GUI 1.0.3:

- validates `conda --version`;
- validates `conda run --help`;
- prioritizes `~/miniforge3/bin/conda` and other standard Conda-family
  installation paths before `PATH`;
- rejects activation-only wrappers and aliases exposed as files;
- repairs a stale invalid Conda path loaded from the GUI settings file;
- updates the Advanced-tab label to `Conda/Miniforge executable`.

The scientific engine remains version 7.1.0 and is unchanged.
