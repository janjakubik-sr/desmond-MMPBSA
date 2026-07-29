# Desmond GRO/XTC → gmx_MMPBSA wxPython miniapp

**GUI version:** 1.0.3  
**Scientific workflow engine:** 7.1.0, based on the validated version 7 mapping/cap-handling code

This package adds a wxPython desktop interface to the open-source workflow that:

1. accepts a full-system GRO/XTC trajectory converted from Desmond and a ligand MOL2;
2. extracts and maps the protein–ligand complex into Amber ff19SB/GAFF2 atom order;
3. handles the receptor-specific terminal/cap conventions supported by workflow version 7;
4. selects every *N*th frame, removes problematic periodic imaging, and fits the complex;
5. builds and validates Amber/GROMACS topologies;
6. runs single-trajectory gmx_MMPBSA with PB or GB solvent;
7. displays and saves the final ΔTOTAL statistics.

## wxPython 4.2.3 font-construction hotfix

wxPython's `Window.SetFont()` expects a `wx.Font` object.  GUI 1.0.0
passed `wx.FontInfo` directly, which raises a `TypeError` on Debian's
wxPython 4.2.3 package.  GUI 1.0.1 now constructs `wx.Font` explicitly
for the title, log label, and monospaced log control.

## wxWidgets parent/sizer hotfix

GUI 1.0.1 created the two Results-tab action buttons with the outer Results
panel as their parent, but managed them with a sizer attached to a nested
`result_actions` panel. wxWidgets correctly rejected that parent/sizer
mismatch during startup. GUI 1.0.2 creates both buttons as children of
`result_actions`, matching the sizer's containing window.


## Miniforge/Conda executable detection

GUI 1.0.3 fixes environment checks on systems where Miniforge is installed
under `~/miniforge3` but `PATH` contains a different helper named `conda`.
An activation-only helper cannot execute `conda run`; it produces errors such
as:

```text
ArgumentError: activate does not accept more than one argument
```

The GUI now validates both `--version` and `run --help`, rejects such helpers,
and searches these locations before the user's `PATH`:

```text
$CONDA_EXE
~/miniforge3/bin/conda
~/miniforge3/condabin/conda
~/mambaforge/...
~/miniconda3/...
~/anaconda3/...
```

A stale invalid path stored in
`~/.config/desmond_gmxmmpbsa_wx/settings.json` is replaced automatically by a
working Miniforge executable. The expected setting for the present Debian
installation is:

```text
$HOME/miniforge3/bin/conda
```

## Result-parser correction

Current gmx_MMPBSA output commonly labels the binding-energy row as Unicode `ΔTOTAL`. The earlier engine searched only for ASCII `DELTA TOTAL`, so a successful calculation could be followed by the false exception:

```text
ValueError: No DELTA TOTAL row found in FINAL_RESULTS_MMPBSA.dat
```

Engine 7.1.0 accepts both `ΔTOTAL` and `DELTA TOTAL`. Summary parsing is also post-processing: if a future output format cannot be parsed, completed gmx_MMPBSA result files are retained and the calculation is not reported as failed.

## Recommended Linux installation

Keep wxPython in the operating-system Python and the scientific tools in the existing Conda environment.

On Debian/Ubuntu:

```bash
sudo apt install python3-wxgtk4.0
```

The scientific Conda environment should already contain:

```text
numpy, mdtraj, networkx, parmed
AmberTools: antechamber, parmchk2, tleap, sander, cpptraj
GROMACS
gmx_MMPBSA
mpi4py/OpenMPI for parallel runs
```

## Run the packaged application

```bash
tar -xzf desmond_gmxMMPBSA_wx_v7_1_3.tar.gz
cd desmond_gmxMMPBSA_wx_v7_1_3
./run_desmond_gmxmmpbsa_wx.sh
```

The launcher searches for an interpreter with `wx` in this order:

1. `$WX_PYTHON`, when set;
2. `/usr/bin/python3`;
3. `python3` from `PATH`.

To specify one explicitly:

```bash
WX_PYTHON=/usr/bin/python3 ./run_desmond_gmxmmpbsa_wx.sh
```

The GUI defaults to executing the molecular-modelling engine through:

```bash
conda run --no-capture-output -n desmond-mmpbsa python -u ...
```

This lets the GUI use system wxPython without installing wxPython into the Conda environment.

## Patch the existing installation

For an existing installation at `$HOME/desmond_gmxMMPBSA`:

```bash
./install_wx_v7_1_3.sh $HOME/desmond_gmxMMPBSA
```

The installer creates timestamped backups, installs the patched 7.1.0 engine, and adds:

```text
$HOME/desmond_gmxMMPBSA/desmond_gmxmmpbsa_wx.py
$HOME/desmond_gmxMMPBSA/run_desmond_gmxmmpbsa_wx.sh
```

Launch it with:

```bash
$HOME/desmond_gmxMMPBSA/run_desmond_gmxmmpbsa_wx.sh
```

An optional desktop-menu entry can be installed from the final application directory:

```bash
./install_desktop_launcher.sh
```

## Calculation tab

The application supports four operations:

- **Complete workflow: prepare + run**
- **Prepare inputs only**
- **Run an existing prepared directory**
- **Summarize an existing result**

For a complete calculation select or enter:

- work/output directory;
- full-system GRO;
- matching full-system XTC;
- ligand MOL2;
- ligand residue name in the GRO;
- integer ligand charge, or leave blank to infer it;
- frame stride;
- PB or GB solvent model;
- number of MPI processes.

The **Detect from GRO/MOL2** button compares the MOL2 molecule name and atom count against residue names in the GRO. Detection remains advisory and the value is editable.

## Advanced tab

Available parameters include:

- MDTraj chunk size;
- preservation versus normalization of supplied ligand charges;
- ACE/NME handling of Desmond hydrogen-capped termini;
- handling of excess internal backbone-N hydrogen sites;
- ff19SB XML path;
- Conda executable and environment name;
- current-Python versus `conda run` execution;
- workflow-engine path;
- additional expert command-line arguments.

## Execution and monitoring

The GUI:

- validates inputs before launch;
- shows the exact command;
- streams stdout/stderr into a live monospaced log;
- tracks gmx_MMPBSA percentage output and calculation stage;
- writes `wx_gui_run.log` into the work directory;
- can terminate the complete subprocess group, including MPI workers;
- stores the last settings under `~/.config/desmond_gmxmmpbsa_wx/settings.json`;
- supports reusable JSON profiles.

## Results tab

After a successful run, the GUI loads `binding_energy_summary.json` or parses `FINAL_RESULTS_MMPBSA.dat` directly. It reports:

```text
solvent model
number of frames
ΔTOTAL average
SD(Prop.)
SD
SEM(Prop.)
SEM
matched result label
```

A result from an older workflow can be repaired without recalculation by choosing **Summarize an existing result** and selecting `FINAL_RESULTS_MMPBSA.dat`.

The command-line equivalent is:

```bash
conda activate desmond-mmpbsa
python app/desmond_gmxmmpbsa.py summarize \
    --result /path/to/FINAL_RESULTS_MMPBSA.dat
```

## Scientific scope

The calculation is an endpoint single-trajectory MM/PBSA or MM/GBSA estimate using conformations sampled in Desmond and rescored with Amber ff19SB plus GAFF2. It is not an OPLS4/Prime MM-GBSA reproduction and does not include an entropy term unless the input workflow is deliberately extended to calculate one.

## wxPython 4.2.3 compatibility hotfix

GUI 1.0.2 retains the GUI 1.0.1 correction and wraps every `wx.FontInfo` specification in a `wx.Font` object before passing it to `Window.SetFont`. This is required by wxPython Phoenix 4.2.3 and fixes the startup exception `Window.SetFont(): argument 1 has unexpected type 'FontInfo'`.
