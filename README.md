# Desmond GRO/XTC → gmx_MMPBSA wxPython miniapp

**GUI version:** 1.0.6  
**Scientific workflow engine:** 9.1.0, retaining the validated version-9 topology-aware preflight and version-8 atom-mapping workflow

This package adds a wxPython desktop interface to the open-source workflow that:

1. accepts a full-system GRO/XTC trajectory converted from Desmond and a ligand MOL2;
2. extracts and maps the protein–ligand complex into Amber ff19SB/GAFF2 atom order;
3. handles Desmond/Schrodinger terminal-cap naming variants, hydrogen-capped chain breaks, and explicitly audited truncated side chains;
4. selects every *N*th frame, removes problematic periodic imaging, and fits the complex;
5. builds and validates Amber/GROMACS topologies;
6. runs single-trajectory gmx_MMPBSA with PB or GB solvent;
7. performs a topology-aware steric preflight with the actual Amber Lennard-Jones parameters and exclusions;
8. runs single-trajectory gmx_MMPBSA and displays/saves the final ΔTOTAL statistics.



## Engine 9.1 work-directory state correction

Engine 9.1 fixes a reproducible state-management defect found by auditing a
complete failed work directory. A successful cap-repair report could remain in
an output folder after a later `all` run regenerated the base PDB/XTC. The old
report still said `validation_passed=true`, so the repair wrapper returned
without reapplying or selecting the validated capfixed coordinates. Both the
standard and structure-built routes then used the newly regenerated,
unrepaired trajectory.

The corrected behavior is:

- repaired coordinates remain explicit `*_capfix.pdb` and `*_capfix.xtc`
  inputs; the base prepared coordinates are not silently replaced;
- `preparation_report.json` records the exact repaired files selected for the
  calculation;
- a fresh `prepare` into an existing directory archives prior repair state and
  final results under `previous_run_state/<timestamp>/`;
- a successful old report is reused only when its output files can be tied to
  the current preparation by content checks;
- compatible older folders are recovered when the active files are identical
  to the recorded `.before_capfix` backups;
- topology-aware preflight attempts synthetic-cap-only repair before a
  frame-wide source-pair heuristic is treated as fatal.

For an already affected directory, use:

```bash
./run_capfixed_recovery.sh /path/to/LIG_mmpbsa 8
```

This creates `/path/to/LIG_mmpbsa/capfixed_recovery`, uses the validated
capfixed PDB/XTC directly, and leaves the original directory unchanged.

## Engine 9 topology-aware steric preflight

Engine 9 addresses calculations that finish all Sander complex, receptor, and
ligand evaluations but fail while gmx_MMPBSA parses an undefined `VDWAALS`
field.  This normally means that a retained Sander mdout contains an overflow
field such as `VDWAALS = *************`.  A simple atom-count check or
covalent-bond check cannot detect every such case.

After LEaP creates `COM.prmtop`, but before the expensive endpoint calculation,
the engine now scans every prepared trajectory frame using:

- the actual Amber per-atom Lennard-Jones radii and well depths;
- the Amber 1-2/1-3 exclusion graph;
- scaled 1-4 interactions and their `SCNB` factors;
- the exact prepared atom order and source-frame mapping.

The preflight writes:

```text
steric_preflight_before.json
steric_clashes_before.tsv
steric_preflight.json
steric_clashes.tsv
steric_preflight.stdout.log
```

Only residues explicitly recorded as **synthetic ACE/NME caps** may be moved
automatically.  Source protein and ligand coordinates are never altered by the
preflight.  If every unsafe pair involves a synthetic cap, the cap is repaired
and the topology-aware scan is repeated.  If an unsafe pair contains a source
protein or ligand atom, the workflow stops before gmx_MMPBSA and reports the
exact prepared frame, original source frame, atom indices, atom labels,
distance, interaction class, and estimated pair Lennard-Jones energy.

If gmx_MMPBSA nevertheless returns an undefined-energy error, the engine scans
the retained thread-specific Sander mdout files and writes:

```text
undefined_energy_report.json
```

This report records the exact mdout filename, line, local energy record, and
context containing the overflow/NaN field.

## Engine 8 atom-mapping extensions

Engine 8.0.0 resolves three additional Desmond/VMD representations found in
the supplied receptor trajectories:

- **ACE/NME methyl hydrogens:** `HH31/HH32/HH33`,
  `1HH3/2HH3/3HH3`, `HA1/HA2/HA3`, and `1HA/2HA/3HA` are mapped to
  ff19SB `H1/H2/H3`.
- **NMA cap alias:** a Desmond `NMA` residue (`N`, `CA`, `H`, and three
  methyl hydrogens) is mapped to the ff19SB `NME` template (`N`, `C`, `H`,
  `H1/H2/H3`).
- **Exact ALA/GLY-like truncations:** when a residue label retains the
  sequence identity but its complete sampled atom set is an exact
  alanine- or glycine-like covalent graph, the engine may use that ff19SB
  surrogate rather than allow LEaP to invent an unsampled side chain. The
  substitution is restricted to an exact atom-count, element, name-mapping,
  and covalent-graph match and is recorded in
  `template_substitutions.tsv` and `preparation_report.json`.

For the supplied M4 structure, terminal `ARG427` contains only the backbone
and a methyl-capped C-beta (`HB1/HB2/HB3`), so its retained chemistry exactly
matches ALA after the carbonyl-H cap is replaced by NME. Engine 8 reports the
`ARG -> ALA` surrogate explicitly; it does not reconstruct a missing arginine
side chain.

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
/home/roshi/miniforge3/bin/conda
```

## Result-parser correction

Current gmx_MMPBSA output commonly labels the binding-energy row as Unicode `ΔTOTAL`. The earlier engine searched only for ASCII `DELTA TOTAL`, so a successful calculation could be followed by the false exception:

```text
ValueError: No DELTA TOTAL row found in FINAL_RESULTS_MMPBSA.dat
```

Engine 8.0.0 retains the 7.1.0 parser correction and accepts both `ΔTOTAL` and `DELTA TOTAL`. Summary parsing is also post-processing: if a future output format cannot be parsed, completed gmx_MMPBSA result files are retained and the calculation is not reported as failed.

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
tar -xzf desmond_gmxMMPBSA_wx_v9_1.tar.gz
cd desmond_gmxMMPBSA_wx_v9_1
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

For an existing installation at `/home/roshi/desmond_gmxMMPBSA`:

```bash
./install_wx_v9_1.sh /home/roshi/desmond_gmxMMPBSA
```

The installer creates timestamped backups, installs engine 9.1.0 and GUI 1.0.6, and adds:

```text
/home/roshi/desmond_gmxMMPBSA/desmond_gmxmmpbsa_wx.py
/home/roshi/desmond_gmxMMPBSA/run_desmond_gmxmmpbsa_wx.sh
```

Launch it with:

```bash
/home/roshi/desmond_gmxMMPBSA/run_desmond_gmxmmpbsa_wx.sh
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
