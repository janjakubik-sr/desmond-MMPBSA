#!/usr/bin/env python3
"""wxPython front end for the Desmond GRO/XTC -> gmx_MMPBSA workflow.

The GUI is deliberately a process launcher rather than a reimplementation of
molecular preparation.  It runs the bundled, validated command-line engine in
an isolated Conda environment and streams its output into a live log window.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
import unicodedata
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import remote_run

try:
    import wx
except ImportError as exc:  # pragma: no cover - depends on the desktop system
    print(
        "wxPython is required for the graphical interface.\n"
        "On Debian/Ubuntu install it with:\n\n"
        "    sudo apt install python3-wxgtk4.0\n\n"
        "Then launch this file with /usr/bin/python3, or set WX_PYTHON in the "
        "provided shell launcher.",
        file=sys.stderr,
    )
    raise SystemExit(1) from exc


GUI_VERSION = "1.3.0"
ENGINE_VERSION = "9.1.1"
APP_TITLE = "Desmond → gmx_MMPBSA"
CONFIG_DIR = Path.home() / ".config" / "desmond_gmxmmpbsa_wx"
SETTINGS_FILE = CONFIG_DIR / "settings.json"

AA_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN",
    "GLU", "GLH", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU",
    "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR",
    "VAL", "ACE", "NME",
}
COMMON_NONLIGAND = AA_NAMES | {
    "HOH", "WAT", "SOL", "SPC", "SPCE", "TIP3", "TIP3P", "TIP4P",
    "NA", "NA+", "SOD", "CL", "CL-", "CLA", "K", "K+", "POT",
    "CA", "MG", "ZN", "POPC", "POPE", "POPG", "DOPC", "DPPC",
    "CHOL", "CLR", "DMSO", "GOL", "PEG", "SO4", "PO4",
}

FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"

ENV_CHECK_SCRIPT = r'''
import importlib.util
import shutil
import sys

modules = ["numpy", "scipy", "mdtraj", "networkx", "parmed"]
programs = ["gmx", "gmx_MMPBSA", "parmchk2", "tleap", "antechamber", "mpirun"]
failed = False
print("Python:", sys.executable)
print("Version:", sys.version.replace("\\n", " "))
print("\\nPython modules")
for name in modules:
    ok = importlib.util.find_spec(name) is not None
    print(f"  {name:14s} {'OK' if ok else 'MISSING'}")
    failed |= not ok
print("\\nExecutables")
for name in programs:
    path = shutil.which(name)
    print(f"  {name:14s} {path or 'MISSING'}")
    failed |= path is None
raise SystemExit(2 if failed else 0)
'''


def parse_binding_energy(result_file: Path) -> dict[str, Any]:
    """Parse either ASCII ``DELTA TOTAL`` or Unicode ``ΔTOTAL`` output."""

    text = result_file.read_text(encoding="utf-8", errors="replace")
    normalized = unicodedata.normalize("NFKC", text).replace("\u2212", "-")
    frame_match = re.search(
        r"Calculations performed using\s+(\d+)\s+complex frames",
        normalized,
        flags=re.I,
    )

    matches: list[tuple[int, list[float], str]] = []
    lines = normalized.splitlines()
    for number, raw in enumerate(lines, start=1):
        line = raw.strip().lstrip("|+").strip()
        match = re.match(
            r"^(?:(?:Δ|DELTA)\s*[_ -]*TOTAL)\b(.*)$",
            line,
            flags=re.I,
        )
        if not match:
            continue
        values = [float(item) for item in re.findall(FLOAT_PATTERN, match.group(1))]
        if values:
            matches.append((number, values, raw))

    if not matches:
        candidates = [line for line in lines if "TOTAL" in line.upper() or "Δ" in line]
        preview = "\n".join(candidates[-12:])
        raise ValueError(
            f"No DELTA TOTAL/ΔTOTAL statistics row was found in {result_file}."
            + (f"\nCandidate lines:\n{preview}" if preview else "")
        )

    line_number, values, source_line = matches[-1]
    if len(values) >= 5:
        average, sd_prop, sd, sem_prop, sem = values[:5]
    elif len(values) == 3:
        average, sd, sem = values
        sd_prop = sem_prop = None
    elif len(values) == 2:
        average, sd = values
        sem = sd_prop = sem_prop = None
    else:
        raise ValueError(
            f"The matched ΔTOTAL row on line {line_number} contains only "
            f"{len(values)} numeric value: {source_line!r}"
        )

    prefix = "\n".join(lines[:line_number])
    pb = [
        m.start()
        for m in re.finditer(
            r"POISSON\s+BOLTZMANN|\bPB\s+CALCULATION\b", prefix, flags=re.I
        )
    ]
    gb = [
        m.start()
        for m in re.finditer(
            r"GENERALIZED\s+BORN|\bGB\s+CALCULATION\b", prefix, flags=re.I
        )
    ]
    model = "unknown"
    if pb or gb:
        model = "PB" if max(pb or [-1]) > max(gb or [-1]) else "GB"

    result: dict[str, Any] = {
        "result_file": str(result_file.resolve()),
        "solvent_model": model,
        "matched_label": "ΔTOTAL" if "Δ" in source_line else "DELTA TOTAL",
        "matched_line_number": line_number,
        "delta_total_kcal_per_mol": {
            "average": average,
            "standard_deviation_propagated": sd_prop,
            "standard_deviation": sd,
            "standard_error_propagated": sem_prop,
            "standard_error": sem,
        },
    }
    if frame_match:
        result["frames"] = int(frame_match.group(1))
    return result


def _resolve_executable_candidate(candidate: str) -> str | None:
    """Resolve an executable name or path without trusting PATH wrappers."""

    candidate = os.path.expandvars(os.path.expanduser(candidate.strip()))
    if not candidate:
        return None

    has_separator = os.sep in candidate or bool(os.altsep and os.altsep in candidate)
    resolved = candidate if has_separator else shutil.which(candidate)
    if not resolved:
        return None

    path = Path(resolved)
    if not path.is_file() or not os.access(path, os.X_OK):
        return None
    return str(path.resolve())


@lru_cache(maxsize=32)
def probe_conda(candidate: str) -> tuple[bool, str | None, str]:
    """Return whether *candidate* is a real Conda/Mamba CLI with ``run``.

    A shell helper named ``conda`` may implement only ``conda activate``.  Such
    helpers are usable interactively but cannot execute ``conda run`` and must
    not be selected by the GUI.
    """

    executable = _resolve_executable_candidate(candidate)
    if not executable:
        return False, None, "file not found or not executable"

    try:
        version = subprocess.run(
            [executable, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=8,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, executable, f"could not execute --version: {exc}"

    version_text = (version.stdout or "").strip()
    if version.returncode != 0 or not re.search(
        r"\b(?:conda|mamba)\b", version_text, flags=re.I
    ):
        detail = version_text or f"exit status {version.returncode}"
        return False, executable, f"not a Conda/Mamba CLI: {detail}"

    try:
        run_help = subprocess.run(
            [executable, "run", "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=8,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, executable, f"could not test 'run': {exc}"

    if run_help.returncode != 0:
        detail = (run_help.stdout or "").strip()
        return False, executable, (
            "does not support 'conda run': "
            + (detail or f"exit status {run_help.returncode}")
        )

    return True, executable, version_text


@lru_cache(maxsize=1)
def discover_conda() -> str:
    """Find a functional Conda/Miniforge executable with ``run`` support.

    Known Miniforge locations are checked before the user's PATH.  This avoids
    selecting interactive activation shims such as ``~/bin/conda``.
    """

    home = Path.home()
    candidates: list[str | None] = [
        os.environ.get("CONDA_EXE"),
        str(home / "miniforge3" / "bin" / "conda"),
        str(home / "miniforge3" / "condabin" / "conda"),
        str(home / "mambaforge" / "bin" / "conda"),
        str(home / "mambaforge" / "condabin" / "conda"),
        str(home / "miniconda3" / "bin" / "conda"),
        str(home / "miniconda3" / "condabin" / "conda"),
        str(home / "anaconda3" / "bin" / "conda"),
        str(home / "anaconda3" / "condabin" / "conda"),
        "/opt/conda/bin/conda",
        shutil.which("conda"),
        os.environ.get("MAMBA_EXE"),
        shutil.which("mamba"),
    ]

    seen: set[str] = set()
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        valid, executable, _detail = probe_conda(candidate)
        if valid and executable:
            return executable
    return ""


def resolve_conda_executable(requested: str | None) -> str:
    """Validate the selected CLI and fall back to a discovered Miniforge CLI."""

    requested = (requested or "").strip()
    requested_error = ""
    if requested:
        valid, executable, detail = probe_conda(requested)
        if valid and executable:
            return executable
        requested_error = f"Selected executable {requested!r}: {detail}."

    discovered = discover_conda()
    if discovered:
        return discovered

    message = (
        "No functional Conda/Miniforge executable with 'conda run' support "
        "was found. Select the real Miniforge binary, commonly:\n\n"
        "    ~/miniforge3/bin/conda\n\n"
        "An interactive shell helper that only performs 'conda activate' is "
        "not suitable."
    )
    if requested_error:
        message = requested_error + "\n\n" + message
    raise FileNotFoundError(message)


def format_command(command: Iterable[str]) -> str:
    return shlex.join([str(item) for item in command])


def open_external(path: Path) -> None:
    path = path.expanduser().resolve()
    if sys.platform.startswith("linux"):
        subprocess.Popen(["xdg-open", str(path)])
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    elif os.name == "nt":  # pragma: no cover
        os.startfile(str(path))  # type: ignore[attr-defined]
    else:  # pragma: no cover
        raise RuntimeError(f"No file opener is configured for {sys.platform}")


def read_mol2_info(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    title = ""
    atom_count: int | None = None
    charge_sum = 0.0

    try:
        molecule_index = next(
            i for i, line in enumerate(lines) if line.strip() == "@<TRIPOS>MOLECULE"
        )
        following = [line.strip() for line in lines[molecule_index + 1 :] if line.strip()]
        if following:
            title = following[0]
        if len(following) > 1:
            fields = following[1].split()
            if fields and fields[0].lstrip("+-").isdigit():
                atom_count = int(fields[0])
    except StopIteration:
        pass

    in_atoms = False
    parsed_atoms = 0
    for line in lines:
        stripped = line.strip()
        if stripped == "@<TRIPOS>ATOM":
            in_atoms = True
            continue
        if in_atoms and stripped.startswith("@<TRIPOS>"):
            break
        if in_atoms and stripped:
            fields = stripped.split()
            if len(fields) >= 9:
                parsed_atoms += 1
                try:
                    charge_sum += float(fields[8])
                except ValueError:
                    pass
    if parsed_atoms:
        atom_count = parsed_atoms
    return {"title": title, "atom_count": atom_count, "charge_sum": charge_sum}


def gro_resname_counts(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        handle.readline()
        count_line = handle.readline().strip()
        if not count_line.isdigit():
            raise ValueError(f"Invalid GRO atom-count line in {path}")
        atom_count = int(count_line)
        for _ in range(atom_count):
            line = handle.readline()
            if len(line) < 15:
                raise ValueError(f"Truncated GRO atom line in {path}")
            counts[line[5:10].strip().upper()] += 1
    return counts


def infer_ligand_resname(gro: Path, mol2: Path) -> tuple[str | None, str]:
    info = read_mol2_info(mol2)
    counts = gro_resname_counts(gro)
    atom_count = info.get("atom_count")
    title = str(info.get("title", "")).strip().upper()

    if title and title in counts and (atom_count is None or counts[title] == atom_count):
        return title, f"Matched MOL2 molecule name {title!r} in the GRO."

    candidates = [
        name
        for name, count in counts.items()
        if name not in COMMON_NONLIGAND and (atom_count is None or count == atom_count)
    ]
    if len(candidates) == 1:
        return candidates[0], (
            f"Matched the MOL2 atom count ({atom_count}) to GRO residue "
            f"{candidates[0]!r}."
        )
    if candidates:
        return None, "Ambiguous candidates: " + ", ".join(sorted(candidates))
    return None, "No unique GRO residue matched the MOL2 molecule name/atom count."


class PathField(wx.Panel):
    """Text field plus a file/directory browse button."""

    def __init__(
        self,
        parent: wx.Window,
        *,
        kind: str = "file",
        wildcard: str = "All files (*.*)|*.*",
        dialog_title: str = "Select path",
    ) -> None:
        super().__init__(parent)
        self.kind = kind
        self.wildcard = wildcard
        self.dialog_title = dialog_title
        self.text = wx.TextCtrl(self)
        self.button = wx.Button(self, label="Browse…")
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.text, 1, wx.EXPAND | wx.RIGHT, 6)
        sizer.Add(self.button, 0)
        self.SetSizer(sizer)
        self.button.Bind(wx.EVT_BUTTON, self.on_browse)

    def on_browse(self, _event: wx.CommandEvent) -> None:
        current = Path(self.text.GetValue()).expanduser() if self.text.GetValue() else None
        if self.kind == "dir":
            default = str(current if current and current.is_dir() else Path.home())
            dialog = wx.DirDialog(
                self,
                self.dialog_title,
                defaultPath=default,
                style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST,
            )
        else:
            default_dir = str(
                current.parent if current and current.parent.is_dir() else Path.home()
            )
            default_file = current.name if current else ""
            dialog = wx.FileDialog(
                self,
                self.dialog_title,
                defaultDir=default_dir,
                defaultFile=default_file,
                wildcard=self.wildcard,
                style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
            )
        with dialog:
            if dialog.ShowModal() == wx.ID_OK:
                self.text.SetValue(dialog.GetPath())

    def get(self) -> str:
        return self.text.GetValue().strip()

    def set(self, value: str | Path | None) -> None:
        self.text.SetValue("" if value is None else str(value))

    def enable(self, enabled: bool) -> None:
        self.text.Enable(enabled)
        self.button.Enable(enabled)


class MainFrame(wx.Frame):
    MODES = [
        ("Complete workflow: prepare + run", "all"),
        ("Prepare inputs only", "prepare"),
        ("Run an existing prepared directory", "run"),
        ("Summarize an existing result", "summarize"),
    ]

    def __init__(self, startup: argparse.Namespace) -> None:
        super().__init__(None, title=f"{APP_TITLE}  v{GUI_VERSION}", size=(1120, 850))
        self.SetMinSize((900, 680))
        self.engine_default = Path(__file__).resolve().parent / "app" / "desmond_gmxmmpbsa.py"
        self.ffxml_default = Path(__file__).resolve().parent / "app" / "protein.ff19SB.xml"
        self.processes: list[dict[str, Any]] = []
        self.process_seq = 0
        self.last_workdir: Path | None = None
        self.last_result: Path | None = None

        panel = wx.Panel(self)
        root = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(root)

        header = wx.BoxSizer(wx.HORIZONTAL)
        title = wx.StaticText(panel, label=APP_TITLE)
        title.SetFont(wx.Font(wx.FontInfo(16).Bold()))
        subtitle = wx.StaticText(
            panel,
            label=(
                f"ff19SB/GAFF2 single-trajectory rescoring · engine {ENGINE_VERSION}"
            ),
        )
        header.Add(title, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 14)
        header.Add(subtitle, 0, wx.ALIGN_CENTER_VERTICAL)
        header.AddStretchSpacer()
        root.Add(header, 0, wx.EXPAND | wx.ALL, 10)

        self.notebook = wx.Notebook(panel)
        self.calculation_panel = wx.Panel(self.notebook)
        self.advanced_panel = wx.Panel(self.notebook)
        self.results_panel = wx.Panel(self.notebook)
        self.notebook.AddPage(self.calculation_panel, "Calculation")
        self.notebook.AddPage(self.advanced_panel, "Advanced")
        self.notebook.AddPage(self.results_panel, "Results")
        root.Add(self.notebook, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 10)

        self._build_calculation_panel()
        self._build_advanced_panel()
        self._build_results_panel()

        progress_box = wx.BoxSizer(wx.HORIZONTAL)
        self.status_text = wx.StaticText(panel, label="Ready")
        self.progress = wx.Gauge(panel, range=100, style=wx.GA_HORIZONTAL)
        progress_box.Add(self.status_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        progress_box.Add(self.progress, 1, wx.ALIGN_CENTER_VERTICAL)
        root.Add(progress_box, 0, wx.EXPAND | wx.ALL, 10)

        button_box = wx.BoxSizer(wx.HORIZONTAL)
        self.run_button = wx.Button(panel, label="Run complete workflow")
        self.stop_button = wx.Button(panel, label="Stop")
        self.check_button = wx.Button(panel, label="Check environment")
        self.show_command_button = wx.Button(panel, label="Show command")
        self.clear_log_button = wx.Button(panel, label="Clear log")
        self.save_log_button = wx.Button(panel, label="Save log…")
        self.open_folder_button = wx.Button(panel, label="Open work directory")
        button_box.Add(self.run_button, 0, wx.RIGHT, 6)
        button_box.Add(self.stop_button, 0, wx.RIGHT, 6)
        button_box.Add(self.check_button, 0, wx.RIGHT, 6)
        button_box.Add(self.show_command_button, 0, wx.RIGHT, 6)
        button_box.AddStretchSpacer()
        button_box.Add(self.open_folder_button, 0, wx.RIGHT, 6)
        button_box.Add(self.clear_log_button, 0, wx.RIGHT, 6)
        button_box.Add(self.save_log_button, 0)
        root.Add(button_box, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)

        log_label = wx.StaticText(panel, label="Live log")
        log_label.SetFont(wx.Font(wx.FontInfo(10).Bold()))
        root.Add(log_label, 0, wx.LEFT | wx.RIGHT, 10)
        self.log = wx.TextCtrl(
            panel,
            style=wx.TE_MULTILINE | wx.TE_READONLY | wx.HSCROLL | wx.TE_RICH2,
        )
        self.log.SetFont(wx.Font(wx.FontInfo(9).Family(wx.FONTFAMILY_TELETYPE)))
        root.Add(self.log, 1, wx.EXPAND | wx.ALL, 10)

        self.CreateStatusBar()
        self.SetStatusText("Ready")
        self._build_menu()
        self.update_controls()
        self.open_folder_button.Disable()

        self.run_button.Bind(wx.EVT_BUTTON, self.on_run)
        self.run_remote.Bind(wx.EVT_CHECKBOX, lambda _evt: self.update_controls())
        self.stop_button.Bind(wx.EVT_BUTTON, self.on_stop)
        self.check_button.Bind(wx.EVT_BUTTON, self.on_check_environment)
        self.show_command_button.Bind(wx.EVT_BUTTON, self.on_show_command)
        self.clear_log_button.Bind(wx.EVT_BUTTON, lambda _evt: self.log.Clear())
        self.save_log_button.Bind(wx.EVT_BUTTON, self.on_save_log)
        self.open_folder_button.Bind(wx.EVT_BUTTON, self.on_open_workdir)
        self.mode_choice.Bind(wx.EVT_CHOICE, self.on_mode_changed)
        self.detect_ligand_button.Bind(wx.EVT_BUTTON, self.on_detect_ligand)
        self.parse_result_button.Bind(wx.EVT_BUTTON, self.on_parse_result)
        self.open_result_button.Bind(wx.EVT_BUTTON, self.on_open_result)
        self.result_file.text.Bind(wx.EVT_TEXT, self.on_result_path_changed)
        self.gro_file.text.Bind(wx.EVT_TEXT, self.on_input_path_changed)
        self.mol2_file.text.Bind(wx.EVT_TEXT, self.on_input_path_changed)
        self.workdir.text.Bind(wx.EVT_TEXT, self.on_workdir_changed)
        self.Bind(wx.EVT_CLOSE, self.on_close)

        self.load_settings(silent=True)
        if startup.engine:
            self.engine_file.set(startup.engine)
        if startup.profile:
            self.load_profile(Path(startup.profile))
        if startup.workdir:
            self.workdir.set(startup.workdir)
        self.on_mode_changed(None)
        self.Centre()

    # ------------------------------ UI construction -----------------------

    def _build_menu(self) -> None:
        menu_bar = wx.MenuBar()
        file_menu = wx.Menu()
        save_profile = file_menu.Append(wx.ID_SAVE, "Save profile…\tCtrl+S")
        load_profile = file_menu.Append(wx.ID_OPEN, "Load profile…\tCtrl+O")
        file_menu.AppendSeparator()
        exit_item = file_menu.Append(wx.ID_EXIT, "Exit")
        tools_menu = wx.Menu()
        check_item = tools_menu.Append(wx.ID_ANY, "Check environment")
        parse_item = tools_menu.Append(wx.ID_ANY, "Parse result file")
        help_menu = wx.Menu()
        about_item = help_menu.Append(wx.ID_ABOUT, "About")
        menu_bar.Append(file_menu, "File")
        menu_bar.Append(tools_menu, "Tools")
        menu_bar.Append(help_menu, "Help")
        self.SetMenuBar(menu_bar)
        self.Bind(wx.EVT_MENU, self.on_save_profile, save_profile)
        self.Bind(wx.EVT_MENU, self.on_load_profile, load_profile)
        self.Bind(wx.EVT_MENU, lambda _evt: self.Close(), exit_item)
        self.Bind(wx.EVT_MENU, self.on_check_environment, check_item)
        self.Bind(wx.EVT_MENU, self.on_parse_result, parse_item)
        self.Bind(wx.EVT_MENU, self.on_about, about_item)

    def _build_calculation_panel(self) -> None:
        panel = self.calculation_panel
        grid = wx.FlexGridSizer(cols=2, vgap=7, hgap=10)
        grid.AddGrowableCol(1, 1)

        self.mode_choice = wx.Choice(panel, choices=[label for label, _ in self.MODES])
        self.mode_choice.SetSelection(0)
        self.run_remote = wx.CheckBox(
            panel, label="Run on a remote server via Slurm (copy to scratch, sbatch, copy back)"
        )
        self.workdir = PathField(panel, kind="dir", dialog_title="Select work/output directory")
        self.gro_file = PathField(
            panel,
            wildcard="GROMACS structure (*.gro)|*.gro|All files (*.*)|*.*",
            dialog_title="Select full-system GRO",
        )
        self.xtc_file = PathField(
            panel,
            wildcard="GROMACS trajectory (*.xtc)|*.xtc|All files (*.*)|*.*",
            dialog_title="Select full-system XTC",
        )
        self.mol2_file = PathField(
            panel,
            wildcard="Tripos MOL2 (*.mol2)|*.mol2|All files (*.*)|*.*",
            dialog_title="Select ligand MOL2",
        )

        ligand_row = wx.Panel(panel)
        ligand_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ligand_resname = wx.TextCtrl(ligand_row, size=(110, -1))
        self.detect_ligand_button = wx.Button(ligand_row, label="Detect from GRO/MOL2")
        ligand_sizer.Add(self.ligand_resname, 0, wx.RIGHT, 7)
        ligand_sizer.Add(self.detect_ligand_button, 0)
        ligand_row.SetSizer(ligand_sizer)

        parameter_row = wx.Panel(panel)
        parameter_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ligand_charge = wx.TextCtrl(parameter_row, size=(70, -1))
        self.ligand_charge.SetHint("infer")
        self.stride = wx.SpinCtrl(parameter_row, min=1, max=100000, initial=10, size=(90, -1))
        self.solvent_model = wx.Choice(parameter_row, choices=["PB", "GB"])
        self.solvent_model.SetSelection(0)
        self.mpi_processes = wx.SpinCtrl(
            parameter_row,
            min=1,
            max=max(256, os.cpu_count() or 1),
            initial=max(1, min(8, os.cpu_count() or 1)),
            size=(90, -1),
        )
        for label, control in (
            ("Charge", self.ligand_charge),
            ("Stride", self.stride),
            ("Model", self.solvent_model),
            ("MPI", self.mpi_processes),
        ):
            parameter_sizer.Add(wx.StaticText(parameter_row, label=label), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
            parameter_sizer.Add(control, 0, wx.RIGHT, 12)
        parameter_row.SetSizer(parameter_sizer)

        rows = [
            ("Operation", self.mode_choice),
            ("Remote execution", self.run_remote),
            ("Work/output directory", self.workdir),
            ("Full-system GRO", self.gro_file),
            ("Full-system XTC", self.xtc_file),
            ("Ligand MOL2", self.mol2_file),
            ("Ligand residue name", ligand_row),
            ("Calculation parameters", parameter_row),
        ]
        for label, control in rows:
            text = wx.StaticText(panel, label=label)
            text.SetMinSize((175, -1))
            grid.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 3)
            grid.Add(control, 1, wx.EXPAND)

        note = wx.StaticText(
            panel,
            label=(
                "The output directory may be typed even if it does not yet exist. "
                "For a prepared-directory run, select the directory containing "
                "preparation_report.json."
            ),
        )
        note.Wrap(850)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(grid, 0, wx.EXPAND | wx.ALL, 12)
        sizer.Add(note, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)
        panel.SetSizer(sizer)

    def _build_advanced_panel(self) -> None:
        panel = self.advanced_panel
        grid = wx.FlexGridSizer(cols=2, vgap=7, hgap=10)
        grid.AddGrowableCol(1, 1)

        self.chunk_size = wx.SpinCtrl(panel, min=1, max=10000, initial=20, size=(100, -1))
        self.preserve_ligand_charges = wx.CheckBox(
            panel, label="Preserve supplied MOL2 charges exactly"
        )
        self.cap_mode = wx.Choice(panel, choices=["ace-nme", "error"])
        self.cap_mode.SetSelection(0)
        self.internal_h_mode = wx.Choice(panel, choices=["planar-remove", "error"])
        self.internal_h_mode.SetSelection(0)
        self.ffxml_file = PathField(
            panel,
            wildcard="OpenMM force-field XML (*.xml)|*.xml|All files (*.*)|*.*",
            dialog_title="Select ff19SB XML",
        )
        self.ffxml_file.set(self.ffxml_default)
        self.conda_file = PathField(
            panel,
            wildcard="Executable files (*)|*|All files (*.*)|*.*",
            dialog_title="Select Conda executable",
        )
        self.conda_file.set(discover_conda())
        self.conda_env = wx.TextCtrl(panel, value="desmond-mmpbsa")
        self.use_current_python = wx.CheckBox(
            panel,
            label="Use the Python/environment that launched this GUI instead of conda run",
        )
        self.engine_file = PathField(
            panel,
            wildcard="Python program (*.py)|*.py|All files (*.*)|*.*",
            dialog_title="Select workflow engine",
        )
        self.engine_file.set(self.engine_default)
        self.extra_arguments = wx.TextCtrl(panel)
        self.extra_arguments.SetHint("Optional expert CLI arguments")

        self.remote_host = wx.TextCtrl(panel, value=remote_run.HOST_DEFAULT)
        self.remote_user = wx.TextCtrl(panel)
        self.remote_user.SetHint("optional; empty = SSH config default")
        self.remote_partition = wx.TextCtrl(panel, value=remote_run.PARTITION_DEFAULT)
        self.remote_walltime = wx.TextCtrl(panel, value=remote_run.WALLTIME_DEFAULT)
        self.remote_scratch = wx.TextCtrl(panel, value=remote_run.SCRATCH_BASE_DEFAULT)
        self.remote_conda = wx.TextCtrl(panel, value=remote_run.REMOTE_CONDA_DEFAULT)
        self.remote_stage_rsync = wx.CheckBox(
            panel,
            label="Copy inputs/results via rsync over SSH (no shared Data mount)",
        )

        rows = [
            ("MDTraj chunk size", self.chunk_size),
            ("Ligand charges", self.preserve_ligand_charges),
            ("Hydrogen-cap handling", self.cap_mode),
            ("Internal duplicate N-H", self.internal_h_mode),
            ("Protein ff19SB XML", self.ffxml_file),
            ("Conda/Miniforge executable", self.conda_file),
            ("Conda environment", self.conda_env),
            ("Python execution", self.use_current_python),
            ("Workflow engine", self.engine_file),
            ("Additional arguments", self.extra_arguments),
            ("Remote SSH host", self.remote_host),
            ("Remote SSH user", self.remote_user),
            ("Remote Slurm partition", self.remote_partition),
            ("Remote Slurm walltime", self.remote_walltime),
            ("Remote scratch base", self.remote_scratch),
            ("Remote conda executable", self.remote_conda),
            ("Remote staging", self.remote_stage_rsync),
        ]
        for label, control in rows:
            text = wx.StaticText(panel, label=label)
            text.SetMinSize((190, -1))
            grid.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 3)
            grid.Add(control, 1, wx.EXPAND)

        explanation = wx.StaticText(
            panel,
            label=(
                "Recommended Linux setup: launch this GUI with system Python/wxPython, "
                "then execute the scientific engine through 'conda run -n "
                "desmond-mmpbsa'. This keeps wxPython separate from the molecular-"
                "modelling environment."
            ),
        )
        explanation.Wrap(850)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(grid, 0, wx.EXPAND | wx.ALL, 12)
        sizer.Add(explanation, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)
        panel.SetSizer(sizer)

    def _build_results_panel(self) -> None:
        panel = self.results_panel
        grid = wx.FlexGridSizer(cols=2, vgap=7, hgap=10)
        grid.AddGrowableCol(1, 1)
        self.result_file = PathField(
            panel,
            wildcard="gmx_MMPBSA result (*.dat)|*.dat|All files (*.*)|*.*",
            dialog_title="Select FINAL_RESULTS_MMPBSA.dat",
        )
        result_actions = wx.Panel(panel)
        self.parse_result_button = wx.Button(result_actions, label="Parse result")
        self.open_result_button = wx.Button(result_actions, label="Open result file")
        action_sizer = wx.BoxSizer(wx.HORIZONTAL)
        action_sizer.Add(self.parse_result_button, 0, wx.RIGHT, 7)
        action_sizer.Add(self.open_result_button, 0)
        result_actions.SetSizer(action_sizer)

        self.summary_fields: dict[str, wx.TextCtrl] = {}
        labels = [
            ("model", "Solvent model"),
            ("frames", "Frames"),
            ("average", "ΔTOTAL average (kcal/mol)"),
            ("sd_prop", "SD propagated"),
            ("sd", "SD"),
            ("sem_prop", "SEM propagated"),
            ("sem", "SEM"),
            ("label", "Matched result label"),
        ]
        grid.Add(wx.StaticText(panel, label="Result DAT"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.result_file, 1, wx.EXPAND)
        grid.Add(wx.StaticText(panel, label="Actions"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(result_actions, 0, wx.EXPAND)
        for key, label in labels:
            control = wx.TextCtrl(panel, style=wx.TE_READONLY)
            self.summary_fields[key] = control
            grid.Add(wx.StaticText(panel, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(control, 1, wx.EXPAND)

        note = wx.StaticText(
            panel,
            label=(
                "The parser accepts both the older ASCII label 'DELTA TOTAL' and "
                "the Unicode label 'ΔTOTAL' used by current gmx_MMPBSA output."
            ),
        )
        note.Wrap(850)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(grid, 0, wx.EXPAND | wx.ALL, 12)
        sizer.Add(note, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)
        panel.SetSizer(sizer)
        self.open_result_button.Disable()

    # ------------------------------ settings ------------------------------

    def current_mode(self) -> str:
        selection = self.mode_choice.GetSelection()
        return self.MODES[selection][1]

    def remote_active(self) -> bool:
        return self.run_remote.GetValue() and self.current_mode() in {"all", "run"}

    def local_workflow_running(self) -> bool:
        return any(
            p["task"] == "workflow" and not p["is_remote"] for p in self.processes
        )

    def update_controls(self) -> None:
        running = bool(self.processes)
        if running:
            # A further run is allowed only for remote submissions while
            # no local workflow occupies this machine.
            allowed = self.remote_active() and not self.local_workflow_running()
        else:
            allowed = True
        self.run_button.Enable(allowed)
        self.check_button.Enable(not running)
        self.stop_button.Enable(running)

    def settings_dict(self) -> dict[str, Any]:
        return {
            "mode": self.current_mode(),
            "workdir": self.workdir.get(),
            "gro": self.gro_file.get(),
            "xtc": self.xtc_file.get(),
            "mol2": self.mol2_file.get(),
            "ligand_resname": self.ligand_resname.GetValue().strip(),
            "ligand_charge": self.ligand_charge.GetValue().strip(),
            "stride": self.stride.GetValue(),
            "solvent_model": self.solvent_model.GetStringSelection().lower(),
            "np": self.mpi_processes.GetValue(),
            "chunk": self.chunk_size.GetValue(),
            "preserve_ligand_charges": self.preserve_ligand_charges.GetValue(),
            "cap_mode": self.cap_mode.GetStringSelection(),
            "internal_h_mode": self.internal_h_mode.GetStringSelection(),
            "ffxml": self.ffxml_file.get(),
            "conda": self.conda_file.get(),
            "conda_env": self.conda_env.GetValue().strip(),
            "use_current_python": self.use_current_python.GetValue(),
            "engine": self.engine_file.get(),
            "extra_arguments": self.extra_arguments.GetValue(),
            "result_file": self.result_file.get(),
            "run_remote": self.run_remote.GetValue(),
            "remote_host": self.remote_host.GetValue().strip(),
            "remote_user": self.remote_user.GetValue().strip(),
            "remote_partition": self.remote_partition.GetValue().strip(),
            "remote_walltime": self.remote_walltime.GetValue().strip(),
            "remote_scratch": self.remote_scratch.GetValue().strip(),
            "remote_conda": self.remote_conda.GetValue().strip(),
            "remote_rsync_stage": self.remote_stage_rsync.GetValue(),
        }

    def apply_settings(self, data: dict[str, Any]) -> None:
        mode = str(data.get("mode", "all"))
        index = next((i for i, (_, value) in enumerate(self.MODES) if value == mode), 0)
        self.mode_choice.SetSelection(index)
        self.workdir.set(data.get("workdir", ""))
        self.gro_file.set(data.get("gro", ""))
        self.xtc_file.set(data.get("xtc", ""))
        self.mol2_file.set(data.get("mol2", ""))
        self.ligand_resname.SetValue(str(data.get("ligand_resname", "")))
        self.ligand_charge.SetValue(str(data.get("ligand_charge", "")))
        self.stride.SetValue(int(data.get("stride", 10)))
        model = str(data.get("solvent_model", "pb")).upper()
        self.solvent_model.SetSelection(0 if model == "PB" else 1)
        self.mpi_processes.SetValue(int(data.get("np", max(1, min(8, os.cpu_count() or 1)))))
        self.chunk_size.SetValue(int(data.get("chunk", 20)))
        self.preserve_ligand_charges.SetValue(bool(data.get("preserve_ligand_charges", False)))
        cap_mode = str(data.get("cap_mode", "ace-nme"))
        self.cap_mode.SetSelection(max(0, self.cap_mode.FindString(cap_mode)))
        internal_mode = str(data.get("internal_h_mode", "planar-remove"))
        self.internal_h_mode.SetSelection(max(0, self.internal_h_mode.FindString(internal_mode)))
        self.ffxml_file.set(data.get("ffxml", self.ffxml_default))
        saved_conda = str(data.get("conda", "")).strip()
        try:
            selected_conda = resolve_conda_executable(saved_conda)
        except FileNotFoundError:
            selected_conda = saved_conda or discover_conda()
        self.conda_file.set(selected_conda)
        self.conda_env.SetValue(str(data.get("conda_env", "desmond-mmpbsa")))
        self.use_current_python.SetValue(bool(data.get("use_current_python", False)))
        self.engine_file.set(data.get("engine", self.engine_default))
        self.extra_arguments.SetValue(str(data.get("extra_arguments", "")))
        self.result_file.set(data.get("result_file", ""))
        self.run_remote.SetValue(bool(data.get("run_remote", False)))
        self.remote_host.SetValue(
            str(data.get("remote_host", remote_run.HOST_DEFAULT)))
        self.remote_user.SetValue(str(data.get("remote_user", "")))
        self.remote_partition.SetValue(
            str(data.get("remote_partition", remote_run.PARTITION_DEFAULT)))
        self.remote_walltime.SetValue(
            str(data.get("remote_walltime", remote_run.WALLTIME_DEFAULT)))
        self.remote_scratch.SetValue(
            str(data.get("remote_scratch", remote_run.SCRATCH_BASE_DEFAULT)))
        self.remote_conda.SetValue(
            str(data.get("remote_conda", remote_run.REMOTE_CONDA_DEFAULT)))
        self.remote_stage_rsync.SetValue(
            bool(data.get("remote_rsync_stage", False)))
        self.on_mode_changed(None)

    def save_settings(self) -> None:
        CONFIG_DIR.mkdir(parents=True, exist_ok=True)
        SETTINGS_FILE.write_text(
            json.dumps(self.settings_dict(), indent=2) + "\n", encoding="utf-8"
        )

    def load_settings(self, *, silent: bool = False) -> None:
        if not SETTINGS_FILE.is_file():
            return
        try:
            self.apply_settings(json.loads(SETTINGS_FILE.read_text(encoding="utf-8")))
        except Exception as exc:
            if not silent:
                wx.MessageBox(str(exc), "Could not load settings", wx.OK | wx.ICON_WARNING)

    def save_profile(self, path: Path) -> None:
        path.write_text(json.dumps(self.settings_dict(), indent=2) + "\n", encoding="utf-8")

    def load_profile(self, path: Path) -> None:
        self.apply_settings(json.loads(path.read_text(encoding="utf-8")))

    # ------------------------------ events --------------------------------

    def on_mode_changed(self, _event: wx.CommandEvent | None) -> None:
        mode = self.current_mode()
        prepare_inputs = mode in {"all", "prepare"}
        self.gro_file.enable(prepare_inputs)
        self.xtc_file.enable(prepare_inputs)
        self.mol2_file.enable(prepare_inputs)
        self.ligand_resname.Enable(prepare_inputs)
        self.detect_ligand_button.Enable(prepare_inputs)
        self.ligand_charge.Enable(prepare_inputs)
        self.stride.Enable(prepare_inputs)
        self.solvent_model.Enable(prepare_inputs)
        self.chunk_size.Enable(prepare_inputs)
        self.preserve_ligand_charges.Enable(prepare_inputs)
        self.cap_mode.Enable(prepare_inputs)
        self.internal_h_mode.Enable(prepare_inputs)
        self.ffxml_file.enable(prepare_inputs)
        self.mpi_processes.Enable(mode in {"all", "run"})
        self.workdir.enable(mode != "summarize")
        self.run_remote.Enable(mode in {"all", "run"})
        if mode == "all":
            self.run_button.SetLabel("Run complete workflow")
        elif mode == "prepare":
            self.run_button.SetLabel("Prepare inputs")
        elif mode == "run":
            self.run_button.SetLabel("Run gmx_MMPBSA")
        else:
            self.run_button.SetLabel("Summarize result")
            self.notebook.SetSelection(2)
        self.update_controls()
        self.Layout()

    def on_input_path_changed(self, _event: wx.CommandEvent) -> None:
        gro_value = self.gro_file.get()
        if gro_value:
            gro = Path(gro_value).expanduser()
            if gro.is_file():
                if not self.xtc_file.get():
                    candidate = gro.with_suffix(".xtc")
                    if candidate.is_file():
                        self.xtc_file.set(candidate)
                if not self.workdir.get():
                    self.workdir.set(gro.parent / "LIG_mmpbsa")
        if self.gro_file.get() and self.mol2_file.get() and not self.ligand_resname.GetValue().strip():
            try:
                inferred, _message = infer_ligand_resname(
                    Path(self.gro_file.get()).expanduser(),
                    Path(self.mol2_file.get()).expanduser(),
                )
                if inferred:
                    self.ligand_resname.SetValue(inferred)
            except Exception:
                pass
        mol2_value = self.mol2_file.get()
        if mol2_value and not self.ligand_charge.GetValue().strip():
            try:
                info = read_mol2_info(Path(mol2_value).expanduser())
                self.ligand_charge.SetValue(str(round(float(info["charge_sum"]))))
            except Exception:
                pass

    def on_workdir_changed(self, _event: wx.CommandEvent) -> None:
        value = self.workdir.get()
        if value:
            result = Path(value).expanduser() / "FINAL_RESULTS_MMPBSA.dat"
            self.result_file.set(result)
            self.last_workdir = Path(value).expanduser()
            self.open_folder_button.Enable(True)

    def on_result_path_changed(self, _event: wx.CommandEvent) -> None:
        value = self.result_file.get()
        self.open_result_button.Enable(bool(value and Path(value).expanduser().is_file()))

    def on_detect_ligand(self, _event: wx.CommandEvent) -> None:
        try:
            gro = Path(self.gro_file.get()).expanduser().resolve()
            mol2 = Path(self.mol2_file.get()).expanduser().resolve()
            if not gro.is_file() or not mol2.is_file():
                raise FileNotFoundError("Select existing GRO and MOL2 files first.")
            inferred, message = infer_ligand_resname(gro, mol2)
            if inferred:
                self.ligand_resname.SetValue(inferred)
                info = read_mol2_info(mol2)
                if not self.ligand_charge.GetValue().strip():
                    self.ligand_charge.SetValue(str(round(float(info["charge_sum"]))))
                wx.MessageBox(message, "Ligand detected", wx.OK | wx.ICON_INFORMATION)
            else:
                wx.MessageBox(message, "Ligand detection", wx.OK | wx.ICON_WARNING)
        except Exception as exc:
            wx.MessageBox(str(exc), "Ligand detection failed", wx.OK | wx.ICON_ERROR)

    def on_run(self, _event: wx.CommandEvent) -> None:
        if self.processes and not (
            self.remote_active() and not self.local_workflow_running()
        ):
            wx.MessageBox(
                "Another process is already running. Only remote "
                "submissions may run concurrently.",
                APP_TITLE,
                wx.OK | wx.ICON_WARNING,
            )
            return
        mode = self.current_mode()
        if mode == "summarize":
            self.on_parse_result(None)
            return
        try:
            command, cwd, workdir = self.build_engine_command()
        except Exception as exc:
            wx.MessageBox(str(exc), "Invalid calculation setup", wx.OK | wx.ICON_ERROR)
            return

        if mode in {"all", "prepare"} and workdir.exists() and any(workdir.iterdir()):
            answer = wx.MessageBox(
                f"The output directory is not empty:\n\n{workdir}\n\n"
                "Prepared files and stale result files may be replaced. Continue?",
                "Existing output directory",
                wx.YES_NO | wx.NO_DEFAULT | wx.ICON_WARNING,
            )
            if answer != wx.YES:
                return

        workdir.mkdir(parents=True, exist_ok=True)
        self.last_workdir = workdir
        self.result_file.set(workdir / "FINAL_RESULTS_MMPBSA.dat")
        is_remote = self.remote_active()
        if is_remote:
            try:
                command = self.build_remote_command(command, workdir, mode)
                cwd = workdir
            except Exception as exc:
                wx.MessageBox(
                    str(exc), "Remote submission setup failed", wx.OK | wx.ICON_ERROR
                )
                return
        try:
            self.start_process(
                command,
                cwd,
                "workflow",
                workdir / "wx_gui_run.log",
                workdir=workdir,
                is_remote=is_remote,
            )
        except Exception as exc:
            wx.MessageBox(str(exc), "Could not start workflow", wx.OK | wx.ICON_ERROR)

    def on_check_environment(self, _event: wx.CommandEvent | None) -> None:
        if self.processes:
            wx.MessageBox("A process is already running.", APP_TITLE, wx.OK | wx.ICON_WARNING)
            return
        try:
            if self.remote_active():
                command = self.remote_environment_check_command()
            else:
                command = self.environment_check_command()
        except Exception as exc:
            wx.MessageBox(str(exc), "Environment check", wx.OK | wx.ICON_ERROR)
            return
        try:
            self.start_process(command, Path.home(), "environment", None)
        except Exception as exc:
            wx.MessageBox(str(exc), "Could not start environment check", wx.OK | wx.ICON_ERROR)

    def on_show_command(self, _event: wx.CommandEvent) -> None:
        try:
            if self.current_mode() == "summarize":
                text = f"Internal parser: {self.result_file.get()}"
            else:
                command, cwd, workdir = self.build_engine_command(validate_paths=False)
                if self.remote_active():
                    text = self.remote_summary_text(command, cwd, workdir)
                else:
                    text = f"Working directory:\n{cwd}\n\nCommand:\n{format_command(command)}"
            dialog = wx.MessageDialog(self, text, "Command preview", wx.OK | wx.ICON_INFORMATION)
            dialog.SetSize((900, 400))
            with dialog:
                dialog.ShowModal()
        except Exception as exc:
            wx.MessageBox(str(exc), "Command preview", wx.OK | wx.ICON_ERROR)

    def on_stop(self, _event: wx.CommandEvent) -> None:
        if not self.processes:
            return
        self.append_log("\n[GUI] Termination requested for all running processes.\n")
        for record in list(self.processes):
            record["stop_requested"] = True
            self.cancel_remote_job(record)
            popen = record["popen"]
            if popen.poll() is not None:
                continue
            try:
                if os.name == "posix":
                    os.killpg(os.getpgid(popen.pid), signal.SIGTERM)
                else:  # pragma: no cover
                    popen.terminate()
            except ProcessLookupError:
                pass
        wx.CallLater(5000, self.force_kill_if_running)

    def cancel_remote_job(self, record: dict[str, Any]) -> None:
        job_id = record.get("remote_job_id")
        if not job_id:
            return
        host = self.remote_host.GetValue().strip()
        self.append_log(f"[GUI] Cancelling Slurm job {job_id} on {host} ...\n")

        def _cancel() -> None:
            try:
                subprocess.run(
                    [
                        "ssh",
                        "-o",
                        "BatchMode=yes",
                        "-o",
                        "ConnectTimeout=15",
                        host,
                        "scancel",
                        job_id,
                    ],
                    timeout=30,
                    capture_output=True,
                )
            except Exception:
                pass

        threading.Thread(target=_cancel, daemon=True).start()

    def force_kill_if_running(self) -> None:
        for record in list(self.processes):
            popen = record["popen"]
            if popen.poll() is not None:
                continue
            self.append_log(
                f"[GUI] Process [#{record['seq']}] did not terminate; "
                "sending SIGKILL.\n"
            )
            try:
                if os.name == "posix":
                    os.killpg(os.getpgid(popen.pid), signal.SIGKILL)
                else:  # pragma: no cover
                    popen.kill()
            except ProcessLookupError:
                pass

    def on_parse_result(self, _event: wx.CommandEvent | None) -> None:
        try:
            path = Path(self.result_file.get()).expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(f"Result file not found: {path}")
            summary = parse_binding_energy(path)
            summary_path = path.with_name("binding_energy_summary.json")
            summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
            self.display_summary(summary)
            self.last_result = path
            self.open_result_button.Enable(True)
            self.notebook.SetSelection(2)
            self.SetStatusText(f"Parsed {path.name}")
            self.append_log(
                "\n[GUI] Result summary parsed successfully:\n"
                + json.dumps(summary, indent=2)
                + "\n"
            )
        except Exception as exc:
            wx.MessageBox(str(exc), "Result parsing failed", wx.OK | wx.ICON_ERROR)

    def on_open_workdir(self, _event: wx.CommandEvent) -> None:
        value = self.workdir.get()
        if not value:
            return
        path = Path(value).expanduser()
        if not path.exists():
            wx.MessageBox(f"Directory does not exist: {path}", APP_TITLE, wx.OK | wx.ICON_WARNING)
            return
        open_external(path)

    def on_open_result(self, _event: wx.CommandEvent) -> None:
        path = Path(self.result_file.get()).expanduser()
        if path.is_file():
            open_external(path)

    def on_save_log(self, _event: wx.CommandEvent) -> None:
        with wx.FileDialog(
            self,
            "Save GUI log",
            defaultFile="desmond_gmxMMPBSA_wx.log",
            wildcard="Log files (*.log)|*.log|Text files (*.txt)|*.txt|All files (*.*)|*.*",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        ) as dialog:
            if dialog.ShowModal() == wx.ID_OK:
                Path(dialog.GetPath()).write_text(self.log.GetValue(), encoding="utf-8")

    def on_save_profile(self, _event: wx.CommandEvent) -> None:
        with wx.FileDialog(
            self,
            "Save calculation profile",
            defaultFile="gmxMMPBSA_profile.json",
            wildcard="JSON profile (*.json)|*.json|All files (*.*)|*.*",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        ) as dialog:
            if dialog.ShowModal() == wx.ID_OK:
                self.save_profile(Path(dialog.GetPath()))

    def on_load_profile(self, _event: wx.CommandEvent) -> None:
        with wx.FileDialog(
            self,
            "Load calculation profile",
            wildcard="JSON profile (*.json)|*.json|All files (*.*)|*.*",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
        ) as dialog:
            if dialog.ShowModal() == wx.ID_OK:
                try:
                    self.load_profile(Path(dialog.GetPath()))
                except Exception as exc:
                    wx.MessageBox(str(exc), "Profile load failed", wx.OK | wx.ICON_ERROR)

    def on_about(self, _event: wx.CommandEvent) -> None:
        wx.MessageBox(
            f"{APP_TITLE}\n\n"
            f"GUI version: {GUI_VERSION}\n"
            f"Bundled workflow engine: {ENGINE_VERSION}\n\n"
            "Runs the validated Desmond GRO/XTC → ff19SB/GAFF2 → "
            "gmx_MMPBSA single-trajectory workflow in a Conda environment.\n"
            "Calculations can optionally be submitted to a remote server "
            "via Slurm (staging to /scratch and copy-back to the work "
            "directory are handled automatically, either through the "
            "shared Data mount or via rsync over SSH).\n\n"
            "The GUI itself can use Debian/Ubuntu's system wxPython; the "
            "scientific calculation remains isolated in the selected Conda environment.",
            "About",
            wx.OK | wx.ICON_INFORMATION,
        )

    def on_close(self, event: wx.CloseEvent) -> None:
        if self.processes:
            answer = wx.MessageBox(
                "One or more calculations are still running. Stop them and "
                "close the application?",
                APP_TITLE,
                wx.YES_NO | wx.NO_DEFAULT | wx.ICON_WARNING,
            )
            if answer != wx.YES:
                event.Veto()
                return
            self.on_stop(wx.CommandEvent())
        try:
            self.save_settings()
        except Exception:
            pass
        event.Skip()

    # ------------------------------ command construction ------------------

    def python_prefix(self) -> list[str]:
        engine = Path(self.engine_file.get()).expanduser().resolve()
        if not engine.is_file():
            raise FileNotFoundError(f"Workflow engine not found: {engine}")
        if self.use_current_python.GetValue():
            return [sys.executable, "-u", str(engine)]
        conda_text = self.conda_file.get()
        conda_command = resolve_conda_executable(conda_text)
        if conda_command != conda_text:
            self.conda_file.set(conda_command)
        environment = self.conda_env.GetValue().strip()
        if not environment:
            raise ValueError("Specify the Conda environment name.")
        return [
            conda_command,
            "run",
            "--no-capture-output",
            "-n",
            environment,
            "python",
            "-u",
            str(engine),
        ]

    def build_engine_command(
        self, *, validate_paths: bool = True
    ) -> tuple[list[str], Path, Path]:
        mode = self.current_mode()
        if mode == "summarize":
            raise ValueError("Summarization is handled directly by the GUI.")

        workdir_text = self.workdir.get()
        if not workdir_text:
            raise ValueError("Select or type a work/output directory.")
        workdir = Path(workdir_text).expanduser().resolve()
        if workdir.exists() and not workdir.is_dir():
            raise ValueError(f"Work/output path is not a directory: {workdir}")
        command = self.python_prefix() + [mode]

        if mode in {"all", "prepare"}:
            paths = {
                "GRO": Path(self.gro_file.get()).expanduser().resolve(),
                "XTC": Path(self.xtc_file.get()).expanduser().resolve(),
                "MOL2": Path(self.mol2_file.get()).expanduser().resolve(),
                "ff19SB XML": Path(self.ffxml_file.get()).expanduser().resolve(),
            }
            if validate_paths:
                for label, path in paths.items():
                    if not path.is_file():
                        raise FileNotFoundError(f"{label} file not found: {path}")
            ligand = self.ligand_resname.GetValue().strip()
            if not ligand:
                raise ValueError("Specify the ligand residue name used in the GRO.")
            if not re.fullmatch(r"[A-Za-z0-9_+\-]{1,8}", ligand):
                raise ValueError(
                    "Ligand residue name may contain letters, digits, underscore, + or - only."
                )
            charge_text = self.ligand_charge.GetValue().strip()
            if charge_text:
                try:
                    int(charge_text)
                except ValueError as exc:
                    raise ValueError("Ligand charge must be an integer or blank.") from exc

            command += [
                "--gro", str(paths["GRO"]),
                "--xtc", str(paths["XTC"]),
                "--mol2", str(paths["MOL2"]),
                "--ligand-resname", ligand,
                "--outdir", str(workdir),
                "--ffxml", str(paths["ff19SB XML"]),
                "--stride", str(self.stride.GetValue()),
                "--chunk", str(self.chunk_size.GetValue()),
                "--solvent-model", self.solvent_model.GetStringSelection().lower(),
                "--hydrogen-cap-mode", self.cap_mode.GetStringSelection(),
                "--internal-duplicate-h-mode", self.internal_h_mode.GetStringSelection(),
            ]
            if charge_text:
                command += ["--ligand-charge", charge_text]
            if self.preserve_ligand_charges.GetValue():
                command.append("--preserve-ligand-charges")
            if mode == "all":
                command += ["--np", str(self.mpi_processes.GetValue())]
            cwd = paths["GRO"].parent if paths["GRO"].parent.exists() else Path.home()
        else:
            if validate_paths:
                if not workdir.is_dir():
                    raise FileNotFoundError(f"Prepared work directory not found: {workdir}")
                manifest = workdir / "preparation_report.json"
                if not manifest.is_file():
                    raise FileNotFoundError(
                        f"Prepared-run manifest not found: {manifest}"
                    )
            command += ["--workdir", str(workdir), "--np", str(self.mpi_processes.GetValue())]
            cwd = workdir if workdir.exists() else workdir.parent

        extra = self.extra_arguments.GetValue().strip()
        if extra:
            command += shlex.split(extra)
        return command, cwd, workdir

    def environment_check_command(self) -> list[str]:
        script = ENV_CHECK_SCRIPT
        if self.use_current_python.GetValue():
            return [sys.executable, "-u", "-c", script]
        conda_text = self.conda_file.get()
        conda_command = resolve_conda_executable(conda_text)
        if conda_command != conda_text:
            self.conda_file.set(conda_command)
        environment = self.conda_env.GetValue().strip()
        if not environment:
            raise ValueError("Specify the Conda environment name.")
        return [
            conda_command,
            "run",
            "--no-capture-output",
            "-n",
            environment,
            "python",
            "-u",
            "-c",
            script,
        ]

    def remote_fields(self) -> dict[str, str]:
        return {
            "host": self.remote_host.GetValue().strip(),
            "user": self.remote_user.GetValue().strip(),
            "partition": self.remote_partition.GetValue().strip(),
            "walltime": self.remote_walltime.GetValue().strip(),
            "scratch_base": self.remote_scratch.GetValue().strip(),
            "conda": self.remote_conda.GetValue().strip(),
            "stage_mode": (
                remote_run.STAGE_MODE_RSYNC
                if self.remote_stage_rsync.GetValue()
                else remote_run.STAGE_MODE_MOUNT
            ),
        }

    def build_remote_command(
        self, local_command: list[str], workdir: Path, mode: str
    ) -> list[str]:
        """Generate the submission scripts and return the local driver call."""
        if self.use_current_python.GetValue():
            raise ValueError(
                "'Use the Python/environment that launched this GUI' cannot be "
                "honoured on the remote server; uncheck it so the remote "
                "conda environment is used."
            )
        conda_env = self.conda_env.GetValue().strip()
        if not conda_env:
            raise ValueError("Specify the Conda environment name.")
        fields = self.remote_fields()
        engine = Path(self.engine_file.get()).expanduser().resolve()
        driver = remote_run.write_remote_scripts(
            workdir,
            mode=mode,
            local_command=local_command,
            engine_file=engine,
            result_dir=workdir,
            host=fields["host"],
            partition=fields["partition"],
            walltime=fields["walltime"],
            scratch_base=fields["scratch_base"],
            conda=fields["conda"],
            conda_env=conda_env,
            cpus=self.mpi_processes.GetValue(),
            user=fields["user"],
            stage_mode=fields["stage_mode"],
        )
        return ["bash", str(driver)]

    def remote_environment_check_command(self) -> list[str]:
        fields = self.remote_fields()
        remote_run.validate_remote_fields(**fields)
        environment = self.conda_env.GetValue().strip()
        if not environment:
            raise ValueError("Specify the Conda environment name.")
        probe = (
            'echo "[REMOTE] Executables on $(hostname):"; '
            "for p in gmx sbatch squeue sacct; do "
            "printf '  %-10s %s\\n' \"$p\" \"$(command -v $p || echo MISSING)\"; "
            "done; "
            f"{fields['conda']} run --no-capture-output -n "
            f"{shlex.quote(environment)} python -u -c "
            f"{shlex.quote(ENV_CHECK_SCRIPT)}"
        )
        return [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=15",
            remote_run._ssh_target(fields["host"], fields["user"]),
            probe,
        ]

    def remote_summary_text(self, command: list[str], cwd: Path, workdir: Path) -> str:
        fields = self.remote_fields()
        staging = (
            "rsync over SSH (inputs pushed, results pulled)"
            if fields["stage_mode"] == remote_run.STAGE_MODE_RSYNC
            else "shared Data mount on the server"
        )
        return (
            f"Remote execution on {fields['host']} via Slurm\n\n"
            f"SSH user:   {fields['user'] or '(SSH config default)'}\n"
            f"Partition:  {fields['partition']}\n"
            f"Walltime:   {fields['walltime']}\n"
            f"CPUs:       {self.mpi_processes.GetValue()}\n"
            f"Scratch:    {fields['scratch_base']}/wxmmpbsa_XXXXXXXXXX\n"
            f"Remote env: {fields['conda']} run -n {self.conda_env.GetValue().strip()}\n"
            f"Staging:    {staging}\n\n"
            f"Inputs and the engine are staged to scratch; results are copied "
            f"back to:\n{workdir}\n\n"
            f"Local driver command:\n{format_command(command)}\n\n"
            f"Generated scripts:\n"
            f"  {workdir / remote_run.DRIVER_SCRIPT_NAME}\n"
            f"  {workdir / remote_run.JOB_SCRIPT_NAME}"
        )

    # ------------------------------ subprocess management ----------------

    def start_process(
        self,
        command: list[str],
        cwd: Path,
        task: str,
        log_path: Path | None,
        workdir: Path | None = None,
        is_remote: bool = False,
    ) -> None:
        self.process_seq += 1
        seq = self.process_seq
        self.progress.SetValue(0)
        self.status_text.SetLabel(f"Starting [#{seq}] …")
        self.SetStatusText("Starting process")
        self.open_folder_button.Enable(bool(self.workdir.get()))
        header = (
            "\n" + "=" * 88 + "\n"
            + f"[GUI] [#{seq}] {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            + f"[GUI] [#{seq}] cwd: {cwd}\n"
            + f"[GUI] [#{seq}] command: {format_command(command)}\n"
            + "=" * 88 + "\n"
        )
        self.append_log(header)

        log_handle: Any = None
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_handle = log_path.open("at", encoding="utf-8")
            log_handle.write(header)
            log_handle.flush()

        environment = os.environ.copy()
        environment["PYTHONUNBUFFERED"] = "1"
        try:
            popen = subprocess.Popen(
                command,
                cwd=str(cwd),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                start_new_session=(os.name == "posix"),
                env=environment,
            )
        except Exception:
            if log_handle is not None:
                log_handle.close()
            raise

        record: dict[str, Any] = {
            "seq": seq,
            "popen": popen,
            "thread": None,
            "task": task,
            "log_handle": log_handle,
            "stop_requested": False,
            "remote_job_id": None,
            "workdir": workdir,
            "is_remote": is_remote,
        }
        self.processes.append(record)
        record["thread"] = threading.Thread(
            target=self.read_process_output, args=(record,), daemon=True
        )
        record["thread"].start()
        self.update_controls()

    def read_process_output(self, record: dict[str, Any]) -> None:
        popen = record["popen"]
        assert popen.stdout is not None
        for line in popen.stdout:
            wx.CallAfter(self.handle_process_line, record, line)
            log_handle = record.get("log_handle")
            if log_handle is not None:
                log_handle.write(line)
                log_handle.flush()
        return_code = popen.wait()
        wx.CallAfter(self.process_finished, record, return_code)

    def handle_process_line(self, record: dict[str, Any], line: str) -> None:
        job_match = re.search(r"\[REMOTE\] Slurm JobID:\s*(\d+)", line)
        if job_match:
            record["remote_job_id"] = job_match.group(1)
        if len(self.processes) > 1 and line.strip():
            self.append_log(f"[#{record['seq']}] " + line)
        else:
            self.append_log(line)
        if not self.processes or record is not self.processes[-1]:
            return
        progress_match = re.search(r"\b(\d{1,3})%\|", line)
        if progress_match:
            self.progress.SetValue(min(100, int(progress_match.group(1))))
        lower = line.lower()
        stages = [
            ("preparation completed", "Preparation completed"),
            ("topology validation passed", "Topology validated"),
            ("calculating complex contribution", "PB/GB: complex"),
            ("calculating receptor contribution", "PB/GB: receptor"),
            ("calculating ligand contribution", "PB/GB: ligand"),
            ("parsing results", "Parsing results"),
            ("binding-energy summary", "Summary generated"),
        ]
        for needle, label in stages:
            if needle in lower:
                self.status_text.SetLabel(f"[#{record['seq']}] {label}")
                self.SetStatusText(f"[#{record['seq']}] {label}")
                if "calculating" in needle:
                    self.progress.SetValue(0)
                break

    def process_finished(self, record: dict[str, Any], return_code: int) -> None:
        seq = record["seq"]
        task = record["task"]
        stopped = record["stop_requested"]
        if record in self.processes:
            self.processes.remove(record)
        log_handle = record.get("log_handle")
        if log_handle is not None:
            log_handle.close()
            record["log_handle"] = None
        self.update_controls()
        self.progress.SetValue(100 if return_code == 0 else 0)

        if stopped:
            message = f"[#{seq}] Stopped (exit code {return_code})"
            self.status_text.SetLabel(message)
            self.SetStatusText(message)
            self.append_log(f"\n[GUI] {message}\n")
            return

        if return_code == 0:
            done = "Environment check passed" if task == "environment" else "Completed successfully"
            message = f"[#{seq}] {done}"
            self.status_text.SetLabel(message)
            self.SetStatusText(message)
            self.append_log(f"\n[GUI] {message}.\n")
            if task == "workflow":
                self.load_summary_after_run(workdir=record.get("workdir"))
        else:
            message = f"[#{seq}] Process failed with exit code {return_code}"
            self.status_text.SetLabel(message)
            self.SetStatusText(message)
            self.append_log(f"\n[GUI] {message}.\n")
            # If gmx_MMPBSA produced a complete result but an older engine's
            # summary parser returned nonzero, show the scientifically valid
            # result rather than hiding it.
            if task == "workflow":
                self.load_summary_after_run(silent=True, workdir=record.get("workdir"))
            wx.Bell()

    def load_summary_after_run(
        self, *, silent: bool = False, workdir: Path | None = None
    ) -> None:
        if workdir is None:
            workdir = self.last_workdir
        if workdir is None:
            return
        result = workdir / "FINAL_RESULTS_MMPBSA.dat"
        if not result.is_file():
            return
        try:
            summary_path = workdir / "binding_energy_summary.json"
            if summary_path.is_file():
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
            else:
                summary = parse_binding_energy(result)
                summary_path.write_text(
                    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
                )
            self.result_file.set(result)
            self.last_result = result
            self.display_summary(summary)
            self.open_result_button.Enable(True)
            self.notebook.SetSelection(2)
        except Exception as exc:
            if not silent:
                wx.MessageBox(str(exc), "Could not load result summary", wx.OK | wx.ICON_WARNING)

    def display_summary(self, summary: dict[str, Any]) -> None:
        stats = summary.get("delta_total_kcal_per_mol", {})
        if not isinstance(stats, dict):
            stats = {}

        def shown(value: Any) -> str:
            if value is None or value == "":
                return "n/a"
            if isinstance(value, float):
                return f"{value:.6f}"
            return str(value)

        values = {
            "model": summary.get("solvent_model", "unknown"),
            "frames": summary.get("frames", "n/a"),
            "average": stats.get("average"),
            "sd_prop": stats.get("standard_deviation_propagated"),
            "sd": stats.get("standard_deviation"),
            "sem_prop": stats.get("standard_error_propagated"),
            "sem": stats.get("standard_error"),
            "label": summary.get("matched_label", "n/a"),
        }
        for key, value in values.items():
            self.summary_fields[key].SetValue(shown(value))

    def append_log(self, text: str) -> None:
        self.log.AppendText(text.replace("\r", ""))
        self.log.ShowPosition(self.log.GetLastPosition())


class WorkflowApp(wx.App):
    def __init__(self, startup: argparse.Namespace) -> None:
        self.startup = startup
        super().__init__(redirect=False)

    def OnInit(self) -> bool:  # noqa: N802 - wxPython API
        frame = MainFrame(self.startup)
        frame.Show()
        self.SetTopWindow(frame)
        return True


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="wxPython GUI for Desmond gmx_MMPBSA")
    parser.add_argument("--profile", help="Load a saved JSON profile")
    parser.add_argument("--engine", help="Override the bundled workflow engine")
    parser.add_argument("--workdir", help="Initial work/output directory")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    app = WorkflowApp(args)
    app.MainLoop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
