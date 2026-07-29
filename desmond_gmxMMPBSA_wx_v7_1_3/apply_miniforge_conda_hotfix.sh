#!/usr/bin/env bash
set -euo pipefail

TARGET_ROOT=${1:-"$HOME/desmond_gmxMMPBSA"}
TARGET_ROOT=$(readlink -f -- "$TARGET_ROOT")
GUI="$TARGET_ROOT/desmond_gmxmmpbsa_wx.py"
STAMP=$(date +%Y%m%d_%H%M%S)

if [[ ! -f "$GUI" ]]; then
    printf 'ERROR: GUI source not found: %s\n' "$GUI" >&2
    exit 1
fi

cp -a -- "$GUI" "${GUI}.backup_${STAMP}"
printf 'Backed up %s -> %s\n' "$GUI" "${GUI}.backup_${STAMP}"

python3 - "$GUI" <<'PY'
from pathlib import Path
import sys

p = Path(sys.argv[1])
text = p.read_text(encoding="utf-8")

if 'GUI_VERSION = "1.0.3"' in text and 'def is_usable_conda(' in text:
    print("The Miniforge/Conda hotfix is already present.")
    raise SystemExit(0)

text = text.replace('GUI_VERSION = "1.0.2"', 'GUI_VERSION = "1.0.3"', 1)
text = text.replace('GUI_VERSION = "1.0.1"', 'GUI_VERSION = "1.0.3"', 1)

old = '''def discover_conda() -> str:
    candidates: list[str | None] = [
        os.environ.get("CONDA_EXE"),
        shutil.which("conda"),
        str(Path.home() / "miniforge3" / "bin" / "conda"),
        str(Path.home() / "miniconda3" / "bin" / "conda"),
        str(Path.home() / "anaconda3" / "bin" / "conda"),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    return "conda"
'''
new = '''def normalize_executable(candidate: str | Path | None) -> str | None:
    """Resolve an executable name/path without trusting shell aliases or functions."""

    if candidate is None:
        return None
    text = os.path.expandvars(os.path.expanduser(str(candidate).strip()))
    if not text:
        return None
    if os.sep not in text:
        found = shutil.which(text)
        if not found:
            return None
        text = found
    path = Path(text)
    if not path.is_file() or not os.access(path, os.X_OK):
        return None
    return str(path.resolve())


def is_usable_conda(candidate: str | Path | None) -> bool:
    """Return True only for a real Conda executable supporting ``conda run``."""

    executable = normalize_executable(candidate)
    if not executable:
        return False
    try:
        version = subprocess.run(
            [executable, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=8,
            check=False,
        )
        if version.returncode != 0 or not re.search(
            r"\\bconda\\s+\\d", version.stdout or "", flags=re.I
        ):
            return False
        run_help = subprocess.run(
            [executable, "run", "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=8,
            check=False,
        )
        return run_help.returncode == 0 and "conda run" in (run_help.stdout or "").lower()
    except (OSError, subprocess.SubprocessError):
        return False


def conda_candidates() -> list[str]:
    """Return likely Conda executables, preferring canonical local installs."""

    home = Path.home()
    candidates: list[str | None] = [
        os.environ.get("CONDA_EXE"),
        str(home / "miniforge3" / "bin" / "conda"),
        str(home / "mambaforge" / "bin" / "conda"),
        str(home / "miniconda3" / "bin" / "conda"),
        str(home / "anaconda3" / "bin" / "conda"),
        shutil.which("conda"),
    ]
    unique: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        executable = normalize_executable(candidate)
        if executable and executable not in seen:
            seen.add(executable)
            unique.append(executable)
    return unique


def discover_conda() -> str:
    """Find a genuine Conda/Miniforge executable, not an activation wrapper."""

    for candidate in conda_candidates():
        if is_usable_conda(candidate):
            return candidate
    return ""


def resolve_conda_executable(preferred: str | Path | None = None) -> str:
    """Resolve a selected Conda path, falling back to a valid local install."""

    preferred_executable = normalize_executable(preferred)
    if preferred_executable and is_usable_conda(preferred_executable):
        return preferred_executable
    discovered = discover_conda()
    if discovered:
        return discovered
    shown = str(preferred).strip() if preferred is not None else ""
    detail = f" The selected path was: {shown}" if shown else ""
    raise FileNotFoundError(
        "No usable Conda/Miniforge executable was found." + detail + "\\n\\n"
        "Select the real executable, normally:\\n"
        "    ~/miniforge3/bin/conda\\n\\n"
        "Do not select a shell wrapper that only calls 'conda activate'."
    )
'''
if old not in text:
    raise SystemExit("ERROR: expected GUI 1.0.2 Conda-discovery block was not found")
text = text.replace(old, new, 1)

text = text.replace(
    '("Conda executable", self.conda_file),',
    '("Conda/Miniforge executable", self.conda_file),',
    1,
)
text = text.replace(
    "then execute the scientific engine through 'conda run -n ",
    "then execute the scientific engine through Miniforge/Conda 'conda run -n ",
    1,
)

old_apply = '        self.conda_file.set(data.get("conda", discover_conda()))\n'
new_apply = '''        saved_conda = str(data.get("conda", "")).strip()
        try:
            resolved_conda = resolve_conda_executable(saved_conda)
        except FileNotFoundError:
            resolved_conda = saved_conda or discover_conda()
        self.conda_file.set(resolved_conda)
'''
if old_apply not in text:
    raise SystemExit("ERROR: expected saved-Conda settings line was not found")
text = text.replace(old_apply, new_apply, 1)

old_prefix = '''    def python_prefix(self) -> list[str]:
        engine = Path(self.engine_file.get()).expanduser().resolve()
        if not engine.is_file():
            raise FileNotFoundError(f"Workflow engine not found: {engine}")
        if self.use_current_python.GetValue():
            return [sys.executable, "-u", str(engine)]
        conda_text = self.conda_file.get() or discover_conda()
        conda = Path(conda_text).expanduser()
        if conda_text != "conda" and not conda.is_file():
            raise FileNotFoundError(f"Conda executable not found: {conda}")
        conda_command = str(conda if conda_text != "conda" else "conda")
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
'''
new_prefix = '''    def resolved_conda_command(self) -> str:
        preferred = self.conda_file.get()
        resolved = resolve_conda_executable(preferred)
        preferred_normalized = normalize_executable(preferred)
        if preferred_normalized != resolved:
            self.conda_file.set(resolved)
            old = preferred or "<empty>"
            self.append_log(
                f"[GUI] Ignored unusable Conda launcher {old!r}; "
                f"using {resolved!r}.\\n"
            )
        return resolved

    def python_prefix(self) -> list[str]:
        engine = Path(self.engine_file.get()).expanduser().resolve()
        if not engine.is_file():
            raise FileNotFoundError(f"Workflow engine not found: {engine}")
        if self.use_current_python.GetValue():
            return [sys.executable, "-u", str(engine)]
        conda_command = self.resolved_conda_command()
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
'''
if old_prefix not in text:
    raise SystemExit("ERROR: expected python_prefix block was not found")
text = text.replace(old_prefix, new_prefix, 1)

old_env = '''        conda_text = self.conda_file.get() or discover_conda()
        environment = self.conda_env.GetValue().strip()
        if not environment:
            raise ValueError("Specify the Conda environment name.")
        return [
            conda_text,
            "run",
            "--no-capture-output",
            "-n",
            environment,
            "python",
            "-u",
            "-c",
            script,
        ]
'''
new_env = '''        conda_command = self.resolved_conda_command()
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
'''
if old_env not in text:
    raise SystemExit("ERROR: expected environment-check Conda block was not found")
text = text.replace(old_env, new_env, 1)

p.write_text(text, encoding="utf-8")
print(f"Patched {p}")
PY

# Repair the persisted GUI setting immediately when standard Miniforge exists.
SETTINGS="$HOME/.config/desmond_gmxmmpbsa_wx/settings.json"
REAL_CONDA="$HOME/miniforge3/bin/conda"
if [[ -x "$REAL_CONDA" && -f "$SETTINGS" ]]; then
    cp -a -- "$SETTINGS" "${SETTINGS}.backup_${STAMP}"
    python3 - "$SETTINGS" "$REAL_CONDA" <<'PY'
import json
import sys
from pathlib import Path
p = Path(sys.argv[1])
conda = str(Path(sys.argv[2]).resolve())
data = json.loads(p.read_text(encoding="utf-8"))
data["conda"] = conda
p.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
print(f"Updated saved Conda path in {p}: {conda}")
PY
fi

python3 -m py_compile "$GUI"

if [[ -x "$HOME/miniforge3/bin/conda" ]]; then
    "$HOME/miniforge3/bin/conda" --version
    "$HOME/miniforge3/bin/conda" run --help >/dev/null
    printf 'Validated Miniforge executable: %s\n' "$HOME/miniforge3/bin/conda"
fi

printf '\nHotfix installed. Launch with:\n  %s/run_desmond_gmxmmpbsa_wx.sh\n' "$TARGET_ROOT"
