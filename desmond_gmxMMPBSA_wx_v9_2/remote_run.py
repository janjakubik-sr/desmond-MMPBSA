#!/usr/bin/env python3
"""Slurm-side script generator for running the workflow on a remote server.

The GUI executes the local driver script (for example
``fgusrvvs03_remote.sh``) locally.  That driver allocates a scratch
directory on the server, uploads the workflow engine and the job script
over SSH, and finally executes the matching server-side job script (for
example ``fgusrvvs03_job.sh``) on the server.  Both script names are
derived from the SSH target host so the files always identify the server
they belong to.  The job script submits the engine to Slurm, streams
the job log back through the ssh session, and exits with the Slurm job's
exit status -- so the GUI sees one ordinary local process whose stdout
carries the engine's live output.

Two staging modes are supported:

* ``mount`` (alkyone): the server sees the same sshfs-mounted Data
  directory as the local machine, so the job script copies the inputs
  off the mount into scratch and copies the results back onto the
  mounted work directory itself.
* ``rsync`` (e.g. fgusrvvs03): the server has no shared mount, so the
  local driver pushes the input files (or the whole prepared work
  directory) to scratch with rsync before the job and pulls the results
  back with rsync afterwards; scratch is removed on success.
"""

from __future__ import annotations

import re
import shlex
import stat
from pathlib import Path
from typing import Any, Sequence

REMOTE_CONDA_DEFAULT = "~/miniconda3/bin/conda"
HOST_DEFAULT = "alkyone"
PARTITION_DEFAULT = "LocalQ"
WALLTIME_DEFAULT = "48:00:00"
SCRATCH_BASE_DEFAULT = "/scratch/roshi"

STAGE_MODE_MOUNT = "mount"
STAGE_MODE_RSYNC = "rsync"


def script_names(host: str) -> tuple[str, str]:
    """Return (driver script, job script) file names for one SSH host.

    Names are derived from the target host so that a work directory never
    holds a script named after the wrong server (the historical fixed
    ``alkyone_job.sh`` name was carried over to every other host).  The
    default host ``alkyone`` keeps its original file names; every other
    host gets scripts named after itself, e.g. ``fgusrvvs03_remote.sh``
    and ``fgusrvvs03_job.sh``.
    """

    prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", host.strip().lower()).strip("._-")
    if not prefix:
        prefix = "remote"
    return f"{prefix}_remote.sh", f"{prefix}_job.sh"

_WALLTIME_PATTERN = re.compile(r"^(?:\d+-)?\d{1,2}:\d{2}:\d{2}$")
_SAFE_FIELD_PATTERN = re.compile(r"^[A-Za-z0-9._~/-]+$")


def _q(value: Any) -> str:
    """Shell-quote one token."""
    return shlex.quote(str(value))


_PATH_OPTIONS = {"--gro", "--xtc", "--mol2", "--ffxml"}
_WORKDIR_OPTIONS = {"--outdir", "--workdir"}


def remap_engine_command(
    local_command: Sequence[str],
    *,
    remote_conda: str,
    remote_env: str,
) -> tuple[list[str], dict[str, str]]:
    """Rewrite a local engine command for execution inside remote scratch.

    The Slurm job runs with ``--chdir`` set to the remote scratch
    directory, so every staged artefact is addressed by a relative path:
    inputs under ``in/``, results under ``out/`` (mode ``all``) or the
    prepared directory under ``wd/`` (mode ``run``), and the uploaded
    engine under ``engine/``.

    Returns the remote command tokens plus a mapping of remote staging
    names (relative to ``in/``) to the original local mount paths that
    must be copied by the job script.
    """

    tokens = [str(t) for t in local_command]
    try:
        u_index = tokens.index("-u")
    except ValueError as exc:
        raise ValueError("Could not locate 'python -u' in the engine command.") from exc
    if u_index + 1 >= len(tokens):
        raise ValueError("Engine command does not name the workflow script.")
    engine_path = Path(tokens[u_index + 1])

    remote = [
        remote_conda,
        "run",
        "--no-capture-output",
        "-n",
        remote_env,
        "python",
        "-u",
        f"engine/{engine_path.name}",
    ]

    input_map: dict[str, str] = {}
    rest = tokens[u_index + 2 :]
    index = 0
    while index < len(rest):
        token = rest[index]
        if token in _PATH_OPTIONS and index + 1 < len(rest):
            original = rest[index + 1]
            name = Path(original).name
            if name in input_map and input_map[name] != original:
                raise ValueError(
                    f"Two different inputs share the file name {name!r}; "
                    "rename one of them before remote submission."
                )
            input_map[name] = original
            remote += [token, f"in/{name}"]
            index += 2
            continue
        if token in _WORKDIR_OPTIONS and index + 1 < len(rest):
            target = "wd" if token == "--workdir" else "out"
            remote += [token, target]
            index += 2
            continue
        remote.append(token)
        index += 1
    return remote, input_map


def validate_remote_fields(
    *,
    host: str,
    partition: str,
    walltime: str,
    scratch_base: str,
    conda: str,
    user: str = "",
    stage_mode: str = STAGE_MODE_MOUNT,
) -> None:
    """Reject malformed remote-submission fields before touching ssh."""
    if not host or not _SAFE_FIELD_PATTERN.match(host):
        raise ValueError(f"Invalid SSH host name: {host!r}")
    if user and not _SAFE_FIELD_PATTERN.match(user):
        raise ValueError(f"Invalid remote user name: {user!r}")
    if stage_mode not in {STAGE_MODE_MOUNT, STAGE_MODE_RSYNC}:
        raise ValueError(
            f"Unknown staging mode: {stage_mode!r} "
            f"(expected '{STAGE_MODE_MOUNT}' or '{STAGE_MODE_RSYNC}')"
        )
    if not partition or not _SAFE_FIELD_PATTERN.match(partition):
        raise ValueError(f"Invalid Slurm partition: {partition!r}")
    if not _WALLTIME_PATTERN.match(walltime):
        raise ValueError(
            "Invalid Slurm walltime. Use HH:MM:SS or D-HH:MM:SS, e.g. 48:00:00."
        )
    if not scratch_base.startswith("/"):
        raise ValueError(
            f"Remote scratch base must be an absolute path: {scratch_base!r}"
        )
    if not scratch_base or "//" in scratch_base:
        raise ValueError(f"Invalid remote scratch base: {scratch_base!r}")
    if not conda.strip():
        raise ValueError("Specify the conda executable on the remote server.")


def _dq(value: Any) -> str:
    """Escape one token for embedding inside a double-quoted shell string."""
    return (
        str(value)
        .replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("$", "\\$")
        .replace("`", "\\`")
    )


def _ssh_target(host: str, user: str) -> str:
    """Return the ssh destination, prefixing an explicit user when given."""
    return f"{user}@{host}" if user else host


def _push_lines(
    *, mode: str, input_map: dict[str, str], result_dir: Path | None
) -> str:
    """Driver-side rsync commands that upload inputs to remote scratch."""
    if mode == "run":
        if result_dir is None:
            raise ValueError("rsync staging of mode 'run' needs the local workdir.")
        return (
            'echo "[REMOTE] Uploading the prepared directory ..."\n'
            "ssh -o BatchMode=yes \"$SSH_TARGET\" \"mkdir -p '$REMOTE_TMP/wd'\"\n"
            f"rsync -a {_q(result_dir)}/ \"$SSH_TARGET:$REMOTE_TMP/wd/\""
        )
    lines = [
        'echo "[REMOTE] Uploading input files ..."',
        "ssh -o BatchMode=yes \"$SSH_TARGET\" \"mkdir -p '$REMOTE_TMP/in'\"",
    ]
    for staged_name, original in sorted(input_map.items()):
        lines.append(
            f"rsync -a -- {_q(original)} "
            f"\"$SSH_TARGET:$REMOTE_TMP/in/{_dq(staged_name)}\""
        )
    return "\n".join(lines)


def driver_script(
    *,
    host: str,
    scratch_base: str,
    engine_dir: Path,
    workdir: Path,
    user: str = "",
    stage_mode: str = STAGE_MODE_MOUNT,
    mode: str = "all",
    input_map: dict[str, str] | None = None,
    result_dir: Path | None = None,
) -> str:
    """Local wrapper: allocate scratch, upload engine + job script, hand over."""
    target = _ssh_target(host, user)
    driver_name, job_name = script_names(host)
    rsync_stage = stage_mode == STAGE_MODE_RSYNC

    push_block = _push_lines(
        mode=mode, input_map=input_map or {}, result_dir=result_dir
    ) if rsync_stage else ":"

    if rsync_stage:
        if mode == "run":
            pull_src = '"$SSH_TARGET:$REMOTE_TMP/wd/"'
        else:
            pull_src = '"$SSH_TARGET:$REMOTE_TMP/out/"'
        handover = f"""set +e
ssh -o BatchMode=yes "$SSH_TARGET" "bash '${{REMOTE_TMP}}/{job_name}' '${{REMOTE_TMP}}'"
STATUS=$?
set -e
echo "[REMOTE] Remote job finished with exit code: $STATUS"
echo "[REMOTE] Pulling results back to: $RESULT_DIR"
rsync -a {pull_src} "$RESULT_DIR/" || true
for extra in slurm.out job.log wxmmpbsa.sbatch job.exit; do
    rsync -a "$SSH_TARGET:$REMOTE_TMP/$extra" "$RESULT_DIR/" 2>/dev/null || true
done
if [ "$STATUS" -eq 0 ]; then
    echo "[REMOTE] Success; removing scratch $REMOTE_TMP"
    ssh -o BatchMode=yes "$SSH_TARGET" "rm -rf '$REMOTE_TMP'" || true
    exit 0
fi
echo "[REMOTE] Failure; scratch kept for inspection: $REMOTE_TMP"
exit "$STATUS"
"""
    else:
        handover = (
            f'exec ssh -o BatchMode=yes "$SSH_TARGET" '
            f'"bash \'${{REMOTE_TMP}}/{job_name}\' \'${{REMOTE_TMP}}\'"'
        )

    return f"""#!/usr/bin/env bash
# {driver_name} -- local driver generated by the wx GUI.
set -euo pipefail
HOST={_q(host)}
SSH_TARGET={_q(target)}
SCRATCH_BASE={_q(scratch_base)}
ENGINE_DIR={_q(engine_dir)}
JOB_SRC={_q(Path(workdir) / job_name)}
RESULT_DIR={_q(result_dir if rsync_stage else workdir)}

echo "[REMOTE] Allocating scratch directory on $SSH_TARGET ..."
REMOTE_TMP=$(ssh -o BatchMode=yes -o ConnectTimeout=15 "$SSH_TARGET" \\
    "mktemp -d '${{SCRATCH_BASE}}/wxmmpbsa_XXXXXXXXXX'")
echo "[REMOTE] Remote scratch directory: ${{REMOTE_TMP}}"
echo "[REMOTE] Uploading workflow engine and job script ..."
rsync -a --exclude='__pycache__' --exclude='*.backup_*' \\
    "$ENGINE_DIR/" "${{SSH_TARGET}}:${{REMOTE_TMP}}/engine/"
rsync -a "$JOB_SRC" "${{SSH_TARGET}}:${{REMOTE_TMP}}/{job_name}"
{push_block}
echo "[REMOTE] Submitting Slurm job on $SSH_TARGET ..."
{handover}
"""


def _engine_line(tokens: Sequence[str]) -> str:
    """Shell-join remote command tokens, expanding a leading '~/' via $HOME."""

    parts: list[str] = []
    for raw in map(str, tokens):
        if raw.startswith("~/"):
            rest = (
                raw[2:]
                .replace("\\", "\\\\")
                .replace('"', '\\"')
                .replace("$", "\\$")
                .replace("`", "\\`")
            )
            parts.append(f'"${{HOME}}/{rest}"')
        else:
            parts.append(shlex.quote(raw))
    return " ".join(parts)


def _staging_lines(input_map: dict[str, str]) -> str:
    lines = [
        'echo "[REMOTE] Staging inputs from the mounted directory ..."',
        'staged_ok=1',
    ]
    for staged_name, original in sorted(input_map.items()):
        lines.append(
            f"if cp -L -- {_q(original)} \"$IN/{_q(staged_name)}\"; then :; else\n"
            f"    echo \"[REMOTE] FATAL: input not available on the server: {_q(original)}\"\n"
            "    staged_ok=0\n"
            "fi"
        )
    lines.append('[ "$staged_ok" = "1" ] || exit 3')
    return "\n".join(lines)


def job_script(
    *,
    mode: str,
    remote_command: Sequence[str],
    input_map: dict[str, str],
    result_dir: Path,
    partition: str,
    walltime: str,
    cpus: int,
    job_name: str = "wxmmpbsa",
    stage_mode: str = STAGE_MODE_MOUNT,
    script_name: str = "wxmmpbsa_job.sh",
) -> str:
    """Server-side script: stage, submit to Slurm, stream, copy back."""
    engine_line = _engine_line(remote_command)
    rsync_stage = stage_mode == STAGE_MODE_RSYNC
    if rsync_stage:
        staging = (
            'echo "[REMOTE] Inputs were uploaded to scratch by the local driver."'
        )
        copyback = 'echo "[REMOTE] The local driver will pull the results back."'
        result_note = "Results will be pulled back to: $RESULT_DIR"
    elif mode == "run":
        staging = (
            'echo "[REMOTE] Staging the prepared directory ..."\n'
            f'rsync -a {_q(result_dir)}/ "$REMOTE_TMP/wd/" \\'
            "\n    || { echo '[REMOTE] FATAL: could not stage the prepared directory'; exit 3; }"
        )
        copyback = (
            'echo "[REMOTE] Copying results back to the mounted directory ..."\n'
            'rsync -a "$REMOTE_TMP/wd/" "$RESULT_DIR/"\n'
            "for extra in slurm.out job.log wxmmpbsa.sbatch job.exit; do\n"
            '    [ -f "$REMOTE_TMP/$extra" ] && cp -f "$REMOTE_TMP/$extra" "$RESULT_DIR/"\n'
            "done"
        )
        result_note = "Results will be copied to: $RESULT_DIR"
    else:
        staging = _staging_lines(input_map)
        copyback = (
            'echo "[REMOTE] Copying results back to the mounted directory ..."\n'
            'rsync -a "$REMOTE_TMP/out/" "$RESULT_DIR/"\n'
            "for extra in slurm.out job.log wxmmpbsa.sbatch job.exit; do\n"
            '    [ -f "$REMOTE_TMP/$extra" ] && cp -f "$REMOTE_TMP/$extra" "$RESULT_DIR/"\n'
            "done"
        )
        result_note = "Results will be copied to: $RESULT_DIR"

    if rsync_stage:
        cleanup = 'exit "$CODE"'
    else:
        cleanup = (
            'if [ "$CODE" = "0" ]; then\n'
            '    echo "[REMOTE] Success; removing scratch $REMOTE_TMP"\n'
            '    rm -rf "$REMOTE_TMP"\n'
            "else\n"
            '    echo "[REMOTE] Failure; scratch kept for inspection: $REMOTE_TMP"\n'
            "fi\n"
            'exit "$CODE"'
        )

    return f"""#!/usr/bin/env bash
# {script_name} -- runs on the server login node.
set -u
REMOTE_TMP="${{1:-}}"
if [ -z "$REMOTE_TMP" ]; then
    echo "[REMOTE] FATAL: scratch directory argument missing" >&2
    exit 2
fi
IN="$REMOTE_TMP/in"
LOG="$REMOTE_TMP/job.log"
RESULT_DIR={_q(result_dir)}
mkdir -p "$IN" "$REMOTE_TMP/out" "$REMOTE_TMP/wd"
echo "[REMOTE] Scratch: $REMOTE_TMP"
echo "[REMOTE] {result_note}"

{staging}

cat > "$REMOTE_TMP/wxmmpbsa.sbatch" <<SBATCH_EOF
#!/bin/bash
#SBATCH --job-name={_q(job_name)}
#SBATCH --partition={_q(partition)}
#SBATCH --time={_q(walltime)}
#SBATCH --nodes=1
#SBATCH --ntasks={int(cpus)}
#SBATCH --cpus-per-task=1
#SBATCH --chdir=${{REMOTE_TMP}}
#SBATCH --output=${{REMOTE_TMP}}/slurm.out

set +e
{engine_line} > "${{REMOTE_TMP}}/job.log" 2>&1
STATUS=\\$?
echo \\$STATUS > "${{REMOTE_TMP}}/job.exit"
exit \\$STATUS
SBATCH_EOF

JOBID=$(sbatch --parsable "$REMOTE_TMP/wxmmpbsa.sbatch") \\
    || {{ echo "[REMOTE] FATAL: sbatch submission failed"; exit 4; }}
echo "[REMOTE] Slurm JobID: $JOBID"
echo "[REMOTE] Streaming job output ..."
touch "$LOG"
tail -n +1 -F "$LOG" &
TAILPID=$!
while [ -n "$(squeue -h -j "$JOBID" 2>/dev/null)" ]; do sleep 10; done
sleep 4
kill "$TAILPID" 2>/dev/null
wait "$TAILPID" 2>/dev/null

if [ -f "$REMOTE_TMP/job.exit" ]; then
    CODE=$(cat "$REMOTE_TMP/job.exit")
else
    CODE=$(sacct -nP -j "$JOBID" --format=ExitCode 2>/dev/null \\
        | head -1 | cut -d: -f1)
fi
CODE=${{CODE:-1}}
echo "[REMOTE] Slurm job finished with exit code: $CODE"

{copyback}

{cleanup}
"""


def write_remote_scripts(
    workdir: Path,
    *,
    mode: str,
    local_command: Sequence[str],
    engine_file: Path,
    result_dir: Path,
    host: str,
    partition: str,
    walltime: str,
    scratch_base: str,
    conda: str,
    conda_env: str,
    cpus: int,
    user: str = "",
    stage_mode: str = STAGE_MODE_MOUNT,
) -> Path:
    """Generate both scripts inside *workdir*; return the driver path."""
    validate_remote_fields(
        host=host, partition=partition, walltime=walltime,
        scratch_base=scratch_base, conda=conda,
        user=user, stage_mode=stage_mode,
    )
    remote_command, input_map = remap_engine_command(
        local_command, remote_conda=conda, remote_env=conda_env
    )
    driver_name, job_name = script_names(host)
    workdir.mkdir(parents=True, exist_ok=True)
    job_path = workdir / job_name
    job_path.write_text(
        job_script(
            mode=mode, remote_command=remote_command, input_map=input_map,
            result_dir=result_dir, partition=partition, walltime=walltime,
            cpus=cpus, stage_mode=stage_mode, script_name=job_name,
        ),
        encoding="utf-8",
    )
    driver_path = workdir / driver_name
    driver_path.write_text(
        driver_script(
            host=host, scratch_base=scratch_base,
            engine_dir=Path(engine_file).resolve().parent, workdir=workdir,
            user=user, stage_mode=stage_mode, mode=mode,
            input_map=input_map, result_dir=result_dir,
        ),
        encoding="utf-8",
    )
    for path in (driver_path, job_path):
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return driver_path
