#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

candidates=()
if [[ -n "${WX_PYTHON:-}" ]]; then
    candidates+=("$WX_PYTHON")
fi
candidates+=("/usr/bin/python3")
if command -v python3 >/dev/null 2>&1; then
    candidates+=("$(command -v python3)")
fi

for python_bin in "${candidates[@]}"; do
    [[ -x "$python_bin" ]] || continue
    if "$python_bin" -c 'import wx' >/dev/null 2>&1; then
        exec "$python_bin" "$ROOT/desmond_gmxmmpbsa_wx.py" "$@"
    fi
done

cat >&2 <<'EOF'
No Python interpreter with wxPython was found.

On Debian/Ubuntu install the distribution package:

    sudo apt install python3-wxgtk4.0

Then rerun this launcher.  To use a specific interpreter:

    WX_PYTHON=/path/to/python ./run_desmond_gmxmmpbsa_wx.sh
EOF
exit 1
