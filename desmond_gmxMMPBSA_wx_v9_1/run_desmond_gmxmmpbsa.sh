#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHON_BIN=${PYTHON_BIN:-python}
exec "$PYTHON_BIN" -u "$ROOT/app/desmond_gmxmmpbsa.py" "$@"
