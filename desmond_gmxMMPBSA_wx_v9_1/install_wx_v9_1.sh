#!/usr/bin/env bash
set -euo pipefail

SOURCE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
TARGET_ROOT=${1:-"$HOME/desmond_gmxMMPBSA"}
TARGET_ROOT=$(readlink -f -- "$TARGET_ROOT")
TARGET_APP="$TARGET_ROOT/app"
STAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "$TARGET_APP"

backup_if_present() {
    local path=$1
    if [[ -e "$path" ]]; then
        cp -a -- "$path" "${path}.backup_${STAMP}"
        printf 'Backed up %s -> %s\n' "$path" "${path}.backup_${STAMP}"
    fi
}

backup_if_present "$TARGET_APP/desmond_gmxmmpbsa.py"
backup_if_present "$TARGET_APP/repair_synthetic_caps.py"
backup_if_present "$TARGET_APP/steric_preflight.py"
backup_if_present "$TARGET_APP/protein.ff19SB.xml"
backup_if_present "$TARGET_ROOT/desmond_gmxmmpbsa_wx.py"
backup_if_present "$TARGET_ROOT/run_desmond_gmxmmpbsa_wx.sh"
backup_if_present "$TARGET_ROOT/run_desmond_gmxmmpbsa.sh"

install -m 0755 "$SOURCE_ROOT/app/desmond_gmxmmpbsa.py" "$TARGET_APP/desmond_gmxmmpbsa.py"
install -m 0755 "$SOURCE_ROOT/app/repair_synthetic_caps.py" "$TARGET_APP/repair_synthetic_caps.py"
install -m 0755 "$SOURCE_ROOT/app/steric_preflight.py" "$TARGET_APP/steric_preflight.py"
install -m 0644 "$SOURCE_ROOT/app/protein.ff19SB.xml" "$TARGET_APP/protein.ff19SB.xml"
install -m 0755 "$SOURCE_ROOT/desmond_gmxmmpbsa_wx.py" "$TARGET_ROOT/desmond_gmxmmpbsa_wx.py"
install -m 0755 "$SOURCE_ROOT/run_desmond_gmxmmpbsa_wx.sh" "$TARGET_ROOT/run_desmond_gmxmmpbsa_wx.sh"
install -m 0755 "$SOURCE_ROOT/run_desmond_gmxmmpbsa.sh" "$TARGET_ROOT/run_desmond_gmxmmpbsa.sh"

printf '\nInstalled engine 9.1.0 and wxPython GUI 1.0.6 into:\n  %s\n\n' "$TARGET_ROOT"
printf 'Launch with:\n  %s/run_desmond_gmxmmpbsa_wx.sh\n' "$TARGET_ROOT"
