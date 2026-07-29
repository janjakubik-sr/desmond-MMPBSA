#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
APP_DIR="$HOME/.local/share/applications"
DESKTOP_FILE="$APP_DIR/desmond-gmxmmpbsa.desktop"
mkdir -p "$APP_DIR"
cat > "$DESKTOP_FILE" <<EOF
[Desktop Entry]
Type=Application
Name=Desmond gmx_MMPBSA
Comment=Prepare Desmond GRO/XTC trajectories and run gmx_MMPBSA
Exec=$ROOT/run_desmond_gmxmmpbsa_wx.sh
Terminal=false
Categories=Education;Science;Chemistry;
StartupNotify=true
EOF
chmod 0644 "$DESKTOP_FILE"
printf 'Installed desktop launcher: %s\n' "$DESKTOP_FILE"
