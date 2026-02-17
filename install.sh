#!/bin/bash
# =============================================================================
# STARsolo CLI — installer
# =============================================================================
# Adds the starsolo command to your PATH by symlinking into a bin directory.
#
# Usage:
#   ./install.sh                    # installs to ~/.local/bin  (default)
#   ./install.sh /usr/local/bin     # installs to /usr/local/bin (needs sudo)
#   PREFIX=~/tools ./install.sh     # installs to ~/tools/bin
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TARGET_DIR="${1:-${PREFIX:+${PREFIX}/bin}}"
TARGET_DIR="${TARGET_DIR:-$HOME/.local/bin}"

mkdir -p "$TARGET_DIR"

STARSOLO_BIN="$SCRIPT_DIR/bin/starsolo"

# Make everything executable
chmod +x "$STARSOLO_BIN"
chmod +x "$SCRIPT_DIR"/lib/*.sh 2>/dev/null || true
chmod +x "$SCRIPT_DIR"/scripts/*.sh 2>/dev/null || true

# Create symlink
ln -sf "$STARSOLO_BIN" "$TARGET_DIR/starsolo"

echo "Installed starsolo → $TARGET_DIR/starsolo"

# Check if target is in PATH
if ! echo "$PATH" | tr ':' '\n' | grep -qx "$TARGET_DIR"; then
    echo ""
    echo "NOTE: $TARGET_DIR is not in your PATH."
    echo "Add it by appending this line to your ~/.bashrc or ~/.bash_profile:"
    echo ""
    echo "  export PATH=\"$TARGET_DIR:\$PATH\""
    echo ""
fi

echo ""
echo "Verify with:  starsolo --version"
