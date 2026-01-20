#!/bin/bash
# Statusline Dispatcher - loads the current template
# Reads template name from statusline-state and executes it

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATE_FILE="$SCRIPT_DIR/statusline-state"
TEMPLATE_DIR="$SCRIPT_DIR/statuslines"

# Read current template (default: adaptive)
CURRENT=$(cat "$STATE_FILE" 2>/dev/null | tr -d '[:space:]')
[ -z "$CURRENT" ] && CURRENT="adaptive"

# Execute the template
TEMPLATE="$TEMPLATE_DIR/${CURRENT}.sh"
if [ -f "$TEMPLATE" ]; then
    exec bash "$TEMPLATE"
else
    # Fallback: try adaptive, then minimal, then basic output
    if [ -f "$TEMPLATE_DIR/adaptive.sh" ]; then
        exec bash "$TEMPLATE_DIR/adaptive.sh"
    elif [ -f "$TEMPLATE_DIR/minimal.sh" ]; then
        exec bash "$TEMPLATE_DIR/minimal.sh"
    else
        echo "Claude Flow V3"
    fi
fi
