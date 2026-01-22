#!/bin/bash
# Stop hook: Persist state when Claude stops responding
# Receives JSON via stdin with stop_hook_active

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Check if already in a stop hook to avoid recursion
STOP_ACTIVE=$(echo "$INPUT" | jq -r '.stop_hook_active // false' 2>/dev/null)

if [ "$STOP_ACTIVE" = "true" ]; then
    # Already in a stop hook, skip to avoid infinite loops
    exit 0
fi

cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."

# End session and persist state
npx @claude-flow/cli@latest hooks session-end --save-state --export-metrics --consolidate || {
    echo "[stop] Warning: session-end hook failed" >&2
}
