#!/bin/bash
# Session end hook: Final cleanup and state persistence
# Receives JSON via stdin with reason

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract reason
REASON=$(echo "$INPUT" | jq -r '.reason // "unknown"' 2>/dev/null)

cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."

# Log session end in memory
TIMESTAMP=$(date +%s)
npx @claude-flow/cli@latest memory store \
    --namespace sessions \
    --key "session-end-${TIMESTAMP}" \
    --value "Session ended: $REASON" 2>/dev/null || true

# Final session end with full state persistence
npx @claude-flow/cli@latest hooks session-end --save-state --export-metrics --consolidate || {
    echo "[session-end] Warning: session-end hook failed" >&2
}
