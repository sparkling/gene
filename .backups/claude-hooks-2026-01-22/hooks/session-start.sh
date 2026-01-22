#!/bin/bash
# Session start hook: Initialize daemon and restore session
# Receives JSON via stdin with session_id

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract session_id
SESSION_ID=$(echo "$INPUT" | jq -r '.session_id // empty' 2>/dev/null)

cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."

# Start the daemon quietly
npx @claude-flow/cli@latest daemon start --quiet 2>/dev/null || {
    echo "[session-start] Note: daemon may already be running" >&2
}

# Restore session if session_id is available
if [ -n "$SESSION_ID" ]; then
    npx @claude-flow/cli@latest hooks session-restore --session-id "$SESSION_ID" || {
        echo "[session-start] Warning: session restore failed for $SESSION_ID" >&2
    }
fi

# Log session start in memory
TIMESTAMP=$(date +%s)
npx @claude-flow/cli@latest memory store \
    --namespace sessions \
    --key "session-start-${TIMESTAMP}" \
    --value "Session started: ${SESSION_ID:-unknown}" 2>/dev/null || true
