#!/bin/bash
# Pre-command hook: Assess risk before executing bash commands
# Receives JSON via stdin with tool_input.command

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract command from tool_input
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)

if [ -z "$COMMAND" ]; then
    exit 0
fi

# Skip our own hooks commands to avoid recursion
case "$COMMAND" in
    *"@claude-flow/cli"*|*"hooks"*|*"memory store"*)
        exit 0
        ;;
esac

# Log for debugging
# echo "[pre-command] Assessing: ${COMMAND:0:80}..." >&2

# Call claude-flow hooks pre-command
cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."
npx @claude-flow/cli@latest hooks pre-command --command "$COMMAND" || {
    echo "[pre-command] Warning: pre-command hook failed" >&2
}
