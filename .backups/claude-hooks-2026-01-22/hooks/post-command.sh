#!/bin/bash
# Post-command hook: Record command execution outcome for learning
# Receives JSON via stdin with tool_input.command and tool_response

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract command from tool_input
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)

# Check if tool_response exists (indicates success)
HAS_RESPONSE=$(echo "$INPUT" | jq -r 'if .tool_response then "true" else "false" end' 2>/dev/null)

if [ -z "$COMMAND" ]; then
    exit 0
fi

# Skip our own hooks commands to avoid recursion
case "$COMMAND" in
    *"@claude-flow/cli"*|*"hooks"*|*"memory store"*)
        exit 0
        ;;
esac

SUCCESS="true"
if [ "$HAS_RESPONSE" = "false" ]; then
    SUCCESS="false"
fi

# Log for debugging
# echo "[post-command] Recording: ${COMMAND:0:60}... (success=$SUCCESS)" >&2

# Call claude-flow hooks post-command
cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."
npx @claude-flow/cli@latest hooks post-command --command "$COMMAND" --success "$SUCCESS" || {
    echo "[post-command] Warning: post-command hook failed" >&2
}
