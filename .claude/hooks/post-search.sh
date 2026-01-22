#!/bin/bash
# Post-search hook: Record search operations (Read, Glob, Grep)
# Receives JSON via stdin with tool_name and tool_input

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract tool name and relevant input
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // empty' 2>/dev/null)

# Check if tool_response exists (indicates success)
HAS_RESPONSE=$(echo "$INPUT" | jq -r 'if .tool_response then "true" else "false" end' 2>/dev/null)

SUCCESS="true"
if [ "$HAS_RESPONSE" = "false" ]; then
    SUCCESS="false"
fi

# Log for debugging
# echo "[post-search] Recording: $TOOL_NAME (success=$SUCCESS)" >&2

# Call claude-flow hooks post-command for search operations
cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."
npx @claude-flow/cli@latest hooks post-command --command "search:$TOOL_NAME" --success "$SUCCESS" || {
    echo "[post-search] Warning: post-search hook failed for $TOOL_NAME" >&2
}
