#!/bin/bash
# Post-edit hook: Record editing outcome for learning
# Receives JSON via stdin with tool_name, tool_input, and tool_response

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract file_path from tool_input
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

# Check if tool_response exists (indicates success)
HAS_RESPONSE=$(echo "$INPUT" | jq -r 'if .tool_response then "true" else "false" end' 2>/dev/null)

if [ -z "$FILE_PATH" ]; then
    exit 0
fi

# Determine success based on presence of tool_response
SUCCESS="true"
if [ "$HAS_RESPONSE" = "false" ]; then
    SUCCESS="false"
fi

# Log for debugging
# echo "[post-edit] Recording: $FILE_PATH (success=$SUCCESS)" >&2

# Call claude-flow hooks post-edit
cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."
npx @claude-flow/cli@latest hooks post-edit --file "$FILE_PATH" --success "$SUCCESS" || {
    echo "[post-edit] Warning: post-edit hook failed for $FILE_PATH" >&2
}

# Store successful edit pattern in memory
if [ "$SUCCESS" = "true" ]; then
    # Extract file extension for pattern categorization
    EXT="${FILE_PATH##*.}"
    TIMESTAMP=$(date +%s)

    npx @claude-flow/cli@latest memory store \
        --namespace patterns \
        --key "edit-${EXT}-${TIMESTAMP}" \
        --value "Edited file: $FILE_PATH" 2>/dev/null || true
fi
