#!/bin/bash
# Pre-edit hook: Get context before editing files
# Receives JSON via stdin with tool_name and tool_input

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract file_path from tool_input using jq
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

if [ -z "$FILE_PATH" ]; then
    # No file path provided, skip silently
    exit 0
fi

# Log for debugging (comment out in production)
# echo "[pre-edit] Processing: $FILE_PATH" >&2

# Call claude-flow hooks pre-edit
cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."
npx @claude-flow/cli@latest hooks pre-edit --file "$FILE_PATH" || {
    echo "[pre-edit] Warning: pre-edit hook failed for $FILE_PATH" >&2
}
