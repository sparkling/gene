#!/bin/bash
# Pre-task hook: Record task start and get agent suggestions
# Receives JSON via stdin with tool_input.prompt

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract prompt and tool_use_id
PROMPT=$(echo "$INPUT" | jq -r '.tool_input.prompt // empty' 2>/dev/null)
TOOL_USE_ID=$(echo "$INPUT" | jq -r '.tool_use_id // empty' 2>/dev/null)

if [ -z "$PROMPT" ]; then
    exit 0
fi

# Generate task ID from tool_use_id or timestamp
if [ -n "$TOOL_USE_ID" ]; then
    TASK_ID="task-${TOOL_USE_ID}"
else
    TASK_ID="task-$(date +%s)"
fi

# Log for debugging
# echo "[pre-task] Starting: $TASK_ID - ${PROMPT:0:60}..." >&2

# Call claude-flow hooks pre-task
cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."
npx @claude-flow/cli@latest hooks pre-task --task-id "$TASK_ID" --description "$PROMPT" || {
    echo "[pre-task] Warning: pre-task hook failed for $TASK_ID" >&2
}
