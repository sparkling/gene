#!/bin/bash
# Post-task hook: Record task completion for learning
# Receives JSON via stdin with tool_input, tool_response, and tool_use_id

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract tool_use_id for task identification
TOOL_USE_ID=$(echo "$INPUT" | jq -r '.tool_use_id // empty' 2>/dev/null)
PROMPT=$(echo "$INPUT" | jq -r '.tool_input.prompt // empty' 2>/dev/null)

# Check if tool_response exists (indicates success)
HAS_RESPONSE=$(echo "$INPUT" | jq -r 'if .tool_response then "true" else "false" end' 2>/dev/null)

# Generate task ID from tool_use_id or timestamp
if [ -n "$TOOL_USE_ID" ]; then
    TASK_ID="task-${TOOL_USE_ID}"
else
    TASK_ID="task-$(date +%s)"
fi

SUCCESS="true"
if [ "$HAS_RESPONSE" = "false" ]; then
    SUCCESS="false"
fi

# Log for debugging
# echo "[post-task] Completed: $TASK_ID (success=$SUCCESS)" >&2

# Call claude-flow hooks post-task
cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."
npx @claude-flow/cli@latest hooks post-task --task-id "$TASK_ID" --success "$SUCCESS" || {
    echo "[post-task] Warning: post-task hook failed for $TASK_ID" >&2
}

# Store successful task pattern in memory
if [ "$SUCCESS" = "true" ] && [ -n "$PROMPT" ]; then
    TIMESTAMP=$(date +%s)

    npx @claude-flow/cli@latest memory store \
        --namespace patterns \
        --key "task-pattern-${TIMESTAMP}" \
        --value "Task completed: ${PROMPT:0:200}" 2>/dev/null || true
fi
