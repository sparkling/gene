#!/bin/bash
# User prompt hook: Route task and handle slash commands
# Receives JSON via stdin with prompt field

set -euo pipefail

# Read JSON from stdin
INPUT=$(cat)

# Extract prompt
PROMPT=$(echo "$INPUT" | jq -r '.prompt // empty' 2>/dev/null)

if [ -z "$PROMPT" ]; then
    exit 0
fi

cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || cd "$(dirname "$0")/../.."

# Handle /slt slash commands
case "$PROMPT" in
    '/slt')
        bash .claude/slt.sh list 2>/dev/null || true
        exit 0
        ;;
    '/slt n'|'/slt next')
        bash .claude/slt.sh n 2>/dev/null || true
        exit 0
        ;;
    '/slt p'|'/slt prev')
        bash .claude/slt.sh p 2>/dev/null || true
        exit 0
        ;;
    '/slt '*)
        bash .claude/slt.sh "${PROMPT#/slt }" 2>/dev/null || true
        exit 0
        ;;
esac

# Route non-slash-command prompts through claude-flow
# Log for debugging
# echo "[user-prompt] Routing: ${PROMPT:0:60}..." >&2

npx @claude-flow/cli@latest hooks route --task "$PROMPT" || {
    echo "[user-prompt] Warning: route hook failed" >&2
}
