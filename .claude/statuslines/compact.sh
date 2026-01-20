#!/bin/bash
# Compact Statusline - Key metrics on one line
# Format: Model │ dir ⎇ branch │ ctx % │ $cost │ ● daemon

# Read input from stdin
INPUT=$(cat)

# Extract workspace info
CWD=$(echo "$INPUT" | jq -r '.workspace.current_dir // .cwd // "."' 2>/dev/null)
MODEL=$(echo "$INPUT" | jq -r '.model.display_name // "Claude"' 2>/dev/null)
[ -z "$CWD" ] && CWD="."
[ -z "$MODEL" ] && MODEL="Claude"

# Shorten model name for compact display
case "$MODEL" in
    "Opus 4.5"|"Claude Opus 4.5") MODEL="Opus" ;;
    "Sonnet 4.5"|"Claude Sonnet 4.5") MODEL="Sonnet" ;;
    "Haiku 4.5"|"Claude Haiku 4.5") MODEL="Haiku" ;;
    "Claude"*) MODEL=$(echo "$MODEL" | sed 's/Claude //' | cut -d' ' -f1) ;;
esac

# Git info
BRANCH=$(cd "$CWD" 2>/dev/null && git branch --show-current 2>/dev/null || echo "")
DIR=$(basename "$CWD")

# ANSI Colors
RST="\033[0m"; BOLD="\033[1m"; DIM="\033[2m"
GRN="\033[32m"; YEL="\033[33m"; CYN="\033[36m"; RED="\033[31m"

# Data sources
DAEMON_FILE="$CWD/.claude-flow/daemon-state.json"
INTEL_FILE="$CWD/.ruvector/intelligence.json"

# Daemon status
DAEMON_RUNNING="false"
WORKERS_SUCCESS=0
if [ -f "$DAEMON_FILE" ]; then
    DAEMON_RUNNING=$(jq -r '.running // false' "$DAEMON_FILE" 2>/dev/null || echo "false")
    WORKERS_SUCCESS=$(jq -r '[.workers[].successCount] | add // 0' "$DAEMON_FILE" 2>/dev/null || echo "0")
fi

# Pattern count
RV_PATTERNS=0
if [ -f "$INTEL_FILE" ]; then
    RV_PATTERNS=$(jq -r '.stats.total_patterns // (.patterns | length) // 0' "$INTEL_FILE" 2>/dev/null || echo "0")
fi

# Context percentage (from swarm data if available)
SWARM_JSON=$(timeout 0.5s npx @claude-flow/cli@latest hooks statusline --json 2>/dev/null || echo "{}")
CTX_PCT=$(echo "$SWARM_JSON" | jq -r '.context.percentage // 0' 2>/dev/null || echo "0")
COST=$(echo "$SWARM_JSON" | jq -r '.cost.session // "0.00"' 2>/dev/null || echo "0.00")

# Build compact line
printf "${BOLD}$MODEL${RST}"
printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
[ -n "$BRANCH" ] && printf " ${YEL}⎇${RST} ${BRANCH:0:15}"

# Context if significant
[ "$CTX_PCT" -gt 0 ] && printf " ${DIM}│${RST} ctx ${BOLD}${CTX_PCT}%%${RST}"

# Cost if any
[ "$COST" != "0.00" ] && [ "$COST" != "0" ] && printf " ${DIM}│${RST} ${GRN}\$${COST}${RST}"

# Daemon indicator
if [ "$DAEMON_RUNNING" = "true" ]; then
    printf " ${DIM}│${RST} ${GRN}●${RST}"
    [ "$WORKERS_SUCCESS" -gt 0 ] && printf " ${DIM}${WORKERS_SUCCESS}${RST}"
else
    printf " ${DIM}│ ○${RST}"
fi

# Patterns if any
[ "$RV_PATTERNS" -gt 0 ] && printf " ${DIM}│${RST} ${DIM}◆${RV_PATTERNS}${RST}"

echo
