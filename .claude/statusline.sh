#!/bin/bash
# Claude Flow V3 + RuVector Live Statusline

INPUT=$(cat)
MODEL=$(echo "$INPUT" | jq -r '.model.display_name // "Claude"')
CWD=$(echo "$INPUT" | jq -r '.workspace.current_dir // .cwd')
DIR=$(basename "$CWD")
BRANCH=$(cd "$CWD" 2>/dev/null && git branch --show-current 2>/dev/null)

# Colors
RESET="\033[0m"; BOLD="\033[1m"; CYAN="\033[36m"; YELLOW="\033[33m"
GREEN="\033[32m"; MAGENTA="\033[35m"; BLUE="\033[34m"; DIM="\033[2m"; RED="\033[31m"

# Line 1: Model and location
printf "$BOLD$MODEL$RESET in $CYAN$DIR$RESET"
[ -n "$BRANCH" ] && printf " on $YELLOW⎇ $BRANCH$RESET"
echo

# Get live stats from claude-flow (with timeout to prevent blocking)
STATS=$(timeout 2s npx @claude-flow/cli@latest hooks statusline --json 2>/dev/null)

if [ -n "$STATS" ]; then
  # Parse live data
  AGENTS=$(echo "$STATS" | jq -r '.swarm.activeAgents // 0')
  MAX_AGENTS=$(echo "$STATS" | jq -r '.swarm.maxAgents // 50')
  PATTERNS=$(echo "$STATS" | jq -r '.v3Progress.patternsLearned // 0')
  SESSIONS=$(echo "$STATS" | jq -r '.v3Progress.sessionsCompleted // 0')
  MEMORY_ENTRIES=$(echo "$STATS" | jq -r '.memory.entries // 0')
  MEMORY_SIZE=$(echo "$STATS" | jq -r '.memory.sizeMb // 0')
  SEC_STATUS=$(echo "$STATS" | jq -r '.security.status // "OK"')
  ROUTING_SAVINGS=$(echo "$STATS" | jq -r '.routing.costSavings // "0%"')

  # Line 2: Swarm status (1 agent = daemon only, treat as idle)
  if [ "$AGENTS" -gt 1 ]; then
    printf "$GREEN⬡$RESET Swarm: $AGENTS/$MAX_AGENTS agents"
  else
    printf "$DIM⬡ Swarm: idle$RESET"
  fi

  # Memory
  if [ "$MEMORY_ENTRIES" -gt 0 ] 2>/dev/null; then
    printf " | $BLUE◆$RESET Memory: ${MEMORY_ENTRIES} entries"
  fi

  # Routing savings
  [ "$ROUTING_SAVINGS" != "0%" ] && printf " | $GREEN↓$RESET $ROUTING_SAVINGS saved"
  echo

  # Line 3: Intelligence
  printf "$MAGENTA🧠 RuVector$RESET"
  [ "$PATTERNS" -gt 0 ] 2>/dev/null && printf " $GREEN◆$RESET ${PATTERNS} patterns" || printf " $DIM◇ learning$RESET"
  [ "$SESSIONS" -gt 0 ] 2>/dev/null && printf " | $DIM#${SESSIONS} sessions$RESET"

  # Security status
  if [ "$SEC_STATUS" = "CLEAN" ] || [ "$SEC_STATUS" = "OK" ]; then
    printf " | $GREEN✓$RESET Secure"
  elif [ "$SEC_STATUS" = "PENDING" ]; then
    printf " | $YELLOW⚠$RESET Security: scan needed"
  else
    printf " | $RED✗$RESET Security: $SEC_STATUS"
  fi
  echo
else
  # Fallback to static file
  INTEL_FILE="$CWD/.ruvector/intelligence.json"
  if [ -f "$INTEL_FILE" ]; then
    PATTERN_COUNT=$(jq -r '.patterns | length // 0' "$INTEL_FILE" 2>/dev/null)
    MEMORY_COUNT=$(jq -r '.memories | length // 0' "$INTEL_FILE" 2>/dev/null)
    printf "$MAGENTA🧠 RuVector$RESET"
    [ "$PATTERN_COUNT" -gt 0 ] 2>/dev/null && printf " $GREEN◆$RESET $PATTERN_COUNT patterns"
    [ "$MEMORY_COUNT" -gt 0 ] 2>/dev/null && printf " $BLUE⬡$RESET $MEMORY_COUNT memories"
    echo
  else
    printf "$DIM🧠 RuVector: initializing...$RESET\n"
  fi
fi
