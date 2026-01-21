#!/bin/bash
# swarm - Swarm status statusline (full width, 1 line)
# Category: Agents & Coordination
# Shows: Model │ Registry agents │ Active tasks │ Topology │ Git status

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_swarm
load_security
load_git

MODE=$(detect_mode)

# Get terminal width for full-width layout
TERM_WIDTH=${COLUMNS:-120}

case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    printf " ${DIM}│${RST} ⬡ $REGISTRY_AGENTS registered"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    ;;
  swarm)
    # Active swarm with tasks
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${GRN}$REGISTRY_AGENTS${RST} registered"
    [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS tasks${RST}"
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    [ "$SWARM_COORD" = "true" ] && printf " ${DIM}│${RST} ${GRN}●${RST}hive"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    [ "$UNCOMMITTED" -gt 0 ] && printf " ${DIM}+$UNCOMMITTED${RST}"
    ;;
  *)
    # Normal mode - show model and git info prominently
    printf "${BOLD}$MODEL_SHORT${RST}"
    # Show running agents if any, otherwise registry count
    if [ "$TASK_AGENTS" -gt 0 ]; then
      printf " ${DIM}│${RST} ${GRN}⬡ $TASK_AGENTS running${RST}"
    elif [ "$REGISTRY_AGENTS" -gt 0 ]; then
      printf " ${DIM}│${RST} ⬡ $REGISTRY_AGENTS registered"
    fi
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    [ "$UNCOMMITTED" -gt 0 ] && printf " ${DIM}+$UNCOMMITTED${RST}"
    [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
    [ "$REMOVED" -gt 0 ] && printf " ${RED}-$REMOVED${RST}"
    ;;
esac
echo
