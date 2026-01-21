#!/bin/bash
# adaptive - Context-aware auto-switching statusline (2-4 lines)
# Category: Adaptive (Context-Aware)
# Auto-detects: alert > swarm > active > idle modes
# Shows relevant info based on current activity

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
load_daemon
load_swarm
load_security
load_cost

MODE=$(detect_mode)

# =============================================================================
# LINE 1: Identity (always shown)
# =============================================================================
case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    printf " ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    [ "$CRITICAL_CVES" -gt 0 ] && printf " (${RED}${BOLD}$CRITICAL_CVES critical${RST})"
    ;;
  swarm)
    printf "${GRN}⬡ SWARM${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    printf " ${DIM}│${RST} ${GRN}$TASK_AGENTS${RST} running"
    [ "$REGISTRY_AGENTS" -gt "$TASK_AGENTS" ] && printf " ${DIM}($REGISTRY_AGENTS registered)${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    ;;
esac
echo

# =============================================================================
# LINE 2: Activity context
# =============================================================================
case "$MODE" in
  alert)
    # Security details
    printf "Run: ${CYN}npx @claude-flow/cli@latest security scan${RST}"
    [ "$HIGH_CVES" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$HIGH_CVES high${RST}"
    echo
    ;;
  swarm)
    # Swarm coordination info
    printf "$SWARM_TOPOLOGY"
    [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS tasks${RST}"
    [ "$SWARM_COORD" = "true" ] && printf " ${DIM}│${RST} ${GRN}●${RST}hive"
    # Daemon status
    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf " ${DIM}│${RST} ${GRN}●${RST}daemon"
      [ "$WORKERS_RUNNING" -gt 0 ] && printf " ${DIM}($WORKERS_RUNNING workers)${RST}"
    fi
    echo
    ;;
  *)
    # Normal mode - git and cost info
    HAS_INFO="false"

    # Git changes
    if [ "$UNCOMMITTED" -gt 0 ]; then
      printf "${YEL}$UNCOMMITTED uncommitted${RST}"
      [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
      [ "$REMOVED" -gt 0 ] && printf " ${RED}-$REMOVED${RST}"
      HAS_INFO="true"
    fi

    # Cost if meaningful
    if [ "$COST_SESSION" != "0.00" ]; then
      [ "$HAS_INFO" = "true" ] && printf " ${DIM}│${RST}"
      printf " ${GRN}\$$COST_SESSION${RST}"
      HAS_INFO="true"
    fi

    # Daemon status
    if [ "$DAEMON_RUNNING" = "true" ]; then
      [ "$HAS_INFO" = "true" ] && printf " ${DIM}│${RST}"
      printf " ${GRN}●${RST}daemon"
      HAS_INFO="true"
    fi

    # Task agents if any but not in swarm mode
    if [ "$TASK_AGENTS" -gt 0 ] && [ "$TASK_AGENTS" -le 1 ]; then
      [ "$HAS_INFO" = "true" ] && printf " ${DIM}│${RST}"
      printf " ⬡ $TASK_AGENTS agent"
      HAS_INFO="true"
    fi

    [ "$HAS_INFO" = "true" ] && echo || echo "${DIM}Ready${RST}"
    ;;
esac
