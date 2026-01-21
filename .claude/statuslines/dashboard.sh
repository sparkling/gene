#!/bin/bash
# dashboard - Full-width comprehensive dashboard (3 lines)
# Category: Dashboard (Comprehensive)
# Line 1: Model │ Project │ Branch │ Mode indicator
# Line 2: Cost │ Context │ Tokens │ Git changes
# Line 3: Memory │ Patterns │ Daemon status

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
load_daemon
load_intel
load_vectors
load_cost
load_security
load_swarm

MODE=$(detect_mode)

# Line 1: Identity + Mode
case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} "
    printf "${BOLD}$MODEL_SHORT${RST} in ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " ${YEL}⎇ $BRANCH${RST}"
    printf " ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    ;;
  swarm)
    printf "${GRN}⬡ SWARM${RST} ${DIM}│${RST} "
    printf "${BOLD}$MODEL_SHORT${RST} in ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " ${YEL}⎇ $BRANCH${RST}"
    printf " ${DIM}│${RST} ${GRN}$REGISTRY_AGENTS${RST}/$SWARM_MAX registered"
    [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS tasks${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} in ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " on ${YEL}⎇ $BRANCH${RST}"
    # Daemon indicator
    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf " ${DIM}│${RST} ${GRN}●${RST} daemon"
    else
      printf " ${DIM}│${RST} ${DIM}○${RST}"
    fi
    ;;
esac
echo

# Line 2: Cost + Context + Tokens + Git
printf "${GRN}💰${RST} \$$COST_SESSION"
printf " ${DIM}│${RST} ctx ${BOLD}$CTX_PCT%%${RST}"
printf " ${DIM}│${RST} $(fmt_num $TOKENS_IN)↓ $(fmt_num $TOKENS_OUT)↑"

# Git changes
if [ "$UNCOMMITTED" -gt 0 ]; then
  printf " ${DIM}│${RST} $UNCOMMITTED uncommitted"
fi
if [ "$ADDED" -gt 0 ] || [ "$REMOVED" -gt 0 ]; then
  printf " ${DIM}│${RST}"
  [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
  [ "$REMOVED" -gt 0 ] && printf "${RED}-$REMOVED${RST}"
fi
echo

# Line 3: Memory + Learning (only show if data exists)
if [ "$DB_EXISTS" = "true" ] || [ "$RV_PATTERNS" -gt 0 ] || [ "$WORKERS_TOTAL" -gt 0 ]; then
  [ "$DB_EXISTS" = "true" ] && printf "${BLU}💾${RST} $TOTAL_DB_SIZE"
  [ "$RV_PATTERNS" -gt 0 ] && printf " ${DIM}│${RST} ${MAG}◆${RST}$RV_PATTERNS patterns"
  [ "$RV_TRAJECTORIES" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}↝${RST}$RV_TRAJECTORIES trajectories"
  [ "$WORKERS_TOTAL" -gt 0 ] && printf " ${DIM}│${RST} $WORKERS_SUCCESS/$WORKERS_TOTAL workers"
  echo
fi
