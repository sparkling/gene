#!/bin/bash
# dashboard - Comprehensive 3-line dashboard
# Category: Dashboard (Comprehensive)
# Line 1: Identity + mode
# Line 2: Cost + context
# Line 3: Vectors + learning

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
load_progress

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
    printf " ${DIM}│${RST} ${GRN}$SWARM_AGENTS${RST}/$SWARM_MAX"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} in ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " on ${YEL}⎇ $BRANCH${RST}"

    # Daemon indicator
    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf " ${DIM}│${RST} ${GRN}●${RST}"
    else
      printf " ${DIM}│ ○${RST}"
    fi
    ;;
esac
echo

# Line 2: Cost + Context + Git
printf "${GRN}💰${RST} \$$COST_SESSION"
printf " ${DIM}│${RST} ctx ${BOLD}$CTX_PCT%%${RST}"
printf " ${DIM}│${RST} $(fmt_num $TOKENS_IN)↓ $(fmt_num $TOKENS_OUT)↑"

# Git changes
if [ "$ADDED" -gt 0 ] || [ "$REMOVED" -gt 0 ]; then
  printf " ${DIM}│${RST}"
  [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
  [ "$REMOVED" -gt 0 ] && printf "${RED}-$REMOVED${RST}"
fi
echo

# Line 3: Vectors + Learning
printf "${BLU}📊${RST} $(fmt_num $USER_VECTORS) user"
printf " ${DIM}│${RST} $(fmt_num $OPS_VECTORS) ops"
printf " ${DIM}│${RST} ${MAG}◆${RST}$RV_PATTERNS patterns"
printf " ${DIM}│${RST} ${YEL}↝${RST}$RV_TRAJECTORIES traj"

# V3 Progress if available
if [ "$V3_PROGRESS" -gt 0 ]; then
  printf " ${DIM}│${RST} V3 ${BOLD}$V3_PROGRESS%%${RST}"
fi
echo
