#!/bin/bash
# learning - Learning metrics statusline (2 lines)
# Category: Knowledge Base (Vector/Learning)
# Line 1: Model + patterns
# Line 2: Trajectories, Q-values, sessions

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_intel
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs"
    [ "$RV_ERRORS" -gt 0 ] && printf " ${DIM}│${RST} ${RED}$RV_ERRORS errors${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡ SWARM${RST}: $SWARM_AGENTS agents ${DIM}│${RST} ${MAG}◆${RST}$RV_PATTERNS patterns"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${MAG}🧠${RST} Learning"
    printf " ${DIM}│${RST} ${MAG}◆${RST} ${BOLD}$RV_PATTERNS${RST} patterns"
    printf " ${DIM}│${RST} ${BLU}⬡${RST} ${BOLD}$RV_MEMORIES${RST} memories"
    ;;
esac
echo

# Line 2: Learning details
printf "${YEL}↝${RST} ${BOLD}$RV_TRAJECTORIES${RST} trajectories"

if [ "$RV_BEST_ALGO" != "none" ] && [ -n "$RV_BEST_ALGO" ]; then
  printf " ${DIM}│${RST} ${CYN}$RV_BEST_ALGO${RST}"
  bar=$(progress_bar $RV_CONVERGENCE)
  printf " $bar ${BOLD}${RV_CONVERGENCE}%%${RST}"
fi

[ "$RV_SESSIONS" -gt 0 ] && printf " ${DIM}│${RST} #$RV_SESSIONS sessions"
[ "$RV_ERRORS" -gt 0 ] && printf " ${DIM}│${RST} ${RED}⚠$RV_ERRORS${RST}"
echo
