#!/bin/bash
# focus - Minimalist statusline: Model + branch (1 line)
# Category: Minimalist (Focus Mode)
# Purpose: Distraction-free coding with git awareness

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠${RST} ${BOLD}$MODEL_SHORT${RST}"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${YEL}⎇ $BRANCH${RST}"
    printf " ${RED}$TOTAL_CVES CVEs${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST}"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${YEL}⎇ $BRANCH${RST}"
    printf " ${DIM}($SWARM_AGENTS)${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${YEL}⎇ $BRANCH${RST}"
    ;;
esac
echo
