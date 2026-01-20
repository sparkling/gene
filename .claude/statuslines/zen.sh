#!/bin/bash
# zen - Minimalist statusline: Model only (1 line)
# Category: Minimalist (Focus Mode)
# Purpose: Distraction-free coding, minimal cognitive load

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠${RST} ${BOLD}$MODEL_SHORT${RST}"
    [ "$CRITICAL_CVES" -gt 0 ] && printf " ${RED}$CRITICAL_CVES CVEs${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}($SWARM_AGENTS agents)${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    ;;
esac
echo
