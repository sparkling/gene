#!/bin/bash
# progress - Project progress statusline (1 line)
# Category: Project Progress
# Format: Model │ domains │ DDD % │ sessions

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_progress
load_intel
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    printf " ${DIM}│${RST} 📈 $V3_DOMAINS/$V3_TOTAL_DOMAINS domains"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} 📈 $V3_PROGRESS%%"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${CYN}📈${RST} $V3_DOMAINS/$V3_TOTAL_DOMAINS domains"
    printf " ${DIM}│${RST} DDD ${BOLD}$DDD_PROGRESS%%${RST}"
    [ "$RV_SESSIONS" -gt 0 ] && printf " ${DIM}│${RST} $RV_SESSIONS sessions"
    ;;
esac
echo
