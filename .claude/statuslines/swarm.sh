#!/bin/bash
# swarm - Swarm summary statusline (1 line)
# Category: Swarm (Multi-Agent)
# Format: Model │ agents/max │ topology

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_swarm
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    printf " ${DIM}│${RST} ⬡ $SWARM_AGENTS/$SWARM_MAX"
    ;;
  swarm)
    printf "${GRN}⬡ SWARM${RST} ${DIM}│${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${GRN}${BOLD}$SWARM_AGENTS${RST}/$SWARM_MAX agents"
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS tasks${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ⬡ ${BOLD}$SWARM_AGENTS${RST}/$SWARM_MAX agents"
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    ;;
esac
echo
