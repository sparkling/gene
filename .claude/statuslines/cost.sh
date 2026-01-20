#!/bin/bash
# cost - Cost summary statusline (1 line)
# Category: Cost & Efficiency
# Format: Model │ $cost │ tokens in │ tokens out │ ctx %

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_cost
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    printf " ${DIM}│${RST} ${GRN}\$$COST_SESSION${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} ${GRN}\$$COST_SESSION${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${GRN}\$$COST_SESSION${RST}"
    printf " ${DIM}│${RST} $(fmt_num $TOKENS_IN) in"
    printf " ${DIM}│${RST} $(fmt_num $TOKENS_OUT) out"
    [ "$CTX_PCT" -gt 0 ] && printf " ${DIM}│${RST} ${BOLD}$CTX_PCT%%${RST} ctx"
    ;;
esac
echo
