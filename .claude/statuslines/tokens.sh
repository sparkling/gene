#!/bin/bash
# tokens - Token breakdown statusline (2 lines)
# Category: Cost & Efficiency
# Line 1: Cost + context
# Line 2: Cache stats, efficiency

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_cost
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs"
    printf " ${DIM}│${RST} Session: ${GRN}\$$COST_SESSION${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡ SWARM${RST}: $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} ${GRN}\$$COST_SESSION${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${GRN}💰${RST} Session: ${BOLD}\$$COST_SESSION${RST}"
    printf " ${DIM}│${RST} Context: ${BOLD}$CTX_PCT%%${RST}"
    ;;
esac
echo

# Line 2: Token breakdown
printf "${DIM}Tokens:${RST}"
printf " ${CYN}↓${RST}$(fmt_num $TOKENS_IN)"
printf " ${MAG}↑${RST}$(fmt_num $TOKENS_OUT)"

# Cache if any
if [ "$CACHE_WRITE" -gt 0 ] || [ "$CACHE_READ" -gt 0 ]; then
  printf " ${DIM}│${RST} Cache:"
  [ "$CACHE_WRITE" -gt 0 ] && printf " ${GRN}W$(fmt_num $CACHE_WRITE)${RST}"
  [ "$CACHE_READ" -gt 0 ] && printf " ${BLU}R$(fmt_num $CACHE_READ)${RST}"
fi

# Efficiency (cache read / total input)
if [ "$TOKENS_IN" -gt 0 ] && [ "$CACHE_READ" -gt 0 ]; then
  EFFICIENCY=$((CACHE_READ * 100 / TOKENS_IN))
  printf " ${DIM}│${RST} ${BOLD}$EFFICIENCY%%${RST} cached"
fi
echo
