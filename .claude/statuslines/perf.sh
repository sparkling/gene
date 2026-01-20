#!/bin/bash
# perf - Performance statusline (1 line)
# Category: Performance
# Format: Model │ HNSW ms │ Flash speedup │ memory saved

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_perf
load_intel
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    printf " ${DIM}│${RST} ⚡ HNSW ${HNSW_MS}ms"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} ⚡ ${HNSW_MS}ms"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${YEL}⚡${RST} HNSW ${BOLD}${HNSW_MS}ms${RST}"
    printf " ${DIM}│${RST} Flash ${BOLD}${FLASH_SPEEDUP}x${RST}"
    [ "$MEM_SAVED" -gt 0 ] && printf " ${DIM}│${RST} ${GRN}$MEM_SAVED%%${RST} saved"
    ;;
esac
echo
