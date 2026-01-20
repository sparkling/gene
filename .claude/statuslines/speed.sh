#!/bin/bash
# speed - Speed metrics statusline (2 lines)
# Category: Performance
# Line 1: Search/attention
# Line 2: Worker durations

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_perf
load_daemon
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs"
    printf " ${DIM}│${RST} Performance: HNSW ${HNSW_MS}ms"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡ SWARM${RST}: $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} ⚡ HNSW ${HNSW_MS}ms ${DIM}│${RST} Flash ${FLASH_SPEEDUP}x"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${YEL}⚡${RST} Performance"
    printf " ${DIM}│${RST} HNSW ${BOLD}${HNSW_MS}ms${RST}"
    printf " ${DIM}│${RST} Flash ${BOLD}${FLASH_SPEEDUP}x${RST}"
    ;;
esac
echo

# Line 2: Detailed metrics
printf "${DIM}Metrics:${RST}"
printf " Cache ${BOLD}$CACHE_HIT%%${RST} hit"
[ "$MEM_SAVED" -gt 0 ] && printf " ${DIM}│${RST} Memory ${GRN}-$MEM_SAVED%%${RST}"

# Worker timing if daemon running
if [ "$DAEMON_RUNNING" = "true" ] && [ -f "$DAEMON_FILE" ]; then
  AVG_DURATION=$(jq -r '[.workers[].avgDuration // 0] | add / length | floor' "$DAEMON_FILE" 2>/dev/null || echo "0")
  [ "$AVG_DURATION" -gt 0 ] && printf " ${DIM}│${RST} Workers avg ${BOLD}${AVG_DURATION}ms${RST}"
fi
echo
