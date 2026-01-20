#!/bin/bash
# daemon - Daemon status statusline (1 line)
# Category: Operations (Daemon/Workers)
# Format: Model │ daemon status │ runs │ workers

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_daemon
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf " ${DIM}│${RST} ${GRN}●${RST} daemon"
    fi
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} $SWARM_AGENTS agents"
    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf " ${DIM}│${RST} ${GRN}●${RST} daemon"
    fi
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"

    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf " ${DIM}│${RST} ${GRN}●${RST} daemon"
      [ "$WORKERS_SUCCESS" -gt 0 ] && printf " ${DIM}│${RST} ${BOLD}$WORKERS_SUCCESS${RST} runs"
      printf " ${DIM}│${RST} ${BOLD}$WORKERS_RUNNING${RST} workers"
    else
      printf " ${DIM}│${RST} ${RED}○${RST} daemon ${DIM}(stopped)${RST}"
    fi
    ;;
esac
echo
