#!/bin/bash
# workers - Worker details statusline (2-3 lines)
# Category: Operations (Daemon/Workers)
# Line 1: Model + daemon
# Line 2-3: Per-worker stats

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_daemon
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs"
    [ "$CRITICAL_CVES" -gt 0 ] && printf " (${RED}$CRITICAL_CVES critical${RST})"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡ SWARM${RST}: $SWARM_AGENTS agents"
    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf " ${DIM}│${RST} ${GRN}●${RST} daemon"
    fi
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} "
    if [ "$DAEMON_RUNNING" = "true" ]; then
      printf "${GRN}●${RST} Daemon Active"
      printf " ${DIM}│${RST} ${BOLD}$WORKERS_SUCCESS${RST}/${WORKERS_TOTAL} runs"
    else
      printf "${RED}○${RST} Daemon Stopped"
    fi
    ;;
esac
echo

# Line 2+: Worker details
if [ "$DAEMON_RUNNING" = "true" ] && [ -f "$DAEMON_FILE" ]; then
  # Get worker stats as compact list
  WORKERS=$(jq -r '.workers | to_entries[] | "\(.key):\(.value.runCount)/\(.value.successCount)"' "$DAEMON_FILE" 2>/dev/null | head -6)

  if [ -n "$WORKERS" ]; then
    printf "${DIM}Workers:${RST} "
    echo "$WORKERS" | while read worker; do
      name=$(echo "$worker" | cut -d: -f1)
      stats=$(echo "$worker" | cut -d: -f2)
      runs=$(echo "$stats" | cut -d/ -f1)
      success=$(echo "$stats" | cut -d/ -f2)

      if [ "$runs" -gt 0 ]; then
        if [ "$success" = "$runs" ]; then
          printf "${GRN}$name${RST}:$stats "
        else
          printf "${YEL}$name${RST}:$stats "
        fi
      else
        printf "${DIM}$name:0${RST} "
      fi
    done
    echo
  fi

  # Active workers if any
  ACTIVE=$(jq -r '.workers | to_entries[] | select(.value.running == true) | .key' "$DAEMON_FILE" 2>/dev/null | tr '\n' ' ')
  [ -n "$ACTIVE" ] && printf "${GRN}Running:${RST} $ACTIVE\n"
else
  printf "${DIM}No worker activity${RST}\n"
fi
