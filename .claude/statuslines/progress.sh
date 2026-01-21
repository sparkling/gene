#!/bin/bash
# progress - Project activity statusline (1 line)
# Category: Project Status
# Shows: Model │ Git activity │ Session stats │ Files changed

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
load_cost
load_intel
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    [ "$UNCOMMITTED" -gt 0 ] && printf " ${DIM}+$UNCOMMITTED files${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} $REGISTRY_AGENTS agents"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    [ "$UNCOMMITTED" -gt 0 ] && printf " ${DIM}+$UNCOMMITTED${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    [ "$UNCOMMITTED" -gt 0 ] && printf " +$UNCOMMITTED"
    [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
    [ "$REMOVED" -gt 0 ] && printf " ${RED}-$REMOVED${RST}"
    [ "$RV_SESSIONS" -gt 0 ] && printf " ${DIM}│${RST} $RV_SESSIONS sessions"
    [ "$RV_PATTERNS" -gt 0 ] && printf " ${DIM}│${RST} $RV_PATTERNS patterns"
    ;;
esac
echo
