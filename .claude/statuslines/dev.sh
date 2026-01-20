#!/bin/bash
# dev - Developer daily statusline (1 line)
# Category: Developer Daily (Code-Centric)
# Format: Model │ dir ⎇ branch │ +lines/-lines │ uncommitted

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${BOLD}$MODEL_SHORT${RST}"
    [ -n "$BRANCH" ] && printf " ${YEL}⎇ $BRANCH${RST}"
    printf " ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " ${YEL}⎇ $BRANCH${RST}"
    printf " ${DIM}│${RST} ${GRN}$SWARM_AGENTS${RST} agents"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " ${YEL}⎇ $BRANCH${RST}"

    # Show diff stats if any changes
    if [ "$ADDED" -gt 0 ] || [ "$REMOVED" -gt 0 ]; then
      printf " ${DIM}│${RST}"
      [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
      [ "$REMOVED" -gt 0 ] && printf "${RED}-$REMOVED${RST}"
    fi

    # Uncommitted files
    [ "$UNCOMMITTED" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$UNCOMMITTED uncommitted${RST}"
    ;;
esac
echo
