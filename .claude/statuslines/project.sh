#!/bin/bash
# project - Project overview statusline (2 lines)
# Category: Developer
# Line 1: Model │ Project │ Branch
# Line 2: Git stats │ File changes │ Stash

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
load_cost
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡ SWARM${RST}: $REGISTRY_AGENTS agents"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${CYN}📁${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " ${YEL}⎇ $BRANCH${RST}"
    ;;
esac
echo

# Line 2: Git activity
if [ "$UNCOMMITTED" -gt 0 ] || [ "$ADDED" -gt 0 ] || [ "$REMOVED" -gt 0 ]; then
  printf "Changes:"
  [ "$UNCOMMITTED" -gt 0 ] && printf " ${YEL}$UNCOMMITTED files${RST}"
  [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
  [ "$REMOVED" -gt 0 ] && printf " ${RED}-$REMOVED${RST}"
else
  printf "${DIM}No uncommitted changes${RST}"
fi

[ "$STASH_COUNT" -gt 0 ] && printf " ${DIM}│${RST} $STASH_COUNT stashed"

# Show cost if meaningful
if [ "$COST_SESSION" != "0.00" ]; then
  printf " ${DIM}│${RST} ${GRN}\$$COST_SESSION${RST}"
fi

# Last commit summary if available
if [ -n "$LAST_COMMIT" ]; then
  printf " ${DIM}│${RST} ${DIM}$LAST_COMMIT${RST}"
fi
echo
