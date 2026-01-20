#!/bin/bash
# git - Git-detailed statusline (2 lines)
# Category: Developer Daily (Code-Centric)
# Line 1: Model + dir + branch
# Line 2: Changes, stash, last commit

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
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
    printf "${GRN}⬡ SWARM${RST}: $SWARM_AGENTS/$SWARM_MAX agents │ $SWARM_TOPOLOGY"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} in ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " on ${YEL}⎇ $BRANCH${RST}"
    ;;
esac
echo

# Line 2: Git details
if [ -n "$GIT_STATUS" ]; then
  # Changes
  printf "${DIM}Changes:${RST}"
  [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
  [ "$REMOVED" -gt 0 ] && printf " ${RED}-$REMOVED${RST}"
  printf " ${DIM}│${RST} ${YEL}$UNCOMMITTED${RST} files"

  # Stash
  [ "$STASH_COUNT" -gt 0 ] && printf " ${DIM}│${RST} ${MAG}⊡$STASH_COUNT stash${RST}"

  # Last commit (truncated)
  [ -n "$LAST_COMMIT" ] && printf " ${DIM}│ $LAST_COMMIT${RST}"
else
  printf "${DIM}Working tree clean${RST}"
  [ -n "$LAST_COMMIT" ] && printf " ${DIM}│ $LAST_COMMIT${RST}"
fi
echo
