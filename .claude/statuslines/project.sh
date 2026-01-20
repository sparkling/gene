#!/bin/bash
# project - Project status statusline (2 lines)
# Category: Project Progress
# Line 1: Domains
# Line 2: Codebase structure

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_progress
load_git
load_intel
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs"
    printf " ${DIM}│${RST} Project: $V3_PROGRESS%% complete"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡ SWARM${RST}: $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} Project: $V3_PROGRESS%%"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${CYN}📁${RST} Project ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf " ${YEL}⎇ $BRANCH${RST}"
    ;;
esac
echo

# Line 2: Progress details
bar=$(progress_bar $V3_PROGRESS)
printf "V3: $bar ${BOLD}$V3_PROGRESS%%${RST}"
printf " ${DIM}│${RST} Domains: ${BOLD}$V3_DOMAINS${RST}/$V3_TOTAL_DOMAINS"
printf " ${DIM}│${RST} DDD: ${BOLD}$DDD_PROGRESS%%${RST}"

# Learning stats
[ "$RV_PATTERNS" -gt 0 ] && printf " ${DIM}│${RST} ${MAG}◆${RST}$RV_PATTERNS"
echo
