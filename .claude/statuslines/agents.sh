#!/bin/bash
# agents - Agent details statusline (2-3 lines)
# Category: Swarm (Multi-Agent)
# Line 1: Topology
# Line 2-3: Active agents

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_swarm
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    load_security
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs"
    [ "$SWARM_AGENTS" -gt 0 ] && printf " ${DIM}│${RST} $SWARM_AGENTS agents active"
    ;;
  swarm)
    printf "${GRN}⬡ SWARM ACTIVE${RST}"
    printf " ${DIM}│${RST} ${GRN}${BOLD}$SWARM_AGENTS${RST}/$SWARM_MAX agents"
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    [ "$SWARM_COORD" = "true" ] && printf " ${DIM}│${RST} ${GRN}●${RST} coordinated"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${CYN}⬡${RST} Swarm"
    printf " ${DIM}│${RST} $SWARM_AGENTS/$SWARM_MAX agents"
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    ;;
esac
echo

# Line 2: Active agents
if [ "$SWARM_AGENTS" -gt 0 ]; then
  printf "${DIM}Active:${RST}"

  if [ -n "$ACTIVE_AGENT_LIST" ]; then
    # Show agent types with icons
    for agent in $ACTIVE_AGENT_LIST; do
      case "$agent" in
        researcher)  printf " ${MAG}🔍${RST}$agent" ;;
        coder)       printf " ${GRN}💻${RST}$agent" ;;
        tester)      printf " ${YEL}🧪${RST}$agent" ;;
        reviewer)    printf " ${CYN}👀${RST}$agent" ;;
        architect)   printf " ${BLU}🏗️${RST}$agent" ;;
        security*)   printf " ${RED}🔒${RST}$agent" ;;
        *)           printf " ${DIM}●${RST}$agent" ;;
      esac
    done
  else
    printf " ${DIM}(no details)${RST}"
  fi

  [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS tasks${RST}"
else
  printf "${DIM}No active agents │ Swarm idle${RST}"
fi
echo

# Line 3: Coordination status (only if active)
if [ "$SWARM_AGENTS" -gt 1 ] && [ "$SWARM_COORD" = "true" ]; then
  printf "${GRN}✓${RST} Coordination active ${DIM}│${RST} Topology: $SWARM_TOPOLOGY\n"
fi
