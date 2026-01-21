#!/bin/bash
# agents - Agent registry details (full width, 2-3 lines)
# Category: Agents & Coordination
# Line 1: Status + registry count + topology
# Line 2: Agent types registered
# Line 3: Coordination status (if active)

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_swarm
load_security
load_git
load_cost

MODE=$(detect_mode)

# Line 1: Status header
case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST}: $TOTAL_CVES CVEs detected"
    [ "$REGISTRY_AGENTS" -gt 0 ] && printf " ${DIM}│${RST} $REGISTRY_AGENTS agents in registry"
    ;;
  swarm)
    printf "${GRN}⬡ SWARM${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${GRN}$REGISTRY_AGENTS${RST}/$SWARM_MAX in registry"
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS tasks queued${RST}"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ⬡ $REGISTRY_AGENTS/$SWARM_MAX"
    printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
    [ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${CYN}$BRANCH${RST}"
    [ "$UNCOMMITTED" -gt 0 ] && printf " ${DIM}+$UNCOMMITTED${RST}"
    ;;
esac
echo

# Line 2: Registered agent types
if [ "$REGISTRY_AGENTS" -gt 0 ]; then
  printf "${DIM}Registry:${RST}"

  if [ -n "$ACTIVE_AGENT_LIST" ]; then
    for agent in $ACTIVE_AGENT_LIST; do
      case "$agent" in
        researcher)  printf " ${MAG}🔍${RST}researcher" ;;
        coder)       printf " ${GRN}💻${RST}coder" ;;
        tester)      printf " ${YEL}🧪${RST}tester" ;;
        reviewer)    printf " ${CYN}👀${RST}reviewer" ;;
        architect)   printf " ${BLU}🏗️${RST}architect" ;;
        security*)   printf " ${RED}🔒${RST}security" ;;
        analyst)     printf " ${MAG}📊${RST}analyst" ;;
        planner)     printf " ${BLU}📋${RST}planner" ;;
        *)           printf " ${DIM}●${RST}$agent" ;;
      esac
    done
  else
    printf " ${DIM}(types not available)${RST}"
  fi

  # Show context usage if available
  if [ "$CTX_PCT" -gt 0 ]; then
    printf " ${DIM}│${RST} ctx:${CTX_PCT}%%"
  fi
else
  printf "${DIM}No agents registered │ Use: agent spawn -t coder${RST}"
fi
echo

# Line 3: Coordination status and tasks (only if meaningful)
if [ "$SWARM_COORD" = "true" ]; then
  printf "${GRN}●${RST} Hive-mind active"
  printf " ${DIM}│${RST} Status: $SWARM_STATUS"
  [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS pending tasks${RST}"
  echo
elif [ "$SWARM_TASKS" -gt 0 ]; then
  printf "${YEL}○${RST} $SWARM_TASKS tasks queued ${DIM}│${RST} No coordination active\n"
fi
