#!/bin/bash
# full - Full-width AI engineer dashboard (3 lines, uses terminal width)
# Category: Dashboard (Comprehensive)
# Line 1: Model │ Project │ Git │ Mode
# Line 2: Context │ Cost │ Tokens │ Cache
# Line 3: Agents │ Tasks │ Memory │ Daemon
#
# This template shows the most actionable info for AI engineers:
# - Context window usage (critical for long conversations)
# - Cost tracking (budget awareness)
# - Git status (uncommitted work)
# - Agent/task status (swarm coordination)

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_git
load_cost
load_swarm
load_daemon
load_intel
load_vectors
load_security

MODE=$(detect_mode)
TERM_WIDTH=${COLUMNS:-120}

# ============== LINE 1: Identity + Git + Mode ==============
case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    printf " ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs ($CRITICAL_CVES critical)${RST}"
    ;;
  swarm)
    printf "${GRN}⬡ SWARM${RST} ${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    printf " ${DIM}│${RST} ${GRN}$REGISTRY_AGENTS${RST} agents"
    [ "$SWARM_TASKS" -gt 0 ] && printf " ${YEL}$SWARM_TASKS${RST} tasks"
    [ "$SWARM_COORD" = "true" ] && printf " ${GRN}●${RST}hive"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${CYN}$DIR${RST}"
    [ -n "$BRANCH" ] && printf "/${YEL}$BRANCH${RST}"
    if [ "$UNCOMMITTED" -gt 0 ]; then
      printf " ${DIM}│${RST} ${YEL}$UNCOMMITTED uncommitted${RST}"
    fi
    if [ "$ADDED" -gt 0 ] || [ "$REMOVED" -gt 0 ]; then
      [ "$ADDED" -gt 0 ] && printf " ${GRN}+$ADDED${RST}"
      [ "$REMOVED" -gt 0 ] && printf " ${RED}-$REMOVED${RST}"
    fi
    [ "$STASH_COUNT" -gt 0 ] && printf " ${DIM}($STASH_COUNT stashed)${RST}"
    ;;
esac
echo

# ============== LINE 2: Context + Cost + Tokens ==============
# Context usage (CRITICAL for AI engineers)
if [ "$CTX_PCT" -gt 80 ]; then
  printf "${RED}ctx:${CTX_PCT}%%${RST}"
elif [ "$CTX_PCT" -gt 60 ]; then
  printf "${YEL}ctx:${CTX_PCT}%%${RST}"
else
  printf "ctx:${CTX_PCT}%%"
fi

# Cost tracking
printf " ${DIM}│${RST} ${GRN}\$$COST_SESSION${RST}"

# Token counts
printf " ${DIM}│${RST} $(fmt_num $TOKENS_IN)↓ $(fmt_num $TOKENS_OUT)↑"

# Cache efficiency (if available)
if [ "$CACHE_READ" -gt 0 ] || [ "$CACHE_WRITE" -gt 0 ]; then
  printf " ${DIM}│${RST} cache:$(fmt_num $CACHE_READ)r/$(fmt_num $CACHE_WRITE)w"
fi

# Model info for context
printf " ${DIM}│${RST} $MODEL"
echo

# ============== LINE 3: Agents + Memory + Daemon ==============
# Show this line only if there's meaningful data
HAS_LINE3="false"

# Agents: show RUNNING (Task tool) if any, otherwise show REGISTRY (MCP)
if [ "$TASK_AGENTS" -gt 0 ]; then
  # Show actual running Task tool agents
  printf "${GRN}⬡ $TASK_AGENTS running${RST}"
  [ "$REGISTRY_AGENTS" -gt 0 ] && printf " ${DIM}($REGISTRY_AGENTS registered)${RST}"
  [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS queued${RST}"
  printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
  HAS_LINE3="true"
elif [ "$REGISTRY_AGENTS" -gt 0 ] || [ "$SWARM_TASKS" -gt 0 ]; then
  # No running agents, show registry count
  printf "⬡ $REGISTRY_AGENTS registered"
  [ "$SWARM_TASKS" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}$SWARM_TASKS queued${RST}"
  printf " ${DIM}│${RST} $SWARM_TOPOLOGY"
  HAS_LINE3="true"
fi

# Memory stats
if [ "$DB_EXISTS" = "true" ]; then
  [ "$HAS_LINE3" = "true" ] && printf " ${DIM}│${RST}"
  printf " 💾 $TOTAL_DB_SIZE"
  HAS_LINE3="true"
fi

# Patterns/Learning
if [ "$RV_PATTERNS" -gt 0 ]; then
  [ "$HAS_LINE3" = "true" ] && printf " ${DIM}│${RST}"
  printf " ${MAG}◆${RST}$RV_PATTERNS"
  HAS_LINE3="true"
fi

# Daemon status
if [ "$DAEMON_RUNNING" = "true" ]; then
  [ "$HAS_LINE3" = "true" ] && printf " ${DIM}│${RST}"
  printf " ${GRN}●${RST}daemon"
  [ "$WORKERS_RUNNING" -gt 0 ] && printf " ${YEL}$WORKERS_RUNNING${RST} active"
  [ -n "$DAEMON_UPTIME" ] && printf " ${DIM}up:$DAEMON_UPTIME${RST}"
  HAS_LINE3="true"
fi

[ "$HAS_LINE3" = "true" ] && echo
