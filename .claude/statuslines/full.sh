#!/bin/bash
# Full Statusline - All data visible (5-6 lines)
# Shows: identity, mode, daemon, intelligence, database, metrics

set -o pipefail

# Read input from stdin
INPUT=$(cat)

# Extract workspace info
CWD=$(echo "$INPUT" | jq -r '.workspace.current_dir // .cwd // "."' 2>/dev/null)
MODEL=$(echo "$INPUT" | jq -r '.model.display_name // "Claude"' 2>/dev/null)
[ -z "$CWD" ] && CWD="."
[ -z "$MODEL" ] && MODEL="Claude"

# Git info
BRANCH=$(cd "$CWD" 2>/dev/null && git branch --show-current 2>/dev/null || echo "")
GIT_STATUS=$(cd "$CWD" 2>/dev/null && git status --porcelain 2>/dev/null | head -1)
DIR=$(basename "$CWD")

# =============================================================================
# ANSI Colors
# =============================================================================
RST="\033[0m"; BOLD="\033[1m"; DIM="\033[2m"
RED="\033[31m"; GRN="\033[32m"; YEL="\033[33m"
BLU="\033[34m"; MAG="\033[35m"; CYN="\033[36m"

# =============================================================================
# Data Sources
# =============================================================================
DAEMON_FILE="$CWD/.claude-flow/daemon-state.json"
INTEL_FILE="$CWD/.ruvector/intelligence.json"
USER_DB="$CWD/data/user.db"
OPS_DB="$CWD/data/operational.db"
PROGRESS_FILE="$CWD/.claude-flow/metrics/v3-progress.json"

# =============================================================================
# Helper: Safe JSON extraction
# =============================================================================
jq_safe() {
    local file="$1" query="$2" default="$3"
    [ -f "$file" ] && result=$(jq -r "$query" "$file" 2>/dev/null) && [ "$result" != "null" ] && [ -n "$result" ] && echo "$result" || echo "$default"
}

# =============================================================================
# Extract All Data
# =============================================================================

# --- Daemon/Worker Data ---
DAEMON_RUNNING="false"
WORKERS_TOTAL=0
WORKERS_SUCCESS=0
if [ -f "$DAEMON_FILE" ]; then
    DAEMON_RUNNING=$(jq_safe "$DAEMON_FILE" '.running' "false")
    WORKERS_TOTAL=$(jq_safe "$DAEMON_FILE" '[.workers[].runCount] | add // 0' "0")
    WORKERS_SUCCESS=$(jq_safe "$DAEMON_FILE" '[.workers[].successCount] | add // 0' "0")
fi

# --- Swarm Data (from CLI with timeout) ---
SWARM_JSON=$(timeout 1s npx @claude-flow/cli@latest hooks statusline --json 2>/dev/null || echo "{}")
SWARM_AGENTS=$(echo "$SWARM_JSON" | jq -r '.swarm.activeAgents // 0' 2>/dev/null || echo "0")
SWARM_MAX=$(echo "$SWARM_JSON" | jq -r '.swarm.maxAgents // 15' 2>/dev/null || echo "15")
CTX_PCT=$(echo "$SWARM_JSON" | jq -r '.context.percentage // 0' 2>/dev/null || echo "0")
COST=$(echo "$SWARM_JSON" | jq -r '.cost.session // "0.00"' 2>/dev/null || echo "0.00")

# --- RuVector Intelligence Data ---
RV_PATTERNS=0; RV_MEMORIES=0; RV_TRAJECTORIES=0; RV_ERRORS=0
RV_BEST_ALGO="none"; RV_CONVERGENCE=0; RV_SESSIONS=0
if [ -f "$INTEL_FILE" ]; then
    RV_PATTERNS=$(jq_safe "$INTEL_FILE" '.stats.total_patterns // (.patterns | length) // 0' "0")
    RV_MEMORIES=$(jq_safe "$INTEL_FILE" '.stats.total_memories // (.memories | length) // 0' "0")
    RV_TRAJECTORIES=$(jq_safe "$INTEL_FILE" '.stats.total_trajectories // (.trajectories | length) // 0' "0")
    RV_ERRORS=$(jq_safe "$INTEL_FILE" '.stats.total_errors // 0' "0")
    RV_SESSIONS=$(jq_safe "$INTEL_FILE" '.stats.session_count // 0' "0")
    RV_BEST_ALGO=$(jq -r '.learning.stats // {} | to_entries | map(select(.value.updates > 0)) | sort_by(-.value.updates) | .[0].key // "none"' "$INTEL_FILE" 2>/dev/null || echo "none")
    if [ "$RV_BEST_ALGO" != "none" ]; then
        conv=$(jq_safe "$INTEL_FILE" ".learning.stats.\"$RV_BEST_ALGO\".convergenceScore // 0" "0")
        RV_CONVERGENCE=$(echo "$conv * 100" | bc 2>/dev/null | cut -d. -f1 || echo "0")
        [ -z "$RV_CONVERGENCE" ] && RV_CONVERGENCE=0
    fi
fi

# --- Database Stats ---
USER_VECTORS=0; USER_SIZE="0K"; OPS_VECTORS=0; OPS_SIZE="0K"
if [ -f "$USER_DB" ]; then
    USER_SIZE=$(du -h "$USER_DB" 2>/dev/null | cut -f1)
    bytes=$(stat -c%s "$USER_DB" 2>/dev/null || echo "0")
    USER_VECTORS=$((bytes / 1536))
fi
if [ -f "$OPS_DB" ]; then
    OPS_SIZE=$(du -h "$OPS_DB" 2>/dev/null | cut -f1)
    bytes=$(stat -c%s "$OPS_DB" 2>/dev/null || echo "0")
    OPS_VECTORS=$((bytes / 768))
fi

# --- V3 Progress ---
V3_PROGRESS=0
if [ -f "$PROGRESS_FILE" ]; then
    V3_PROGRESS=$(jq_safe "$PROGRESS_FILE" '.overall.percentage // 0' "0")
fi

# --- Git file changes ---
ADDED=0; REMOVED=0
if [ -n "$GIT_STATUS" ]; then
    DIFF_STAT=$(cd "$CWD" 2>/dev/null && git diff --shortstat 2>/dev/null || echo "")
    ADDED=$(echo "$DIFF_STAT" | grep -oE '[0-9]+ insertion' | grep -oE '[0-9]+' || echo "0")
    REMOVED=$(echo "$DIFF_STAT" | grep -oE '[0-9]+ deletion' | grep -oE '[0-9]+' || echo "0")
    [ -z "$ADDED" ] && ADDED=0
    [ -z "$REMOVED" ] && REMOVED=0
fi

# =============================================================================
# LINE 1: Identity + Git
# =============================================================================
printf "${BOLD}$MODEL${RST} in ${CYN}$DIR${RST}"
[ -n "$BRANCH" ] && printf " on ${YEL}⎇ $BRANCH${RST}"
echo

# =============================================================================
# LINE 2: Mode + Daemon + Security
# =============================================================================
# Detect mode
MODE="idle"
if [ "$SWARM_AGENTS" -gt 1 ]; then
    MODE="swarm"
elif [ "$RV_TRAJECTORIES" -gt 0 ]; then
    MODE="learning"
elif [ "$USER_VECTORS" -gt 0 ]; then
    MODE="database"
fi

case "$MODE" in
    swarm)    printf "${GRN}⬡${RST} swarm" ;;
    learning) printf "${MAG}🧠${RST} learning" ;;
    database) printf "${BLU}💾${RST} database" ;;
    idle)     printf "${DIM}⬡ idle${RST}" ;;
esac

# Daemon
if [ "$DAEMON_RUNNING" = "true" ]; then
    printf " ${DIM}│${RST} ${GRN}●${RST} daemon"
    [ "$WORKERS_SUCCESS" -gt 0 ] && printf " ${DIM}(${WORKERS_SUCCESS} runs)${RST}"
else
    printf " ${DIM}│ ○ daemon${RST}"
fi

# Security indicator
if [ -z "$GIT_STATUS" ]; then
    printf " ${DIM}│${RST} ${GRN}CLEAN${RST}"
else
    printf " ${DIM}│${RST} ${YEL}MODIFIED${RST}"
fi
echo

# =============================================================================
# LINE 3: Intelligence
# =============================================================================
printf "${MAG}🧠${RST} "
printf "${BOLD}$RV_PATTERNS${RST} patterns"
printf " ${DIM}│${RST} ${BOLD}$RV_MEMORIES${RST} memories"
printf " ${DIM}│${RST} ${BOLD}$RV_TRAJECTORIES${RST} trajectories"
[ "$RV_ERRORS" -gt 0 ] && printf " ${DIM}│${RST} ${RED}⚠${RV_ERRORS}${RST}"
echo

# =============================================================================
# LINE 4: Database
# =============================================================================
printf "${BLU}💾${RST} "
printf "User: ${BOLD}$USER_VECTORS${RST} vectors"
[ "$USER_SIZE" != "0K" ] && printf " ($USER_SIZE)"
printf " ${DIM}│${RST} Ops: ${BOLD}$OPS_VECTORS${RST} vectors"
[ "$OPS_SIZE" != "0K" ] && printf " ($OPS_SIZE)"
echo

# =============================================================================
# LINE 5: Metrics
# =============================================================================
printf "${DIM}💰${RST} \$${COST}"
printf " ${DIM}│${RST} ctx ${BOLD}${CTX_PCT}%%${RST}"

# File changes
if [ "$ADDED" -gt 0 ] || [ "$REMOVED" -gt 0 ]; then
    printf " ${DIM}│${RST}"
    [ "$ADDED" -gt 0 ] && printf " ${GRN}+${ADDED}${RST}"
    [ "$REMOVED" -gt 0 ] && printf " ${RED}-${REMOVED}${RST}"
fi

# V3 Progress
if [ "$V3_PROGRESS" -gt 0 ]; then
    printf " ${DIM}│${RST} V3 ${BOLD}${V3_PROGRESS}%%${RST}"
fi
echo
