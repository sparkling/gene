#!/bin/bash
# Adaptive Claude-Flow V3 + RuVector Statusline
# Shows different information based on current activity mode
# Modes: swarm | learning | database | idle

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
SWARM_COORDINATION=$(echo "$SWARM_JSON" | jq -r '.swarm.coordinationActive // false' 2>/dev/null || echo "false")

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
    USER_VECTORS=$((bytes / 1536))  # ~768 dims * 2 bytes estimate
fi
if [ -f "$OPS_DB" ]; then
    OPS_SIZE=$(du -h "$OPS_DB" 2>/dev/null | cut -f1)
    bytes=$(stat -c%s "$OPS_DB" 2>/dev/null || echo "0")
    OPS_VECTORS=$((bytes / 768))  # ~384 dims * 2 bytes estimate
fi

# --- V3 Progress ---
V3_PROGRESS=0
if [ -f "$PROGRESS_FILE" ]; then
    V3_PROGRESS=$(jq_safe "$PROGRESS_FILE" '.overall.percentage // 0' "0")
fi

# =============================================================================
# Mode Detection
# =============================================================================
# Priority: swarm > learning > database > idle

MODE="idle"
if [ "$SWARM_AGENTS" -gt 1 ]; then
    MODE="swarm"
elif [ "$RV_TRAJECTORIES" -gt 0 ] || [ "$RV_CONVERGENCE" -gt 0 ]; then
    MODE="learning"
elif [ "$USER_VECTORS" -gt 0 ] || [ "$OPS_VECTORS" -gt 0 ]; then
    MODE="database"
fi

# =============================================================================
# LINE 1: Identity (always shown)
# =============================================================================
printf "${BOLD}$MODEL${RST} in ${CYN}$DIR${RST}"
[ -n "$BRANCH" ] && printf " on ${YEL}⎇ $BRANCH${RST}"
echo

# =============================================================================
# LINE 2+: Mode-specific content
# =============================================================================

case "$MODE" in
    swarm)
        # SWARM MODE: Show agents, coordination, tasks
        printf "${GRN}⬡${RST} ${BOLD}Swarm Active${RST}: ${GRN}$SWARM_AGENTS${RST}/$SWARM_MAX agents"

        # Daemon status
        if [ "$DAEMON_RUNNING" = "true" ]; then
            printf " ${DIM}│${RST} ${GRN}●${RST} Daemon"
        else
            printf " ${DIM}│${RST} ${RED}○${RST} Daemon"
        fi

        # Workers if any
        if [ "$WORKERS_TOTAL" -gt 0 ]; then
            printf " ${DIM}│${RST} ⚙ ${WORKERS_SUCCESS}/${WORKERS_TOTAL} workers"
        fi
        echo

        # Third line: Intelligence summary
        printf "${MAG}🧠${RST} "
        [ "$RV_PATTERNS" -gt 0 ] && printf "${MAG}◆${RST}${RV_PATTERNS} patterns " || printf "${DIM}◇ patterns${RST} "
        [ "$RV_MEMORIES" -gt 0 ] && printf "${BLU}⬡${RST}${RV_MEMORIES} memories" || printf "${DIM}⬡ memories${RST}"
        echo
        ;;

    learning)
        # LEARNING MODE: Show patterns, trajectories, convergence, algorithm
        printf "${MAG}🧠${RST} ${BOLD}Learning Active${RST}"

        # Daemon
        if [ "$DAEMON_RUNNING" = "true" ]; then
            printf " ${DIM}│${RST} ${GRN}●${RST} Daemon"
        fi
        echo

        # Intelligence details
        printf "${MAG}◆${RST} ${BOLD}$RV_PATTERNS${RST} patterns"
        printf " ${DIM}│${RST} ${BLU}⬡${RST} ${BOLD}$RV_MEMORIES${RST} memories"
        printf " ${DIM}│${RST} ${YEL}↝${RST} ${BOLD}$RV_TRAJECTORIES${RST} trajectories"
        echo

        # Algorithm + convergence (only show if algorithm exists)
        if [ "$RV_BEST_ALGO" != "none" ] && [ -n "$RV_BEST_ALGO" ]; then
            filled=$((RV_CONVERGENCE / 20))
            empty=$((5 - filled))
            bar=""
            for ((i=0; i<filled; i++)); do bar+="▰"; done
            for ((i=0; i<empty; i++)); do bar+="▱"; done
            printf "${CYN}$RV_BEST_ALGO${RST} $bar ${BOLD}${RV_CONVERGENCE}%%${RST} converged"
            [ "$RV_ERRORS" -gt 0 ] && printf " ${DIM}│${RST} ${RED}⚠${RV_ERRORS}${RST} errors"
            echo
        fi
        [ "$RV_SESSIONS" -gt 0 ] && printf "${DIM}#${RV_SESSIONS} sessions${RST}\n"
        ;;

    database)
        # DATABASE MODE: Show vector counts, sizes, both databases
        printf "${BLU}💾${RST} ${BOLD}Database Active${RST}"

        # Daemon
        if [ "$DAEMON_RUNNING" = "true" ]; then
            printf " ${DIM}│${RST} ${GRN}●${RST} Daemon"
        fi
        echo

        # User database (critical)
        if [ "$USER_VECTORS" -gt 0 ]; then
            printf "${GRN}◆${RST} User: ${BOLD}$USER_VECTORS${RST} vectors ($USER_SIZE)"
        else
            printf "${DIM}◇ User: empty${RST}"
        fi

        printf " ${DIM}│${RST} "

        # Operational database
        if [ "$OPS_VECTORS" -gt 0 ]; then
            printf "${CYN}⚙${RST} Ops: ${BOLD}$OPS_VECTORS${RST} vectors ($OPS_SIZE)"
        else
            printf "${DIM}⚙ Ops: empty${RST}"
        fi
        echo

        # V3 Progress if available
        if [ "$V3_PROGRESS" -gt 0 ]; then
            filled=$((V3_PROGRESS / 20))
            empty=$((5 - filled))
            bar=""
            for ((i=0; i<filled; i++)); do bar+="▰"; done
            for ((i=0; i<empty; i++)); do bar+="▱"; done
            printf "V3: $bar ${BOLD}${V3_PROGRESS}%%${RST}"

            # Add patterns if any
            [ "$RV_PATTERNS" -gt 0 ] && printf " ${DIM}│${RST} ${MAG}◆${RST}${RV_PATTERNS} patterns"
            echo
        fi
        ;;

    idle)
        # IDLE MODE: Minimal info
        printf "${DIM}⬡ Swarm: idle${RST}"

        # Daemon
        if [ "$DAEMON_RUNNING" = "true" ]; then
            printf " ${DIM}│${RST} ${GRN}●${RST} Daemon ${DIM}(${WORKERS_SUCCESS} runs)${RST}"
        else
            printf " ${DIM}│${RST} ${DIM}○ Daemon off${RST}"
        fi
        echo

        # Quick intelligence summary if any data exists
        if [ "$RV_PATTERNS" -gt 0 ] || [ "$USER_VECTORS" -gt 0 ]; then
            printf "${MAG}🧠${RST} "
            [ "$RV_PATTERNS" -gt 0 ] && printf "${DIM}◆${RV_PATTERNS}${RST} "
            [ "$USER_VECTORS" -gt 0 ] && printf "${DIM}💾${USER_VECTORS} vectors${RST}"
            echo
        fi
        ;;
esac
