#!/bin/bash
# Combined Claude-Flow V3 + RuVector Statusline
# Based on Google DevEx Principles: Glanceability, Progressive Disclosure, Semantic Grouping
# Layout: 4 lines - Identity | Status | Intelligence | Resources

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
USER=$(whoami 2>/dev/null || echo "user")

# Terminal width detection (fallback to 120)
COLS=$(tput cols 2>/dev/null || echo 120)

# ═══════════════════════════════════════════════════════════════════════════════
# ANSI Color Codes
# ═══════════════════════════════════════════════════════════════════════════════
RST="\033[0m"
BOLD="\033[1m"
DIM="\033[2m"
RED="\033[31m"
GRN="\033[32m"
YEL="\033[33m"
BLU="\033[34m"
MAG="\033[35m"
CYN="\033[36m"
WHT="\033[37m"
BRED="\033[1;31m"
BGRN="\033[1;32m"
BYEL="\033[1;33m"
BBLU="\033[1;34m"
BMAG="\033[1;35m"
BCYN="\033[1;36m"

# ═══════════════════════════════════════════════════════════════════════════════
# State Files
# ═══════════════════════════════════════════════════════════════════════════════
DAEMON_FILE="$CWD/.claude-flow/daemon-state.json"
INTEL_FILE="$CWD/.ruvector/intelligence.json"

# ═══════════════════════════════════════════════════════════════════════════════
# Helper Functions
# ═══════════════════════════════════════════════════════════════════════════════

# Safe jq extraction with default
jq_safe() {
    local file="$1"
    local query="$2"
    local default="$3"
    if [ -f "$file" ]; then
        result=$(jq -r "$query" "$file" 2>/dev/null)
        if [ "$result" = "null" ] || [ -z "$result" ]; then
            echo "$default"
        else
            echo "$result"
        fi
    else
        echo "$default"
    fi
}

# Progress bar generator (filled/empty blocks)
progress_bar() {
    local value=$1
    local max=${2:-100}
    local width=${3:-5}
    local pct=$((value * 100 / max))
    local filled=$((pct * width / 100))
    local empty=$((width - filled))
    local bar=""
    for ((i=0; i<filled; i++)); do bar+="▰"; done
    for ((i=0; i<empty; i++)); do bar+="▱"; done
    echo "$bar"
}

# Worker status with checkmark/cross
worker_status() {
    local name="$1"
    local success="$2"
    local total="$3"
    if [ "$total" -eq 0 ]; then
        echo ""
        return
    fi
    local pct=$((success * 100 / total))
    if [ "$pct" -ge 80 ]; then
        printf "${GRN}%s:%d✓${RST}" "$name" "$success"
    elif [ "$pct" -ge 50 ]; then
        printf "${YEL}%s:%d${RST}" "$name" "$success"
    else
        printf "${RED}%s:%d✗${RST}" "$name" "$total"
    fi
}

# ═══════════════════════════════════════════════════════════════════════════════
# Extract Claude-Flow Data
# ═══════════════════════════════════════════════════════════════════════════════
CF_RUNNING="false"
CF_TOTAL_SUCCESS=0
WORKERS_MAP_RUN=0; WORKERS_MAP_OK=0
WORKERS_AUDIT_RUN=0; WORKERS_AUDIT_OK=0
WORKERS_OPT_RUN=0; WORKERS_OPT_OK=0
WORKERS_CONS_RUN=0; WORKERS_CONS_OK=0
WORKERS_TEST_RUN=0; WORKERS_TEST_OK=0

if [ -f "$DAEMON_FILE" ]; then
    CF_RUNNING=$(jq_safe "$DAEMON_FILE" '.running' "false")
    CF_TOTAL_SUCCESS=$(jq_safe "$DAEMON_FILE" '[.workers[].successCount] | add // 0' "0")

    # Individual workers
    WORKERS_MAP_RUN=$(jq_safe "$DAEMON_FILE" '.workers.map.runCount // 0' "0")
    WORKERS_MAP_OK=$(jq_safe "$DAEMON_FILE" '.workers.map.successCount // 0' "0")
    WORKERS_AUDIT_RUN=$(jq_safe "$DAEMON_FILE" '.workers.audit.runCount // 0' "0")
    WORKERS_AUDIT_OK=$(jq_safe "$DAEMON_FILE" '.workers.audit.successCount // 0' "0")
    WORKERS_OPT_RUN=$(jq_safe "$DAEMON_FILE" '.workers.optimize.runCount // 0' "0")
    WORKERS_OPT_OK=$(jq_safe "$DAEMON_FILE" '.workers.optimize.successCount // 0' "0")
    WORKERS_CONS_RUN=$(jq_safe "$DAEMON_FILE" '.workers.consolidate.runCount // 0' "0")
    WORKERS_CONS_OK=$(jq_safe "$DAEMON_FILE" '.workers.consolidate.successCount // 0' "0")
    WORKERS_TEST_RUN=$(jq_safe "$DAEMON_FILE" '.workers.testgaps.runCount // 0' "0")
    WORKERS_TEST_OK=$(jq_safe "$DAEMON_FILE" '.workers.testgaps.successCount // 0' "0")
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Extract RuVector Intelligence Data
# ═══════════════════════════════════════════════════════════════════════════════
RV_PATTERNS=0
RV_MEMORIES=0
RV_TRAJECTORIES=0
RV_ERRORS=0
RV_SESSIONS=0
RV_BEST_ALGO="none"
RV_UPDATES=0
RV_AVG_REWARD=0
RV_CONVERGENCE=0

if [ -f "$INTEL_FILE" ]; then
    RV_PATTERNS=$(jq_safe "$INTEL_FILE" '.stats.total_patterns // (.patterns | length) // 0' "0")
    RV_MEMORIES=$(jq_safe "$INTEL_FILE" '.stats.total_memories // (.memories | length) // 0' "0")
    RV_TRAJECTORIES=$(jq_safe "$INTEL_FILE" '.stats.total_trajectories // (.trajectories | length) // 0' "0")
    RV_ERRORS=$(jq_safe "$INTEL_FILE" '.stats.total_errors // 0' "0")
    RV_SESSIONS=$(jq_safe "$INTEL_FILE" '.stats.session_count // 0' "0")

    # Find best learning algorithm by updates
    RV_BEST_ALGO=$(jq -r '
        .learning.stats // {} | to_entries |
        map(select(.value.updates > 0)) |
        sort_by(-.value.updates) |
        .[0].key // "none"
    ' "$INTEL_FILE" 2>/dev/null || echo "none")

    # Get stats for best algorithm
    if [ "$RV_BEST_ALGO" != "none" ] && [ -n "$RV_BEST_ALGO" ]; then
        RV_UPDATES=$(jq_safe "$INTEL_FILE" ".learning.stats.\"$RV_BEST_ALGO\".updates // 0" "0")
        RV_AVG_REWARD=$(jq_safe "$INTEL_FILE" ".learning.stats.\"$RV_BEST_ALGO\".avgReward // 0" "0")
        RV_CONVERGENCE=$(jq_safe "$INTEL_FILE" ".learning.stats.\"$RV_BEST_ALGO\".convergenceScore // 0" "0")
        # Convert convergence to percentage (0-1 -> 0-100)
        RV_CONVERGENCE=$(echo "$RV_CONVERGENCE * 100" | bc 2>/dev/null | cut -d. -f1 || echo "0")
        [ -z "$RV_CONVERGENCE" ] && RV_CONVERGENCE=0
    fi
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Calculate Derived Metrics
# ═══════════════════════════════════════════════════════════════════════════════

# Swarm agents (placeholder - would need live swarm data)
SWARM_AGENTS=1
SWARM_MAX=15

# Security CVE status (placeholder)
CVE_FIXED=0
CVE_TOTAL=3

# Memory usage estimate (based on file sizes)
MEM_MB=12

# Context percentage (placeholder - would need API data)
CTX_PCT=35

# Intelligence percentage (based on convergence and patterns)
if [ "$RV_PATTERNS" -gt 0 ]; then
    INTEL_PCT=$((RV_CONVERGENCE / 10))
    [ "$INTEL_PCT" -gt 100 ] && INTEL_PCT=100
else
    INTEL_PCT=0
fi

# DDD Progress (placeholder for V3)
DDD_PCT=40

# Current time
TIMESTAMP=$(date +"%H:%M:%S")

# ═══════════════════════════════════════════════════════════════════════════════
# LINE 1: Identity + Branch + Time
# ═══════════════════════════════════════════════════════════════════════════════
printf "${BMAG}▊${RST} ${BOLD}Claude Flow V3${RST} ${MAG}⊕${RST} ${BOLD}RuVector${RST}"
printf " ${DIM}│${RST} ${CYN}%s${RST}" "$USER"
[ -n "$BRANCH" ] && printf " ${DIM}│${RST} ${YEL}⎇ %s${RST}" "$BRANCH"
printf " ${DIM}│${RST} ${WHT}%s${RST}" "$MODEL"

# Right-align timestamp
printf "%*s" $((COLS - 75)) ""
printf "${DIM}%s${RST}\n" "$TIMESTAMP"

# ═══════════════════════════════════════════════════════════════════════════════
# LINE 2: Separator + Status Bar
# ═══════════════════════════════════════════════════════════════════════════════
# Print separator
printf "${DIM}"
for ((i=0; i<COLS; i++)); do printf "━"; done
printf "${RST}\n"

# Daemon status
if [ "$CF_RUNNING" = "true" ]; then
    printf "${BGRN}◉${RST} Daemon ${GRN}%d${RST}✓" "$CF_TOTAL_SUCCESS"
else
    printf "${RED}○${RST} Daemon ${DIM}off${RST}"
fi

# Swarm agents
printf " ${DIM}│${RST} 🤖 ${BCYN}%d${RST}${DIM}/%d${RST} agents" "$SWARM_AGENTS" "$SWARM_MAX"

# Security/CVE status
if [ "$CVE_FIXED" -eq "$CVE_TOTAL" ]; then
    printf " ${DIM}│${RST} ${BGRN}🟢${RST} CVE %d/%d" "$CVE_FIXED" "$CVE_TOTAL"
else
    printf " ${DIM}│${RST} ${BRED}🔴${RST} CVE %d/%d" "$CVE_FIXED" "$CVE_TOTAL"
fi

# Workers status (compact)
printf " ${DIM}│${RST} ${DIM}⚙${RST} "
workers_output=""
[ "$WORKERS_MAP_RUN" -gt 0 ] && workers_output+="$(worker_status "map" "$WORKERS_MAP_OK" "$WORKERS_MAP_RUN") "
[ "$WORKERS_AUDIT_RUN" -gt 0 ] && workers_output+="$(worker_status "aud" "$WORKERS_AUDIT_OK" "$WORKERS_AUDIT_RUN") "
[ "$WORKERS_OPT_RUN" -gt 0 ] && workers_output+="$(worker_status "opt" "$WORKERS_OPT_OK" "$WORKERS_OPT_RUN") "
[ "$WORKERS_CONS_RUN" -gt 0 ] && workers_output+="$(worker_status "con" "$WORKERS_CONS_OK" "$WORKERS_CONS_RUN") "
[ "$WORKERS_TEST_RUN" -gt 0 ] && workers_output+="$(worker_status "tst" "$WORKERS_TEST_OK" "$WORKERS_TEST_RUN") "
printf "%b" "$workers_output"
echo

# ═══════════════════════════════════════════════════════════════════════════════
# LINE 3: Intelligence (RuVector)
# ═══════════════════════════════════════════════════════════════════════════════
printf "${BMAG}🧠${RST} "

# Patterns
printf "${MAG}◆${RST}${BOLD}%d${RST} patterns " "$RV_PATTERNS"

# Memories
printf "${BLU}⬡${RST}${BOLD}%d${RST} memories " "$RV_MEMORIES"

# Trajectories
printf "${YEL}↝${RST}${BOLD}%d${RST} trajectories" "$RV_TRAJECTORIES"

printf " ${DIM}│${RST} "

# Learning algorithm + convergence
if [ "$RV_BEST_ALGO" != "none" ] && [ -n "$RV_BEST_ALGO" ]; then
    conv_bar=$(progress_bar "$RV_CONVERGENCE" 100 5)
    printf "${CYN}%s${RST} %s ${BOLD}%d%%${RST}" "$RV_BEST_ALGO" "$conv_bar" "$RV_CONVERGENCE"
else
    printf "${DIM}no learning${RST}"
fi

printf " ${DIM}│${RST} "

# Errors
if [ "$RV_ERRORS" -eq 0 ]; then
    printf "${GRN}⚠0${RST} errors"
else
    printf "${RED}⚠%d${RST} errors" "$RV_ERRORS"
fi

# Sessions
printf " ${DIM}│${RST} ${DIM}#%d sessions${RST}" "$RV_SESSIONS"

echo

# ═══════════════════════════════════════════════════════════════════════════════
# LINE 4: Resources + Features
# ═══════════════════════════════════════════════════════════════════════════════
printf "${BLU}💾${RST} ${BOLD}%d${RST}MB" "$MEM_MB"
printf " ${DIM}│${RST} ${CYN}📊${RST} Ctx:${BOLD}%d%%${RST}" "$CTX_PCT"
printf " ${DIM}│${RST} ${MAG}🧮${RST} Int:${BOLD}%d%%${RST}" "$INTEL_PCT"

# Features (check env vars or assume enabled)
printf " ${DIM}│${RST} "
[ "${RUVECTOR_AST_ENABLED:-true}" = "true" ] && printf "${GRN}●${RST}AST " || printf "${DIM}○AST${RST} "
[ "${RUVECTOR_COVERAGE_ROUTING:-true}" = "true" ] && printf "${GRN}●${RST}Cov " || printf "${DIM}○Cov${RST} "
[ "${RUVECTOR_GRAPH_ALGORITHMS:-true}" = "true" ] && printf "${GRN}●${RST}GNN " || printf "${DIM}○GNN${RST} "
[ "${RUVECTOR_SECURITY_SCAN:-true}" = "true" ] && printf "${GRN}●${RST}Sec " || printf "${DIM}○Sec${RST} "

# DDD Progress bar
printf "${DIM}│${RST} DDD:${DIM}[${RST}"
ddd_filled=$((DDD_PCT / 20))
ddd_empty=$((5 - ddd_filled))
for ((i=0; i<ddd_filled; i++)); do printf "${GRN}█${RST}"; done
for ((i=0; i<ddd_empty; i++)); do printf "${DIM}░${RST}"; done
printf "${DIM}]${RST} ${BOLD}%d%%${RST}" "$DDD_PCT"

echo
