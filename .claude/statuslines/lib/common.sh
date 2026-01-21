#!/bin/bash
# Common functions for all statusline templates
# Source this file: source "$(dirname "$0")/lib/common.sh"

# =============================================================================
# ANSI Colors
# =============================================================================
setup_colors() {
  RST="\033[0m"; BOLD="\033[1m"; DIM="\033[2m"
  RED="\033[31m"; GRN="\033[32m"; YEL="\033[33m"
  BLU="\033[34m"; MAG="\033[35m"; CYN="\033[36m"
}

# =============================================================================
# Data File Paths
# =============================================================================
setup_paths() {
  DAEMON_FILE="$CWD/.claude-flow/daemon-state.json"
  INTEL_FILE="$CWD/.ruvector/intelligence.json"
  SECURITY_FILE="$CWD/.claude-flow/metrics/audit-status.json"
  PROGRESS_FILE="$CWD/.claude-flow/metrics/v3-progress.json"
  PERF_FILE="$CWD/.claude-flow/metrics/performance.json"
  CONFIG_FILE="$CWD/claude-flow.config.json"
  USER_DB="$CWD/data/user.db"
  OPS_DB="$CWD/data/operational.db"
  MEMORY_DB="$CWD/data/memories.db"
}

# =============================================================================
# Safe JSON Extraction
# =============================================================================
jq_safe() {
  local file="$1" query="$2" default="$3"
  [ -f "$file" ] && result=$(jq -r "$query" "$file" 2>/dev/null) && [ "$result" != "null" ] && [ -n "$result" ] && echo "$result" || echo "$default"
}

# =============================================================================
# Task Tool Agent Detection (REAL running agents)
# =============================================================================
count_task_agents() {
  # Real Task tool agents are tracked in /tmp/claude/{project}/tasks/*.output
  # These symlinks point to actual JSONL conversation files
  local project_name=$(basename "$CWD")
  local task_dir="/tmp/claude/-home-claude-src-${project_name}/tasks"

  if [ -d "$task_dir" ]; then
    ls "$task_dir"/*.output 2>/dev/null | wc -l
  else
    echo "0"
  fi
}

# =============================================================================
# Mode Detection (Alert > Swarm > Normal)
# =============================================================================
detect_mode() {
  # Alert mode: CVEs or errors or security issues
  local cves=$(jq_safe "$SECURITY_FILE" '.totalCves // 0' "0")
  local critical=$(jq_safe "$SECURITY_FILE" '.critical // 0' "0")
  local errors=$(jq_safe "$INTEL_FILE" '.stats.total_errors // 0' "0")

  [ "$cves" -gt 0 ] || [ "$critical" -gt 0 ] || [ "$errors" -gt 5 ] && echo "alert" && return

  # Swarm mode: REAL Task tool agents > 1 (not MCP registry)
  local task_agents=$(count_task_agents)
  [ "$task_agents" -gt 1 ] && echo "swarm" && return

  echo "normal"
}

# =============================================================================
# Input Processing
# =============================================================================
load_stdin() {
  INPUT=$(cat)
  CWD=$(echo "$INPUT" | jq -r '.workspace.current_dir // .cwd // "."' 2>/dev/null)
  MODEL=$(echo "$INPUT" | jq -r '.model.display_name // "Claude"' 2>/dev/null)
  [ -z "$CWD" ] && CWD="."
  [ -z "$MODEL" ] && MODEL="Claude"
}

# Shorten model name
shorten_model() {
  case "$MODEL" in
    "Opus 4.5"|"Claude Opus 4.5") echo "Opus" ;;
    "Sonnet 4.5"|"Claude Sonnet 4.5") echo "Sonnet" ;;
    "Sonnet 4"|"Claude Sonnet 4") echo "Sonnet4" ;;
    "Haiku 4.5"|"Claude Haiku 4.5") echo "Haiku" ;;
    "Claude"*) echo "$MODEL" | sed 's/Claude //' | cut -d' ' -f1 ;;
    *) echo "$MODEL" ;;
  esac
}

# =============================================================================
# Git Information
# =============================================================================
load_git() {
  BRANCH=$(cd "$CWD" 2>/dev/null && git branch --show-current 2>/dev/null || echo "")
  DIR=$(basename "$CWD")
  GIT_STATUS=$(cd "$CWD" 2>/dev/null && git status --porcelain 2>/dev/null | head -20)

  # Diff stats
  if [ -n "$GIT_STATUS" ]; then
    DIFF_STAT=$(cd "$CWD" 2>/dev/null && git diff --shortstat 2>/dev/null || echo "")
    ADDED=$(echo "$DIFF_STAT" | grep -oE '[0-9]+ insertion' | grep -oE '[0-9]+' || echo "0")
    REMOVED=$(echo "$DIFF_STAT" | grep -oE '[0-9]+ deletion' | grep -oE '[0-9]+' || echo "0")
    [ -z "$ADDED" ] && ADDED=0
    [ -z "$REMOVED" ] && REMOVED=0
  else
    ADDED=0
    REMOVED=0
  fi

  # Count uncommitted files
  UNCOMMITTED=$(echo "$GIT_STATUS" | grep -c '^' || echo "0")

  # Stash count
  STASH_COUNT=$(cd "$CWD" 2>/dev/null && git stash list 2>/dev/null | wc -l || echo "0")

  # Last commit (short)
  LAST_COMMIT=$(cd "$CWD" 2>/dev/null && git log -1 --pretty=format:'%h %s' 2>/dev/null | cut -c1-40 || echo "")
}

# =============================================================================
# Daemon/Worker Data
# =============================================================================
load_daemon() {
  DAEMON_RUNNING="false"
  DAEMON_UPTIME=""
  WORKERS_TOTAL=0
  WORKERS_SUCCESS=0
  WORKERS_RUNNING=0
  WORKER_LIST=""

  if [ -f "$DAEMON_FILE" ]; then
    # FIX: JSON uses .isRunning not .running
    DAEMON_RUNNING=$(jq_safe "$DAEMON_FILE" '.isRunning // .running' "false")
    DAEMON_UPTIME=$(jq_safe "$DAEMON_FILE" '.uptime // ""' "")
    WORKERS_TOTAL=$(jq_safe "$DAEMON_FILE" '[.workers[].runCount] | add // 0' "0")
    WORKERS_SUCCESS=$(jq_safe "$DAEMON_FILE" '[.workers[].successCount] | add // 0' "0")
    # FIX: Workers also use .isRunning
    WORKERS_RUNNING=$(jq_safe "$DAEMON_FILE" '[.workers[] | select(.isRunning == true or .running == true)] | length' "0")
    WORKER_LIST=$(jq -r '.workers | to_entries | map("\(.key):\(.value.runCount)") | join(" ")' "$DAEMON_FILE" 2>/dev/null || echo "")
  fi
}

# =============================================================================
# Intelligence/Learning Data
# =============================================================================
load_intel() {
  RV_PATTERNS=0; RV_MEMORIES=0; RV_TRAJECTORIES=0; RV_ERRORS=0
  RV_BEST_ALGO="none"; RV_CONVERGENCE=0; RV_SESSIONS=0
  RV_Q_AVG=0; RV_UPDATES=0

  if [ -f "$INTEL_FILE" ]; then
    RV_PATTERNS=$(jq_safe "$INTEL_FILE" '.stats.total_patterns // (.patterns | length) // 0' "0")
    RV_MEMORIES=$(jq_safe "$INTEL_FILE" '.stats.total_memories // (.memories | length) // 0' "0")
    RV_TRAJECTORIES=$(jq_safe "$INTEL_FILE" '.stats.total_trajectories // (.trajectories | length) // 0' "0")
    RV_ERRORS=$(jq_safe "$INTEL_FILE" '.stats.total_errors // 0' "0")
    RV_SESSIONS=$(jq_safe "$INTEL_FILE" '.stats.session_count // 0' "0")
    RV_UPDATES=$(jq_safe "$INTEL_FILE" '.stats.total_updates // 0' "0")

    # Best algorithm with convergence
    RV_BEST_ALGO=$(jq -r '.learning.stats // {} | to_entries | map(select(.value.updates > 0)) | sort_by(-.value.updates) | .[0].key // "none"' "$INTEL_FILE" 2>/dev/null || echo "none")
    if [ "$RV_BEST_ALGO" != "none" ]; then
      conv=$(jq_safe "$INTEL_FILE" ".learning.stats.\"$RV_BEST_ALGO\".convergenceScore // 0" "0")
      RV_CONVERGENCE=$(echo "$conv * 100" | bc 2>/dev/null | cut -d. -f1 || echo "0")
      [ -z "$RV_CONVERGENCE" ] && RV_CONVERGENCE=0

      RV_Q_AVG=$(jq_safe "$INTEL_FILE" ".learning.stats.\"$RV_BEST_ALGO\".avgQValue // 0" "0")
    fi
  fi
}

# =============================================================================
# Vector Database Stats (ACCURATE - query actual DB or show file size only)
# =============================================================================
load_vectors() {
  USER_SIZE="0K"; OPS_SIZE="0K"; MEMORY_SIZE="0K"
  DB_EXISTS="false"

  # Only show file sizes - vector counts require actual DB queries
  # which are too slow for statusline (removed fake calculations)
  if [ -f "$USER_DB" ]; then
    USER_SIZE=$(du -h "$USER_DB" 2>/dev/null | cut -f1)
    DB_EXISTS="true"
  fi

  if [ -f "$OPS_DB" ]; then
    OPS_SIZE=$(du -h "$OPS_DB" 2>/dev/null | cut -f1)
    DB_EXISTS="true"
  fi

  if [ -f "$MEMORY_DB" ]; then
    MEMORY_SIZE=$(du -h "$MEMORY_DB" 2>/dev/null | cut -f1)
    DB_EXISTS="true"
  fi

  # Total DB size (accurate metric)
  TOTAL_DB_SIZE="0K"
  if [ "$DB_EXISTS" = "true" ]; then
    local total_bytes=0
    [ -f "$USER_DB" ] && total_bytes=$((total_bytes + $(stat -c%s "$USER_DB" 2>/dev/null || echo 0)))
    [ -f "$OPS_DB" ] && total_bytes=$((total_bytes + $(stat -c%s "$OPS_DB" 2>/dev/null || echo 0)))
    [ -f "$MEMORY_DB" ] && total_bytes=$((total_bytes + $(stat -c%s "$MEMORY_DB" 2>/dev/null || echo 0)))
    if [ "$total_bytes" -ge 1073741824 ]; then
      TOTAL_DB_SIZE=$(echo "scale=1; $total_bytes / 1073741824" | bc 2>/dev/null)"G"
    elif [ "$total_bytes" -ge 1048576 ]; then
      TOTAL_DB_SIZE=$(echo "scale=1; $total_bytes / 1048576" | bc 2>/dev/null)"M"
    elif [ "$total_bytes" -ge 1024 ]; then
      TOTAL_DB_SIZE=$(echo "scale=0; $total_bytes / 1024" | bc 2>/dev/null)"K"
    else
      TOTAL_DB_SIZE="${total_bytes}B"
    fi
  fi
}

# =============================================================================
# Security Data
# =============================================================================
load_security() {
  TOTAL_CVES=0; CRITICAL_CVES=0; HIGH_CVES=0; MEDIUM_CVES=0
  LAST_SCAN="never"; SECURITY_STATUS="unknown"

  if [ -f "$SECURITY_FILE" ]; then
    TOTAL_CVES=$(jq_safe "$SECURITY_FILE" '.totalCves // 0' "0")
    CRITICAL_CVES=$(jq_safe "$SECURITY_FILE" '.critical // 0' "0")
    HIGH_CVES=$(jq_safe "$SECURITY_FILE" '.high // 0' "0")
    MEDIUM_CVES=$(jq_safe "$SECURITY_FILE" '.medium // 0' "0")
    LAST_SCAN=$(jq_safe "$SECURITY_FILE" '.lastScan // "never"' "never")

    if [ "$TOTAL_CVES" -eq 0 ]; then
      SECURITY_STATUS="CLEAN"
    elif [ "$CRITICAL_CVES" -gt 0 ]; then
      SECURITY_STATUS="CRITICAL"
    elif [ "$HIGH_CVES" -gt 0 ]; then
      SECURITY_STATUS="HIGH"
    else
      SECURITY_STATUS="WARN"
    fi
  fi
}

# =============================================================================
# Swarm Data (ACCURATE - distinguishes registry vs running)
# =============================================================================
load_swarm() {
  # REAL running Task tool agents (from /tmp/claude/.../tasks/)
  TASK_AGENTS=0
  # Registry agents (MCP store - shows registered, not necessarily running)
  REGISTRY_AGENTS=0
  # Config values
  SWARM_MAX=15; SWARM_TOPOLOGY="hierarchical"
  # Task tracking
  SWARM_TASKS=0; SWARM_COORD="false"
  ACTIVE_AGENT_LIST=""
  # Status indicators
  SWARM_STATUS="idle"  # idle, active, coordinated

  local AGENTS_STORE="$CWD/.claude-flow/agents/store.json"
  local TASKS_FILE="$CWD/.claude-flow/tasks.json"
  local HIVE_FILE="$CWD/.claude-flow/hive-mind/state.json"

  # Config for topology/max
  if [ -f "$CONFIG_FILE" ]; then
    SWARM_TOPOLOGY=$(jq_safe "$CONFIG_FILE" '.swarm.topology // "hierarchical"' "hierarchical")
    SWARM_MAX=$(jq_safe "$CONFIG_FILE" '.swarm.maxAgents // 15' "15")
  fi

  # PRIMARY: Count REAL Task tool agents (these are actually running)
  TASK_AGENTS=$(count_task_agents)
  [ "$TASK_AGENTS" -gt 1 ] && SWARM_STATUS="active"

  # SECONDARY: Count MCP registry agents (for capacity display)
  if [ -f "$AGENTS_STORE" ]; then
    REGISTRY_AGENTS=$(jq -r '.agents | length // 0' "$AGENTS_STORE" 2>/dev/null || echo "0")
    # Get types for display
    ACTIVE_AGENT_LIST=$(jq -r '.agents | to_entries | .[].value.agentType // empty' "$AGENTS_STORE" 2>/dev/null | sort | uniq | head -5 | tr '\n' ' ')
  fi

  # Count queued tasks
  if [ -f "$TASKS_FILE" ]; then
    SWARM_TASKS=$(jq -r '[.tasks[] | select(.status == "pending" or .status == "in_progress")] | length // 0' "$TASKS_FILE" 2>/dev/null || echo "0")
  fi

  # Check hive-mind coordination state
  if [ -f "$HIVE_FILE" ]; then
    local hive_active=$(jq_safe "$HIVE_FILE" '.active // false' "false")
    [ "$hive_active" = "true" ] && SWARM_COORD="true" && SWARM_STATUS="coordinated"
  fi

  # For display: use TASK_AGENTS as primary, fall back to REGISTRY for capacity
  # SWARM_AGENTS = currently running (Task tool)
  # SWARM_MAX = capacity (from config or registry count)
  SWARM_AGENTS="$TASK_AGENTS"
  [ "$REGISTRY_AGENTS" -gt "$SWARM_MAX" ] && SWARM_MAX="$REGISTRY_AGENTS"
}

# =============================================================================
# Performance Data
# =============================================================================
load_perf() {
  HNSW_MS="0"; FLASH_SPEEDUP="1.0"; CACHE_HIT="0"; MEM_SAVED="0"

  if [ -f "$PERF_FILE" ]; then
    HNSW_MS=$(jq_safe "$PERF_FILE" '.hnsw.avgMs // .hnsw.searchTimeMs // 0' "0")
    FLASH_SPEEDUP=$(jq_safe "$PERF_FILE" '.flashAttention.speedup // 1.0' "1.0")
    # Strip trailing 'x' if present (some JSONs have "1.0x" format)
    FLASH_SPEEDUP="${FLASH_SPEEDUP%x}"
    CACHE_HIT=$(jq_safe "$PERF_FILE" '.cache.hitRate // 0' "0")
    MEM_SAVED=$(jq_safe "$PERF_FILE" '.memory.savedPercent // 0' "0")
  fi
}

# =============================================================================
# V3 Progress Data (NOTE: These are internal dev metrics, not user-facing)
# =============================================================================
load_progress() {
  V3_PROGRESS=0; V3_DOMAINS=0; V3_TOTAL_DOMAINS=5
  DDD_PROGRESS=0

  if [ -f "$PROGRESS_FILE" ]; then
    # FIX: Try multiple field names for compatibility
    V3_PROGRESS=$(jq_safe "$PROGRESS_FILE" '.overall.percentage // .percentage // .progress // 0' "0")
    V3_DOMAINS=$(jq_safe "$PROGRESS_FILE" '.domains.completed // 0' "0")
    V3_TOTAL_DOMAINS=$(jq_safe "$PROGRESS_FILE" '.domains.total // 5' "5")
    # FIX: JSON has .ddd.progress not .ddd.percentage
    DDD_PROGRESS=$(jq_safe "$PROGRESS_FILE" '.ddd.percentage // .ddd.progress // 0' "0")
  fi
}

# =============================================================================
# Cost/Token Data (from stdin or CLI)
# =============================================================================
load_cost() {
  COST_SESSION="0.00"; COST_TOTAL="0.00"
  TOKENS_IN=0; TOKENS_OUT=0; TOKENS_CACHE=0
  CTX_PCT=0; CTX_USED=0; CTX_MAX=200000
  CACHE_WRITE=0; CACHE_READ=0

  # Try stdin first
  COST_SESSION=$(echo "$INPUT" | jq -r '.cost.session // .cost_usd // "0.00"' 2>/dev/null || echo "0.00")
  TOKENS_IN=$(echo "$INPUT" | jq -r '.tokens.input // .input_tokens // 0' 2>/dev/null || echo "0")
  TOKENS_OUT=$(echo "$INPUT" | jq -r '.tokens.output // .output_tokens // 0' 2>/dev/null || echo "0")
  TOKENS_CACHE=$(echo "$INPUT" | jq -r '.tokens.cache // 0' 2>/dev/null || echo "0")
  CTX_PCT=$(echo "$INPUT" | jq -r '.context.percentage // 0' 2>/dev/null || echo "0")
  CACHE_WRITE=$(echo "$INPUT" | jq -r '.cache.write // 0' 2>/dev/null || echo "0")
  CACHE_READ=$(echo "$INPUT" | jq -r '.cache.read // 0' 2>/dev/null || echo "0")
}

# =============================================================================
# Formatting Helpers
# =============================================================================
fmt_cost() { printf "\$%.4f" "$1"; }
fmt_pct() { printf "%.0f%%" "$1"; }
fmt_num() {
  local n="$1"
  if [ "$n" -ge 1000000 ]; then
    printf "%.1fM" "$(echo "$n / 1000000" | bc -l)"
  elif [ "$n" -ge 1000 ]; then
    printf "%.1fK" "$(echo "$n / 1000" | bc -l)"
  else
    printf "%d" "$n"
  fi
}

# Progress bar (5 chars: ▰▱)
progress_bar() {
  local pct="$1"
  local filled=$((pct / 20))
  local empty=$((5 - filled))
  local bar=""
  for ((i=0; i<filled; i++)); do bar+="▰"; done
  for ((i=0; i<empty; i++)); do bar+="▱"; done
  echo "$bar"
}

# =============================================================================
# Alert Header (for alert mode)
# =============================================================================
print_alert_header() {
  local cves="$TOTAL_CVES"
  local critical="$CRITICAL_CVES"
  printf "${RED}⚠️ ALERT${RST}: "
  if [ "$critical" -gt 0 ]; then
    printf "${RED}${BOLD}$critical critical${RST} "
  fi
  printf "$cves CVEs ${DIM}│${RST} Run: ${CYN}security scan${RST}\n"
}

# =============================================================================
# Swarm Header (for swarm mode)
# =============================================================================
print_swarm_header() {
  printf "${GRN}🐝 SWARM${RST}: ${BOLD}$SWARM_AGENTS${RST}/$SWARM_MAX agents │ $SWARM_TOPOLOGY"
  [ "$SWARM_TASKS" -gt 0 ] && printf " │ $SWARM_TASKS tasks"
  printf "\n"
}

# =============================================================================
# Initialize Everything
# =============================================================================
init_statusline() {
  setup_colors
  load_stdin
  setup_paths
}
