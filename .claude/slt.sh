#!/bin/bash
# Fast statusline switcher - handles all /slt operations
# Usage: slt.sh [list|next|n|prev|p|<template>]

STATE=".claude/statusline-state"
DIR=".claude/statuslines"
CURRENT=$(cat "$STATE" 2>/dev/null | tr -d '[:space:]')
[ -z "$CURRENT" ] && CURRENT="adaptive"

# Get sorted template list
TEMPLATES=($(ls "$DIR"/*.sh 2>/dev/null | xargs -n1 basename | sed 's/\.sh$//' | sort))
COUNT=${#TEMPLATES[@]}

# Find current index
for i in $(seq 0 $((COUNT-1))); do
  [ "${TEMPLATES[$i]}" = "$CURRENT" ] && IDX=$i && break
done

# Template descriptions with scenarios
get_info() {
  case "$1" in
    zen)       echo "1|Model only|Deep focus, complex debugging" ;;
    focus)     echo "1|Model + branch|Focus work with git context" ;;
    dev)       echo "1|Dir, branch, diff stats|Daily development" ;;
    git)       echo "2|Full git details|Code reviews, preparing commits" ;;
    vectors)   echo "1|Vector DB counts|Building knowledge base" ;;
    learning)  echo "2|Patterns, trajectories|Training, RAG development" ;;
    daemon)    echo "1|Daemon status|System monitoring" ;;
    workers)   echo "3|Per-worker stats|Debugging workers" ;;
    cost)      echo "1|Cost, tokens|Budget tracking" ;;
    tokens)    echo "2|Token breakdown|Cache optimization" ;;
    secure)    echo "1|Security status|Security awareness" ;;
    audit)     echo "2|CVE breakdown|Security reviews" ;;
    swarm)     echo "1|Agent count, topology|Multi-agent work" ;;
    agents)    echo "3|Active agent types|Swarm debugging" ;;
    perf)      echo "1|HNSW, Flash metrics|Performance tuning" ;;
    speed)     echo "2|Cache, worker timing|Bottleneck analysis" ;;
    progress)  echo "1|Domains, DDD %|Project tracking" ;;
    project)   echo "2|V3 progress bar|Sprint reviews" ;;
    dashboard) echo "3|All key metrics|Full overview" ;;
    full)      echo "6|Everything visible|Comprehensive monitoring" ;;
    adaptive)  echo "2-4|Auto-detects mode|General purpose" ;;
    compact)   echo "1|Key metrics|Quick overview" ;;
    minimal)   echo "1|Model, dir, branch|Simple identity" ;;
    *)         echo "?|Unknown|-" ;;
  esac
}

case "${1:-list}" in
  list|"")
    echo "Statusline Templates (current: $CURRENT)"
    echo ""
    printf "  %-12s %s  %-22s  %s\n" "TEMPLATE" "L" "SHOWS" "BEST FOR"
    printf "  %-12s %s  %-22s  %s\n" "--------" "-" "-----" "--------"
    for t in "${TEMPLATES[@]}"; do
      info=$(get_info "$t")
      lines=$(echo "$info" | cut -d'|' -f1)
      shows=$(echo "$info" | cut -d'|' -f2)
      scenario=$(echo "$info" | cut -d'|' -f3)
      if [ "$t" = "$CURRENT" ]; then
        printf "* %-12s %s  %-22s  %s\n" "$t" "$lines" "$shows" "$scenario"
      else
        printf "  %-12s %s  %-22s  %s\n" "$t" "$lines" "$shows" "$scenario"
      fi
    done
    echo ""
    echo "Switch: .claude/slt.sh <name> | n (next) | p (prev)"
    ;;
  next|n)
    NEW="${TEMPLATES[$(( (IDX + 1) % COUNT ))]}"
    echo "$NEW" > "$STATE"
    info=$(get_info "$NEW")
    shows=$(echo "$info" | cut -d'|' -f2)
    echo "$CURRENT -> $NEW ($shows)"
    ;;
  prev|p)
    NEW="${TEMPLATES[$(( (IDX - 1 + COUNT) % COUNT ))]}"
    echo "$NEW" > "$STATE"
    info=$(get_info "$NEW")
    shows=$(echo "$info" | cut -d'|' -f2)
    echo "$CURRENT -> $NEW ($shows)"
    ;;
  *)
    if [ -f "$DIR/$1.sh" ]; then
      echo "$1" > "$STATE"
      info=$(get_info "$1")
      shows=$(echo "$info" | cut -d'|' -f2)
      echo "$CURRENT -> $1 ($shows)"
    else
      echo "Unknown: $1"
      echo "Available: ${TEMPLATES[*]}"
      exit 1
    fi
    ;;
esac
