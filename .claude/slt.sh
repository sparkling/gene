#!/bin/bash
# Fast statusline switcher - handles all /slt operations
# Usage: slt.sh [list|next|n|prev|p|<template>]

STATE=".claude/statusline-state"
DIR=".claude/statuslines"
CURRENT=$(cat "$STATE" 2>/dev/null | tr -d '[:space:]')
[ -z "$CURRENT" ] && CURRENT="full"

# Get sorted template list (excludes lib directory)
TEMPLATES=($(ls "$DIR"/*.sh 2>/dev/null | xargs -n1 basename | sed 's/\.sh$//' | sort))
COUNT=${#TEMPLATES[@]}

# Find current index
for i in $(seq 0 $((COUNT-1))); do
  [ "${TEMPLATES[$i]}" = "$CURRENT" ] && IDX=$i && break
done

# Template descriptions with scenarios (8 curated templates)
get_info() {
  case "$1" in
    zen)       echo "1|Model only|Deep focus, complex debugging" ;;
    dev)       echo "1|Dir, branch, diff stats|Daily development" ;;
    cost)      echo "1|Cost, tokens, context|Budget tracking" ;;
    daemon)    echo "1|Daemon, workers|System monitoring" ;;
    secure)    echo "1|Security status, CVEs|Security awareness" ;;
    swarm)     echo "1|Running agents, topology|Multi-agent coordination" ;;
    full)      echo "3|All key metrics|Comprehensive monitoring" ;;
    adaptive)  echo "2-3|Auto-detects mode|General purpose (recommended)" ;;
    *)         echo "?|Unknown template|-" ;;
  esac
}

case "${1:-list}" in
  list|"")
    echo "Statusline Templates (current: $CURRENT)"
    echo ""
    printf "  %-12s %s  %-24s  %s\n" "TEMPLATE" "L" "SHOWS" "BEST FOR"
    printf "  %-12s %s  %-24s  %s\n" "--------" "-" "-----" "--------"
    for t in "${TEMPLATES[@]}"; do
      info=$(get_info "$t")
      lines=$(echo "$info" | cut -d'|' -f1)
      shows=$(echo "$info" | cut -d'|' -f2)
      scenario=$(echo "$info" | cut -d'|' -f3)
      if [ "$t" = "$CURRENT" ]; then
        printf "* %-12s %s  %-24s  %s\n" "$t" "$lines" "$shows" "$scenario"
      else
        printf "  %-12s %s  %-24s  %s\n" "$t" "$lines" "$shows" "$scenario"
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
