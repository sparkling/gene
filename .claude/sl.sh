#!/bin/bash
# Fast statusline switcher - handles all /sl operations
# Usage: sl.sh [list|next|n|prev|p|<template>]

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

case "${1:-list}" in
  list|"")
    echo "Statusline templates:"
    for t in "${TEMPLATES[@]}"; do
      [ "$t" = "$CURRENT" ] && echo "* $t" || echo "  $t"
    done
    ;;
  next|n)
    NEW="${TEMPLATES[$(( (IDX + 1) % COUNT ))]}"
    echo "$NEW" > "$STATE"
    echo "$CURRENT -> $NEW"
    ;;
  prev|p)
    NEW="${TEMPLATES[$(( (IDX - 1 + COUNT) % COUNT ))]}"
    echo "$NEW" > "$STATE"
    echo "$CURRENT -> $NEW"
    ;;
  *)
    if [ -f "$DIR/$1.sh" ]; then
      echo "$1" > "$STATE"
      echo "$CURRENT -> $1"
    else
      echo "Unknown: $1"
      echo "Available: ${TEMPLATES[*]}"
      exit 1
    fi
    ;;
esac
