#!/bin/bash
# Minimal Statusline - Identity only (1 line)
# Format: Model in dir on ⎇ branch

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

# ANSI Colors
RST="\033[0m"; BOLD="\033[1m"
CYN="\033[36m"; YEL="\033[33m"

# Single line output
printf "${BOLD}$MODEL${RST} in ${CYN}$DIR${RST}"
[ -n "$BRANCH" ] && printf " on ${YEL}⎇ $BRANCH${RST}"
echo
