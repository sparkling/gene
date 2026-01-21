#!/bin/bash
# vectors - Memory/database statusline (1 line)
# Category: Knowledge Base
# Shows: Model │ Database sizes │ Patterns │ Memories

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_vectors
load_intel
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    [ "$DB_EXISTS" = "true" ] && printf " ${DIM}│${RST} 💾 $TOTAL_DB_SIZE"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${GRN}$REGISTRY_AGENTS agents${RST}"
    [ "$DB_EXISTS" = "true" ] && printf " ${DIM}│${RST} 💾 $TOTAL_DB_SIZE"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    if [ "$DB_EXISTS" = "true" ]; then
      printf " ${DIM}│${RST} ${BLU}💾${RST} user:$USER_SIZE"
      printf " ops:$OPS_SIZE"
      printf " mem:$MEMORY_SIZE"
    else
      printf " ${DIM}│ No databases${RST}"
    fi
    [ "$RV_PATTERNS" -gt 0 ] && printf " ${DIM}│${RST} ${MAG}◆${RST}$RV_PATTERNS patterns"
    [ "$RV_MEMORIES" -gt 0 ] && printf " ${DIM}│${RST} ${YEL}●${RST}$RV_MEMORIES memories"
    ;;
esac
echo
