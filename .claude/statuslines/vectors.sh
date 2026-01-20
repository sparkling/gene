#!/bin/bash
# vectors - Vector database statusline (1 line)
# Category: Knowledge Base (Vector/Learning)
# Format: Model │ vectors count per database

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
    load_security
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    printf " ${DIM}│${RST} 📊 $USER_VECTORS user"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${GRN}$SWARM_AGENTS agents${RST}"
    printf " ${DIM}│${RST} 📊 $USER_VECTORS vectors"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} ${BLU}📊${RST} ${BOLD}$(fmt_num $USER_VECTORS)${RST} user"
    printf " ${DIM}│${RST} ${BOLD}$(fmt_num $OPS_VECTORS)${RST} ops"
    printf " ${DIM}│${RST} ${BOLD}$(fmt_num $RV_MEMORIES)${RST} memories"
    ;;
esac
echo
