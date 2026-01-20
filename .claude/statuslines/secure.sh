#!/bin/bash
# secure - Security status statusline (1 line)
# Category: Security
# Format: Model │ status │ CVE count │ last scan

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_security

MODE=$(detect_mode)

case "$MODE" in
  alert)
    printf "${RED}⚠ ALERT${RST} ${DIM}│${RST} "
    case "$SECURITY_STATUS" in
      CRITICAL) printf "${RED}🔴 CRITICAL${RST}" ;;
      HIGH)     printf "${RED}HIGH${RST}" ;;
      WARN)     printf "${YEL}WARN${RST}" ;;
      *)        printf "$SECURITY_STATUS" ;;
    esac
    printf " ${DIM}│${RST} ${RED}$TOTAL_CVES CVEs${RST}"
    [ "$CRITICAL_CVES" -gt 0 ] && printf " (${RED}${BOLD}$CRITICAL_CVES crit${RST})"
    printf " ${DIM}│${RST} Run: ${CYN}security scan${RST}"
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡${RST} ${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} ${GRN}🔒 $SECURITY_STATUS${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST}"
    printf " ${DIM}│${RST} "
    case "$SECURITY_STATUS" in
      CLEAN)    printf "${GRN}🔒 CLEAN${RST}" ;;
      CRITICAL) printf "${RED}🔴 CRITICAL${RST}" ;;
      HIGH)     printf "${RED}HIGH${RST}" ;;
      WARN)     printf "${YEL}WARN${RST}" ;;
      unknown)  printf "${DIM}🔒 unknown${RST}" ;;
      *)        printf "$SECURITY_STATUS" ;;
    esac
    printf " ${DIM}│${RST} $TOTAL_CVES CVEs"
    [ "$LAST_SCAN" != "never" ] && printf " ${DIM}│${RST} scan: $LAST_SCAN"
    ;;
esac
echo
