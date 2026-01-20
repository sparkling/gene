#!/bin/bash
# audit - Security audit details statusline (2 lines)
# Category: Security
# Line 1: Status
# Line 2: Issue breakdown

source "$(dirname "$0")/lib/common.sh"
init_statusline

MODEL_SHORT=$(shorten_model)
setup_paths
load_security

MODE=$(detect_mode)

# Line 1
case "$MODE" in
  alert)
    printf "${RED}⚠ SECURITY ALERT${RST}"
    case "$SECURITY_STATUS" in
      CRITICAL) printf " ${DIM}│${RST} ${RED}🔴 CRITICAL${RST}" ;;
      HIGH)     printf " ${DIM}│${RST} ${RED}HIGH RISK${RST}" ;;
      *)        printf " ${DIM}│${RST} ${YEL}ISSUES FOUND${RST}" ;;
    esac
    ;;
  swarm)
    load_swarm
    printf "${GRN}⬡ SWARM${RST}: $SWARM_AGENTS agents"
    printf " ${DIM}│${RST} Security: ${GRN}$SECURITY_STATUS${RST}"
    ;;
  *)
    printf "${BOLD}$MODEL_SHORT${RST} ${DIM}│${RST} ${CYN}🔒${RST} Security Audit"
    case "$SECURITY_STATUS" in
      CLEAN)    printf " ${DIM}│${RST} ${GRN}CLEAN${RST}" ;;
      CRITICAL) printf " ${DIM}│${RST} ${RED}CRITICAL${RST}" ;;
      HIGH)     printf " ${DIM}│${RST} ${RED}HIGH${RST}" ;;
      WARN)     printf " ${DIM}│${RST} ${YEL}WARN${RST}" ;;
      *)        printf " ${DIM}│${RST} ${DIM}pending${RST}" ;;
    esac
    ;;
esac
echo

# Line 2: Issue breakdown
if [ "$TOTAL_CVES" -gt 0 ]; then
  printf "${DIM}CVEs:${RST}"
  [ "$CRITICAL_CVES" -gt 0 ] && printf " ${RED}$CRITICAL_CVES critical${RST}"
  [ "$HIGH_CVES" -gt 0 ] && printf " ${YEL}$HIGH_CVES high${RST}"
  [ "$MEDIUM_CVES" -gt 0 ] && printf " ${DIM}$MEDIUM_CVES medium${RST}"

  LOW=$((TOTAL_CVES - CRITICAL_CVES - HIGH_CVES - MEDIUM_CVES))
  [ "$LOW" -gt 0 ] && printf " ${DIM}$LOW low${RST}"

  printf " ${DIM}│${RST} Total: ${BOLD}$TOTAL_CVES${RST}"
elif [ "$LAST_SCAN" != "never" ]; then
  printf "${GRN}✓${RST} No vulnerabilities found ${DIM}│${RST} Last scan: $LAST_SCAN"
else
  printf "${DIM}No security scan data ${DIM}│${RST} Run: ${CYN}npx @claude-flow/cli@latest security scan${RST}"
fi
echo
