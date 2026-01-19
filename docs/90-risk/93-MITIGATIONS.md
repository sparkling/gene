# Risk Mitigations

**Document ID:** 93-MITIGATIONS
**Status:** Final
**Owner:** Strategy
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Mitigation strategies address top threats: (1) Competitive entry — build community moat, launch fast, deepen traditional medicine expertise; (2) Data breach — defense in depth, encryption, SOC 2 certification; (3) Regulatory — educational positioning, conservative claims, legal review process. Key principle: moats over features. Contingency playbooks defined for crisis scenarios.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary defense | Community moat | Network effects hardest to replicate | Jan 2026 |
| Security approach | Defense in depth | Genetic data requires highest protection | Jan 2026 |
| Regulatory posture | Conservative | Avoid triggering enforcement | Jan 2026 |
| Contingency trigger | Early indicators | Proactive vs reactive | Jan 2026 |

---

## Mitigation Framework

### Mitigation Hierarchy

```
MITIGATION PRIORITY
        │
        ├── PREVENT: Stop threats from materializing
        │   └── Build moats, security controls, compliance
        │
        ├── DETECT: Identify threats early
        │   └── Monitoring, early warning indicators
        │
        ├── RESPOND: Act quickly when threats materialize
        │   └── Playbooks, decision frameworks
        │
        └── RECOVER: Return to normal operations
            └── Backup plans, insurance, reserves
```

---

## Competitive Threat Mitigations

### Strategy: Build Defensible Moats

| Moat Type | Description | Investment | Timeline |
|-----------|-------------|------------|----------|
| **Community** | Network effects from user-to-user connections | HIGH | Ongoing |
| **Data depth** | Traditional medicine knowledge base | HIGH | Y1-Y2 |
| **Practitioner** | Relationship switching costs | MEDIUM | Y2-Y3 |
| **Brand** | "THREE WORLDS" positioning ownership | MEDIUM | Ongoing |
| **Content** | SEO compound effects | MEDIUM | Ongoing |

### C1 Mitigation: SelfDecode Adds Traditional Medicine

**Prevention:**
- Launch TCM/Ayurveda before they move
- Establish thought leadership in space
- Partner with traditional medicine experts
- Build community that SelfDecode can't copy

**Detection:**
- Weekly monitoring of SelfDecode job postings
- Google Alerts for SelfDecode + traditional medicine
- Track their feature announcements
- Monitor their partnership announcements

**Response Playbook:**

```
IF SelfDecode announces traditional medicine features THEN:
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ IMMEDIATE (Within 24 hours)                           │
│ □ Internal assessment of their announcement           │
│ □ Draft "depth vs breadth" positioning                │
│ □ Prepare community messaging                         │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ SHORT TERM (Within 1 week)                            │
│ □ Blog post: "Why depth matters in traditional med"   │
│ □ Community reassurance messaging                     │
│ □ Accelerate unique features roadmap                  │
│ □ Practitioner outreach emphasizing our expertise     │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ MEDIUM TERM (Within 1 month)                          │
│ □ Feature comparison content (objective, fact-based)  │
│ □ Double down on community features                   │
│ □ Expert testimonials campaign                        │
│ □ Consider pricing promotion for lock-in              │
└───────────────────────────────────────────────────────┘
```

### C2 Mitigation: Well-Funded Startup Enters

**Prevention:**
- Move fast to establish market position
- Build community that's hard to acquire
- Consider strategic partnerships that deter entry
- Focus on niche (chronic illness) less attractive to VCs

**Detection:**
- Monitor Y Combinator, Techstars cohorts
- Track funding announcements in adjacent spaces
- Watch for talent movement from competitors
- Industry conference presentations

**Response Playbook:**

```
IF well-funded competitor announces THEN:
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ ASSESSMENT (Within 48 hours)                          │
│ □ Analyze their positioning, funding, team            │
│ □ Identify overlap with our target market             │
│ □ Assess timeline to product launch                   │
│ □ Determine if threat is real or hype                 │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ IF SIGNIFICANT THREAT:                                │
│ □ Accelerate community building                       │
│ □ Deepen practitioner relationships                   │
│ □ Consider strategic fundraise                        │
│ □ Explore partnership or acquisition conversations    │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ IF LIMITED THREAT:                                    │
│ □ Monitor quarterly                                   │
│ □ Continue execution                                  │
│ □ Learn from their positioning/messaging              │
└───────────────────────────────────────────────────────┘
```

---

## Security Threat Mitigations

### O1 Mitigation: Data Breach

**Prevention — Defense in Depth:**

| Layer | Control | Status |
|-------|---------|--------|
| Edge | Cloudflare WAF, DDoS protection | Planned |
| Application | JWT auth, RBAC, input validation | Planned |
| Data | AES-256 encryption at rest | Planned |
| Data | Field-level encryption for genetic data | Planned |
| Transit | TLS 1.3 for all connections | Implemented |
| Infrastructure | VPC isolation, security groups | Planned |
| Process | Access reviews, audit logging | Planned |

**Detection:**
- Real-time security monitoring (SIEM)
- Anomaly detection on data access patterns
- Regular penetration testing
- Bug bounty program (future)

**Response Playbook:**

```
IF data breach detected THEN:
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ IMMEDIATE (Within 1 hour)                             │
│ □ Activate incident response team                     │
│ □ Contain breach (isolate systems if needed)          │
│ □ Preserve evidence                                   │
│ □ Initial scope assessment                            │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ SHORT TERM (Within 24 hours)                          │
│ □ Full scope assessment                               │
│ □ Legal counsel engagement                            │
│ □ Regulatory notification preparation                 │
│ □ User notification if required (GDPR: 72 hours)      │
│ □ PR statement preparation                            │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ MEDIUM TERM (Within 1 week)                           │
│ □ Complete incident investigation                     │
│ □ Remediation of vulnerabilities                      │
│ □ User communication (transparency)                   │
│ □ Regulatory filings as required                      │
│ □ Post-mortem and lessons learned                     │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ RECOVERY (Within 1 month)                             │
│ □ Third-party security audit                          │
│ □ Updated security controls                           │
│ □ Trust rebuilding campaign                           │
│ □ Process improvements                                │
└───────────────────────────────────────────────────────┘
```

**Recovery:**
- Cyber insurance coverage
- Legal retainer for breach response
- PR firm on standby for crisis communication
- Reserve fund for remediation

---

## Regulatory Threat Mitigations

### R1 Mitigation: FDA Reclassification

**Prevention:**
- Maintain educational positioning
- Never make diagnostic or treatment claims
- Legal review of all health content
- Document compliance posture

**Content Guidelines:**

| DO | DON'T |
|----|-------|
| "Research suggests..." | "This causes your symptoms" |
| "Studies have associated..." | "You have this condition" |
| "You may want to discuss..." | "Take this to treat..." |
| "Some find benefit from..." | "This will cure..." |

**Detection:**
- Monitor FDA warning letters (weekly)
- Track competitor enforcement actions
- Industry association communications
- Congressional hearing topics

**Response Playbook:**

```
IF FDA issues warning or guidance affecting us THEN:
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ IMMEDIATE (Within 24 hours)                           │
│ □ Legal counsel review of warning/guidance            │
│ □ Assess applicability to our platform                │
│ □ Internal communications to team                     │
│ □ Pause marketing pending review                      │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ SHORT TERM (Within 1 week)                            │
│ □ Comprehensive content audit                         │
│ □ Identify and remediate non-compliant content        │
│ □ Update disclaimers if needed                        │
│ □ User communication if material changes              │
└───────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ IF SEVERE (Classification as medical device):         │
│ □ Assess compliance pathway vs pivot                  │
│ □ Consider practitioner-only model                    │
│ □ Explore partnerships with compliant entities        │
│ □ Prepare for potential service suspension            │
└───────────────────────────────────────────────────────┘
```

### R2 Mitigation: FTC Enforcement

**Prevention:**
- Conservative health claims only
- All claims evidence-linked
- Clear disclosure of limitations
- Legal review process for content

**Claim Review Process:**

```
CONTENT APPROVAL WORKFLOW
        │
        ▼
┌───────────────┐
│ Draft content │
│ (writer)      │
└───────┬───────┘
        │
        ▼
┌───────────────┐
│ Evidence check│
│ (researcher)  │ ──► Must cite peer-reviewed source
└───────┬───────┘
        │
        ▼
┌───────────────┐
│ Compliance    │
│ flag check    │ ──► Automated scan for risky terms
└───────┬───────┘
        │
    ┌───┴───┐
    │       │
    ▼       ▼
┌───────┐ ┌───────┐
│ Clear │ │Flagged│
│(publish)│(legal)│
└───────┘ └───┬───┘
              │
              ▼
        ┌───────────┐
        │ Legal     │
        │ review    │
        └─────┬─────┘
              │
        ┌─────┴─────┐
        │           │
        ▼           ▼
    ┌───────┐  ┌───────┐
    │Approve│  │Reject/│
    │       │  │Revise │
    └───────┘  └───────┘
```

**Risky Terms to Flag:**
- "Cure," "treat," "prevent," "diagnose"
- "Clinically proven," "guaranteed"
- Specific disease names without disclaimer
- "Doctor recommended" without substantiation

---

## Operational Threat Mitigations

### O2 Mitigation: Key Person Risk

**Prevention:**
- Document all critical processes
- Cross-train team members
- Avoid single points of failure
- Key person insurance consideration

**Critical Knowledge Areas:**

| Area | Current Owner | Backup | Status |
|------|---------------|--------|--------|
| Architecture | Founder | CTO hire | Documented |
| Traditional medicine | Founder | Advisory | Partially documented |
| Financial | Founder | Accountant | Documented |
| Legal/Compliance | Legal counsel | Secondary counsel | Documented |
| Community relationships | Founder | Community manager | Needs documentation |

**Documentation Requirements:**
- Architecture Decision Records (ADRs)
- Process playbooks
- Contact and relationship lists
- Access credentials (secure vault)
- Strategic decision history

### O3 Mitigation: Supplier Dependency

**Critical Suppliers:**

| Supplier | Risk | Mitigation |
|----------|------|------------|
| AWS/Cloud | Platform failure | Multi-region, backup provider plan |
| Anthropic/Claude | API changes, pricing | Abstraction layer, fallback models |
| Stripe | Payment processing | Backup processor setup |
| NCBI/PubMed | Data access | Local caching, archival |
| Domain registrar | Domain loss | Premium protection, backup domains |

**Architecture Principles:**
- Abstraction layers for all external APIs
- Multi-provider capability where feasible
- Graceful degradation design
- Data locality and caching

---

## Early Warning System

### Indicator Dashboard

| Category | Indicator | Frequency | Owner |
|----------|-----------|-----------|-------|
| Competitive | SelfDecode job postings | Weekly | Marketing |
| Competitive | Funding announcements | Weekly | CEO |
| Regulatory | FDA warning letters | Weekly | Legal |
| Regulatory | FTC actions | Monthly | Legal |
| Security | Penetration test results | Quarterly | Engineering |
| Security | Anomaly alerts | Real-time | Engineering |
| Market | Customer churn spike | Daily | Customer Success |
| Market | NPS decline | Monthly | Product |

### Escalation Thresholds

| Metric | Yellow | Red | Action |
|--------|--------|-----|--------|
| Daily signups | <50% target | <25% target | Marketing review |
| Error rate | >2% | >5% | Engineering escalation |
| Churn | >5% monthly | >10% monthly | Customer success sprint |
| Security alerts | Any medium | Any high | Immediate investigation |

---

## Contingency Plans

### Scenario: Major Competitor Launches

**If competitor captures >10% of our target market:**
1. Assess differentiation gap
2. Accelerate unique features
3. Double community investment
4. Consider strategic options (partnership, acquisition, pivot)

### Scenario: Regulatory Restriction

**If FDA restricts direct-to-consumer genetics:**
1. Pivot to practitioner-only model
2. Rebrand as educational platform
3. Partner with compliant healthcare entities
4. Explore international markets

### Scenario: Funding Shortfall

**If unable to raise planned funding:**
1. Reduce burn rate (team, marketing)
2. Focus on revenue-generating features
3. Explore revenue-based financing
4. Consider strategic acquisition

### Scenario: Founder Incapacity

**If founder unavailable for extended period:**
1. Advisory board assumes strategic guidance
2. CTO/COO handles operations
3. Legal counsel handles external matters
4. Board identifies interim leadership

---

## Risk Reserve

### Financial Reserve

| Purpose | Amount | Trigger |
|---------|--------|---------|
| Operating runway | 6 months | Always maintain |
| Legal defense | $50,000 | Regulatory action |
| Security incident | $25,000 | Breach response |
| Competitive response | $25,000 | Emergency marketing |

### Relationship Reserve

| Relationship | Purpose | Status |
|--------------|---------|--------|
| PR firm | Crisis communication | On retainer |
| Legal counsel | Regulatory, IP, breach | On retainer |
| Security firm | Incident response | Identified |
| M&A advisor | Strategic options | Identified |

---

## Review and Update

### Quarterly Risk Review

| Activity | Owner | Timing |
|----------|-------|--------|
| Threat assessment update | Strategy | Q1, Q2, Q3, Q4 |
| Mitigation effectiveness review | Operations | Q1, Q2, Q3, Q4 |
| Playbook testing | All | Q2, Q4 |
| Early warning system audit | Analytics | Q1, Q3 |

### Annual Activities

| Activity | Owner | Timing |
|----------|-------|--------|
| Full risk assessment | Strategy | Q1 |
| Security audit | External firm | Q2 |
| Insurance review | Finance | Q3 |
| Contingency plan update | Operations | Q4 |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [90-RISK-OVERVIEW](./90-RISK-OVERVIEW.md) | Risk framework |
| [92-THREATS](./92-THREATS.md) | Threat details |
| [74-COMPLIANCE](../70-operations/74-COMPLIANCE.md) | Compliance mitigations |
| [44-ARCHITECTURE](../40-product/44-ARCHITECTURE.md) | Security architecture |

---

## Open Questions

- [ ] Determine cyber insurance coverage level
- [ ] Establish PR firm crisis retainer
- [ ] Complete business continuity documentation
- [ ] Test incident response playbook

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Strategy | Complete mitigation framework |
