# Risk Overview

**Document ID:** 90-RISK-OVERVIEW
**Status:** Final
**Owner:** Risk Management
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Platform faces 12 key risks across 5 categories: Competitive, Technical, Regulatory, Market, and Operational. Top 3 threats: SelfDecode adding traditional medicine (HIGH), well-funded startup entry (MEDIUM-HIGH), regulatory crackdown (MEDIUM). Mitigation strategy: launch fast to establish network effects, educational positioning to avoid medical device classification, community moat to create switching costs.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Risk framework | 5-category model | Comprehensive coverage of threat landscape | Jan 2026 |
| Top priority risk | Competitive entry | Direct threat to differentiation | Jan 2026 |
| Primary mitigation | Speed to market | Network effects compound early | Jan 2026 |
| Regulatory approach | Educational positioning | Avoids FDA medical device classification | Jan 2026 |

---

## Risk Framework

```
                           RISK CATEGORIES
                                 │
    ┌────────────┬───────────────┼───────────────┬────────────┐
    │            │               │               │            │
COMPETITIVE  TECHNICAL      REGULATORY      MARKET      OPERATIONAL
    │            │               │               │            │
┌───┴───┐    ┌───┴───┐       ┌───┴───┐      ┌───┴───┐    ┌───┴───┐
│ Direct │    │ Data  │       │ FDA   │      │ Demand│    │ Team  │
│ Entry  │    │Quality│       │ Device│      │ Valid │    │ Scale │
├────────┤    ├───────┤       ├───────┤      ├───────┤    ├───────┤
│ Funded │    │ AI    │       │ FTC   │      │ Market│    │ Data  │
│Startup │    │ Costs │       │ Claims│      │ Timing│    │ Ops   │
├────────┤    ├───────┤       ├───────┤      ├───────┤    ├───────┤
│Feature │    │ DNA   │       │ HIPAA │      │Chicken│    │Support│
│ Copy   │    │Parsing│       │ GDPR  │      │ Egg   │    │ Scale │
└────────┘    └───────┘       └───────┘      └───────┘    └───────┘
```

---

## Top 12 Risks (Ranked)

| Rank | Risk | Category | Likelihood | Impact | Score | Trend |
|------|------|----------|------------|--------|-------|-------|
| 1 | SelfDecode adds TCM/Ayurveda | Competitive | HIGH | HIGH | **9** | ↑ |
| 2 | Well-funded startup enters space | Competitive | MEDIUM | HIGH | **7** | → |
| 3 | FDA medical device classification | Regulatory | MEDIUM | HIGH | **7** | → |
| 4 | AI/LLM costs exceed projections | Technical | MEDIUM | MEDIUM | **5** | ↓ |
| 5 | TCM/Ayurveda data quality issues | Technical | MEDIUM | MEDIUM | **5** | → |
| 6 | Practitioner marketplace chicken-egg | Market | HIGH | MEDIUM | **5** | → |
| 7 | FTC health claims enforcement | Regulatory | LOW | HIGH | **4** | → |
| 8 | DNA parsing complexity delays MVP | Technical | MEDIUM | MEDIUM | **4** | ↓ |
| 9 | User demand lower than projected | Market | LOW | HIGH | **4** | → |
| 10 | Key hire failure | Operational | MEDIUM | MEDIUM | **4** | → |
| 11 | Data breach / security incident | Operational | LOW | HIGH | **4** | → |
| 12 | GDPR compliance gaps | Regulatory | LOW | MEDIUM | **3** | → |

**Scoring:** Likelihood × Impact (1-3 scale each, max 9)

---

## Risk Deep-Dives

### Risk #1: SelfDecode Adds Traditional Medicine

| Factor | Assessment |
|--------|------------|
| **Threat** | SelfDecode (200M+ SNPs, DecodyGPT) adds TCM/Ayurveda integration |
| **Likelihood** | HIGH — They have resources, market signal |
| **Impact** | HIGH — Direct competition to our core differentiator |
| **Lead Time** | 6-12 months for basic implementation |
| **Mitigation** | Launch fast, build community moat, deeper integration |

**Warning Signs:**
- SelfDecode announces traditional medicine partnerships
- Job postings for TCM/Ayurveda expertise
- User forum discussions about traditional medicine demand

### Risk #2: Well-Funded Startup Entry

| Factor | Assessment |
|--------|------------|
| **Threat** | Startup with $10M+ funding targets genetics + traditional medicine |
| **Likelihood** | MEDIUM-HIGH — Market opportunity visible |
| **Impact** | HIGH — Could outspend on marketing, talent |
| **Lead Time** | 12-18 months to build competitive product |
| **Mitigation** | Focus on practitioner relationships, community |

**Warning Signs:**
- Y Combinator/a16z investments in adjacent space
- Acquisitions in genetic testing or TCM
- New competitor product announcements

### Risk #3: FDA Medical Device Classification

| Factor | Assessment |
|--------|------------|
| **Threat** | FDA classifies platform as Software as Medical Device (SaMD) |
| **Likelihood** | MEDIUM — Regulatory scrutiny increasing |
| **Impact** | HIGH — Would require 510(k) clearance, massive delays |
| **Lead Time** | Enforcement typically follows warnings |
| **Mitigation** | Educational positioning, clear disclaimers, legal review |

**Warning Signs:**
- FDA warning letters to similar platforms
- Congressional hearings on DTC genetic testing
- FTC enforcement actions in health space

---

## Risk Category Summaries

### Competitive Risks

| Risk | Level | Mitigation |
|------|-------|------------|
| SelfDecode adds TCM | HIGH | Launch fast, community moat |
| Well-funded startup | MEDIUM-HIGH | Practitioner relationships |
| Feature copying | MEDIUM | Network effects, unique data |

### Technical Risks

| Risk | Level | Mitigation |
|------|-------|------------|
| AI cost overruns | MEDIUM | Caching, query optimization |
| Data quality issues | MEDIUM | Validate top 100 herbs first |
| DNA parsing complexity | MEDIUM | Start with 23andMe only |

### Regulatory Risks

| Risk | Level | Mitigation |
|------|-------|------------|
| FDA SaMD classification | MEDIUM | Educational positioning |
| FTC health claims | LOW-MEDIUM | Conservative claims, disclaimers |
| HIPAA compliance | LOW | Standard security practices |
| GDPR compliance | LOW | Privacy-first architecture |

### Market Risks

| Risk | Level | Mitigation |
|------|-------|------------|
| Marketplace chicken-egg | HIGH | Seed with initial practitioners |
| Lower demand than projected | LOW-MEDIUM | Validate with beta users |
| Market timing | LOW | Growing awareness of genetics |

### Operational Risks

| Risk | Level | Mitigation |
|------|-------|------------|
| Key hire failure | MEDIUM | Competitive compensation, mission |
| Data breach | LOW-MEDIUM | Security-first development |
| Support scaling | LOW | Community support, documentation |

---

## Risk Monitoring Dashboard

### Key Risk Indicators (KRIs)

| Indicator | Threshold | Frequency | Owner |
|-----------|-----------|-----------|-------|
| Competitor product launches | Any announcement | Weekly | Product |
| FDA/FTC enforcement actions | Any in category | Monthly | Legal |
| AI cost per query | >$0.05 | Daily | Engineering |
| Data quality errors | >5% error rate | Weekly | Data |
| User complaints (regulatory) | Any | Daily | Support |
| Security incidents | Any | Real-time | Security |

### Risk Review Cadence

| Review Type | Frequency | Participants |
|-------------|-----------|--------------|
| Tactical risk check | Weekly | Product, Engineering |
| Strategic risk review | Monthly | Leadership |
| Regulatory risk assessment | Quarterly | Legal, Compliance |
| Full risk audit | Annual | Board, External |

---

## Mitigation Strategy Summary

### Speed to Market (Primary Strategy)

```
LAUNCH FAST → BUILD COMMUNITY → NETWORK EFFECTS → DEFENSIBLE MOAT
     │              │                  │                │
     │              │                  │                └── Switching costs
     │              │                  └── More users → more value
     │              └── Practitioners, peer groups
     └── Beat competitors to market
```

### Key Mitigation Investments

| Investment | Purpose | Budget Allocation |
|------------|---------|-------------------|
| Fast MVP development | Beat competitors | 60% of Phase 1 |
| Legal/compliance review | Regulatory protection | 5% ongoing |
| Security infrastructure | Breach prevention | 10% of tech |
| Community building | Network moat | 15% of marketing |

---

## Contingency Plans

### If SelfDecode Adds Traditional Medicine

1. **Immediate**: Accelerate TCM/Ayurveda depth
2. **Short-term**: Emphasize community differentiator
3. **Medium-term**: Expand to Kampo, Western herbal
4. **Long-term**: Become practitioner ecosystem

### If FDA Enforcement Begins

1. **Immediate**: Pause marketing claims
2. **Short-term**: Legal review of all content
3. **Medium-term**: Reposition as educational tool
4. **Long-term**: Consider clinical pathway if needed

### If AI Costs Exceed Budget

1. **Immediate**: Implement aggressive caching
2. **Short-term**: Query optimization, batching
3. **Medium-term**: Fine-tuned smaller models
4. **Long-term**: Open-source model deployment

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [91-ASSUMPTIONS](./91-ASSUMPTIONS.md) | Risk depends on assumptions |
| [92-THREATS](./92-THREATS.md) | Detailed threat analysis |
| [93-MITIGATIONS](./93-MITIGATIONS.md) | Full mitigation playbooks |
| [23-COMPETITORS](../20-market/23-COMPETITORS.md) | Competitive threats |

---

## Open Questions

- [ ] Insurance requirements for health-related platform
- [ ] International regulatory variations (EU, Asia)
- [ ] Cybersecurity insurance adequacy
- [ ] Practitioner liability considerations

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Risk Management | Complete risk overview |
