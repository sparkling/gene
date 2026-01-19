# Threat Analysis

**Document ID:** 92-THREATS
**Status:** Final
**Owner:** Strategy
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Primary threats: (1) SelfDecode adds traditional medicine — MEDIUM-HIGH, could erode differentiation; (2) Well-funded startup enters space — MEDIUM-HIGH, resource disadvantage; (3) Regulatory crackdown — MEDIUM, could limit health claims; (4) Data breach — MEDIUM, genetic data is sensitive. Defensive moats: community network effects, traditional medicine data depth, practitioner relationships. Early warning indicators established for all major threats.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary threat focus | Competitive entry | Most likely, most impactful | Jan 2026 |
| Risk tolerance | Moderate | Startup stage, some risk acceptable | Jan 2026 |
| Monitoring approach | Quarterly review | Balance vigilance and overhead | Jan 2026 |

---

## Threat Framework

### Threat Categories

```
THREAT LANDSCAPE
        │
        ├── COMPETITIVE
        │   ├── Existing players expand
        │   ├── Well-funded new entrants
        │   └── Big tech entry
        │
        ├── MARKET
        │   ├── Demand decline
        │   ├── Market consolidation
        │   └── Pricing pressure
        │
        ├── TECHNOLOGY
        │   ├── Platform obsolescence
        │   ├── AI commoditization
        │   └── Data source changes
        │
        ├── REGULATORY
        │   ├── FDA classification
        │   ├── FTC enforcement
        │   └── Privacy law changes
        │
        └── OPERATIONAL
            ├── Data breach
            ├── Key person risk
            └── Supplier dependency
```

### Threat Scoring Matrix

| Factor | Weight | Scale |
|--------|--------|-------|
| Likelihood | 40% | 1-5 (5 = very likely) |
| Impact | 40% | 1-5 (5 = severe) |
| Velocity | 20% | 1-5 (5 = immediate) |

**Risk Score = (Likelihood × 0.4) + (Impact × 0.4) + (Velocity × 0.2)**

---

## Competitive Threats

### C1: SelfDecode Adds Traditional Medicine

| Dimension | Assessment |
|-----------|------------|
| **Description** | SelfDecode adds TCM, Ayurveda, or herbal medicine features |
| **Likelihood** | 4 (Likely within 18 months) |
| **Impact** | 4 (Significant differentiation erosion) |
| **Velocity** | 3 (Would take 6-12 months to implement) |
| **Risk Score** | **3.8 (HIGH)** |

**Analysis:**
- SelfDecode has resources ($12M+ raised)
- They've shown willingness to expand features
- Traditional medicine is logical extension
- They would likely start with Ayurveda doshas (easier)

**Early Warning Indicators:**
- Job postings for TCM/Ayurveda specialists
- Partnership announcements with traditional medicine experts
- Feature announcements or press releases
- Social media hints from leadership

**Our Response:**
- Launch fast, establish first-mover perception
- Build depth beyond surface features
- Community moat harder to replicate than features
- Practitioner relationships as switching costs

### C2: Well-Funded Startup Enters

| Dimension | Assessment |
|-----------|------------|
| **Description** | New startup with $10M+ funding enters genetics + traditional medicine |
| **Likelihood** | 3 (Moderate) |
| **Impact** | 4 (Significant resource disadvantage) |
| **Velocity** | 2 (Would take 12-24 months to build) |
| **Risk Score** | **3.2 (MEDIUM-HIGH)** |

**Analysis:**
- Space is getting attention from investors
- Traditional medicine is trending
- Validated market reduces risk for new entrants
- But: niche market may not attract big funding

**Early Warning Indicators:**
- Funding announcements in adjacent spaces
- Y Combinator or accelerator companies in space
- Talent moving from SelfDecode/StrateGene to new ventures

**Our Response:**
- Focus on community stickiness (network effects)
- Grassroots community relationships
- Don't compete on features, compete on depth
- Consider strategic partnerships or acqui-hire defense

### C3: Big Tech Entry (Google, Apple, Amazon)

| Dimension | Assessment |
|-----------|------------|
| **Description** | Major tech company launches genetics health platform |
| **Likelihood** | 2 (Low in near term) |
| **Impact** | 5 (Existential for mass market) |
| **Velocity** | 3 (Could launch quickly with resources) |
| **Risk Score** | **3.0 (MEDIUM)** |

**Analysis:**
- Google has 23andMe investment history
- Apple expanding into health (HealthKit, Research)
- Amazon Pharmacy shows health interest
- But: genetics is regulatory minefield, not their core

**Early Warning Indicators:**
- Acquisition of genetics companies
- Senior hires in genetics/genomics
- Patent filings in personal genomics
- Healthcare product expansions

**Our Response:**
- Focus on niche (chronic illness) big tech won't prioritize
- Traditional medicine expertise hard to acquire
- Community relationships as moat
- Position for acquisition if needed

### C4: StrateGene Goes Interactive

| Dimension | Assessment |
|-----------|------------|
| **Description** | StrateGene upgrades from static PDFs to interactive platform |
| **Likelihood** | 3 (Moderate) |
| **Impact** | 2 (Limited, they're focused elsewhere) |
| **Velocity** | 2 (Would require significant rebuild) |
| **Risk Score** | **2.4 (MEDIUM)** |

**Analysis:**
- StrateGene owned by supplement company (Seeking Health)
- Business model tied to supplement sales
- Less incentive to build full platform
- Practitioner-focused, not consumer

**Our Response:**
- Beat them on consumer experience
- No supplement conflict of interest
- Traditional medicine (they have none)
- Community features (they have none)

---

## Market Threats

### M1: DTC Genetics Market Decline

| Dimension | Assessment |
|-----------|------------|
| **Description** | Consumer interest in genetic testing declines |
| **Likelihood** | 2 (Low) |
| **Impact** | 4 (Shrinks addressable market) |
| **Velocity** | 1 (Slow decline) |
| **Risk Score** | **2.4 (MEDIUM)** |

**Analysis:**
- 40M+ already tested, stable base
- Privacy concerns could slow new testing
- But: health applications growing even as ancestry plateaus
- Our model works with existing DNA, not new tests

**Early Warning Indicators:**
- 23andMe/Ancestry revenue declines
- Testing kit price increases
- Market research showing declining interest
- Regulatory restrictions on DTC testing

**Our Response:**
- Focus on health interpretation (post-test)
- Partner with testing companies for referrals
- Value from existing data, not new tests
- International expansion for growth

### M2: Market Consolidation

| Dimension | Assessment |
|-----------|------------|
| **Description** | Major acquisition consolidates competitive landscape |
| **Likelihood** | 3 (Moderate) |
| **Impact** | 3 (Changes competitive dynamics) |
| **Velocity** | 2 (Takes time to integrate) |
| **Risk Score** | **2.8 (MEDIUM)** |

**Possible Scenarios:**
- SelfDecode acquires StrateGene
- 23andMe acquires interpretation platform
- Health company (CVS, Optum) enters via acquisition

**Our Response:**
- Be acquisition target (if price is right)
- Or: focus on niche too small for acquirers
- Build defensible moat (community, data)

### M3: Pricing Pressure

| Dimension | Assessment |
|-----------|------------|
| **Description** | Competitors race to bottom on pricing |
| **Likelihood** | 3 (Moderate) |
| **Impact** | 3 (Margin compression) |
| **Velocity** | 2 (Gradual) |
| **Risk Score** | **2.8 (MEDIUM)** |

**Analysis:**
- Free tools already exist (Promethease, Genetic Genie)
- AI could commoditize interpretation
- But: quality and depth still differentiate
- Community value independent of price

**Our Response:**
- Compete on value, not price
- Community features add stickiness
- Practitioner marketplace creates non-commodity revenue
- Premium positioning with "THREE WORLDS"

---

## Technology Threats

### T1: AI Commoditization

| Dimension | Assessment |
|-----------|------------|
| **Description** | LLM capabilities become freely available, eroding AI differentiation |
| **Likelihood** | 4 (Highly likely) |
| **Impact** | 2 (Limited, AI is not our only moat) |
| **Velocity** | 3 (Happening now) |
| **Risk Score** | **3.0 (MEDIUM)** |

**Analysis:**
- OpenAI, Claude, open source models improving
- Anyone can build basic genetics chatbot
- But: our value is curated knowledge + community
- AI is interface, not core differentiator

**Our Response:**
- Differentiate on knowledge quality, not AI
- Traditional medicine expertise remains hard
- Community can't be replicated by AI alone
- Keep AI as feature, not positioning

### T2: Data Source Changes

| Dimension | Assessment |
|-----------|------------|
| **Description** | Key databases change access, pricing, or availability |
| **Likelihood** | 2 (Low) |
| **Impact** | 3 (Operational disruption) |
| **Velocity** | 2 (Would have notice period) |
| **Risk Score** | **2.2 (LOW-MEDIUM)** |

**At Risk:**
- dbSNP / NCBI (government, stable)
- PubMed (government, stable)
- TCMSP (academic, less stable)
- IMPPAT (academic, less stable)

**Our Response:**
- Diversify data sources
- Local caching / archival
- Build proprietary curated database over time
- Contribute to open databases

### T3: Platform Obsolescence

| Dimension | Assessment |
|-----------|------------|
| **Description** | Technology stack becomes outdated, requiring major rebuild |
| **Likelihood** | 2 (Low in near term) |
| **Impact** | 3 (Technical debt, slow development) |
| **Velocity** | 1 (Slow) |
| **Risk Score** | **2.0 (LOW)** |

**Our Response:**
- Use mainstream, well-supported technologies
- Avoid exotic or bleeding-edge dependencies
- Regular dependency updates
- Modular architecture for component replacement

---

## Regulatory Threats

### R1: FDA Reclassification

| Dimension | Assessment |
|-----------|------------|
| **Description** | FDA classifies genetic health platforms as medical devices |
| **Likelihood** | 2 (Low but non-zero) |
| **Impact** | 5 (Would require significant compliance or pivot) |
| **Velocity** | 2 (Regulatory process is slow) |
| **Risk Score** | **2.8 (MEDIUM)** |

**Analysis:**
- FDA has been permissive with DTC genetics
- 21st Century Cures Act provides exemptions
- But: regulation could tighten after adverse events
- Traditional medicine claims especially risky

**Early Warning Indicators:**
- FDA warning letters to competitors
- Congressional hearings on DTC genetics
- Adverse event reports linking genetics to harm
- Industry association guidance changes

**Our Response:**
- Maintain educational positioning
- Avoid diagnostic or treatment claims
- Legal review of all content
- Be able to pivot to practitioner-only if needed

### R2: FTC Enforcement

| Dimension | Assessment |
|-----------|------------|
| **Description** | FTC cracks down on health claims in genetics space |
| **Likelihood** | 3 (Moderate) |
| **Impact** | 3 (Would require content review, possible fines) |
| **Velocity** | 2 (Investigation takes time) |
| **Risk Score** | **2.8 (MEDIUM)** |

**Analysis:**
- FTC has increased health claim enforcement
- Genetic health claims are scrutinized
- Traditional medicine claims especially risky
- Competitors making bold claims could trigger industry action

**Early Warning Indicators:**
- FTC actions against health supplement companies
- Warning letters to genetics companies
- Complaints from medical associations
- Media coverage of "genetics scams"

**Our Response:**
- Conservative health claims
- All claims evidence-linked
- Legal review process
- Clear disclaimers

### R3: Privacy Law Changes

| Dimension | Assessment |
|-----------|------------|
| **Description** | New privacy laws restrict genetic data collection/use |
| **Likelihood** | 3 (Moderate) |
| **Impact** | 3 (Compliance burden, possible feature restrictions) |
| **Velocity** | 2 (Legislation takes time) |
| **Risk Score** | **2.8 (MEDIUM)** |

**Potential Changes:**
- State genetic privacy laws (beyond GINA)
- GDPR-style restrictions on genetic data
- Right to deletion extending to derived insights
- Restrictions on cross-border data transfer

**Our Response:**
- Privacy-first architecture from day one
- Exceed current requirements
- Participate in industry discussions
- Be prepared to geo-fence if needed

---

## Operational Threats

### O1: Data Breach

| Dimension | Assessment |
|-----------|------------|
| **Description** | Unauthorized access to user genetic data |
| **Likelihood** | 2 (Low with proper security) |
| **Impact** | 5 (Existential reputational damage) |
| **Velocity** | 5 (Immediate impact) |
| **Risk Score** | **3.4 (MEDIUM-HIGH)** |

**Analysis:**
- Genetic data is uniquely sensitive
- Breach would be major news story
- User trust impossible to rebuild
- Regulatory fines could be significant

**Our Response:**
- Field-level encryption for genetic data
- SOC 2 Type II certification target
- Regular security audits
- Incident response plan
- Cyber insurance

### O2: Key Person Risk

| Dimension | Assessment |
|-----------|------------|
| **Description** | Founder/key team member becomes unavailable |
| **Likelihood** | 2 (Low) |
| **Impact** | 4 (Significant disruption) |
| **Velocity** | 5 (Immediate) |
| **Risk Score** | **3.0 (MEDIUM)** |

**Our Response:**
- Document processes and decisions
- Cross-train team members
- Key person insurance consideration
- Build team depth as we scale

### O3: Supplier Dependency

| Dimension | Assessment |
|-----------|------------|
| **Description** | Critical supplier (cloud, API, data) fails or changes terms |
| **Likelihood** | 2 (Low) |
| **Impact** | 3 (Operational disruption) |
| **Velocity** | 2 (Usually have notice) |
| **Risk Score** | **2.2 (LOW-MEDIUM)** |

**Critical Dependencies:**
- Cloud provider (AWS/GCP)
- Anthropic Claude API
- Stripe payments
- NCBI/PubMed APIs

**Our Response:**
- Multi-cloud architecture consideration
- API abstraction layer for portability
- Backup payment processor
- Data caching for external APIs

---

## Threat Summary Matrix

| Threat | Score | Likelihood | Impact | Velocity | Priority |
|--------|-------|------------|--------|----------|----------|
| C1: SelfDecode adds TCM | 3.8 | HIGH | HIGH | MEDIUM | **P0** |
| O1: Data breach | 3.4 | LOW | CRITICAL | IMMEDIATE | **P0** |
| C2: Well-funded entrant | 3.2 | MEDIUM | HIGH | LOW | **P1** |
| T1: AI commoditization | 3.0 | HIGH | LOW | HIGH | **P2** |
| C3: Big tech entry | 3.0 | LOW | CRITICAL | MEDIUM | **P2** |
| O2: Key person risk | 3.0 | LOW | HIGH | IMMEDIATE | **P1** |
| R1: FDA reclassification | 2.8 | LOW | CRITICAL | LOW | **P1** |
| R2: FTC enforcement | 2.8 | MEDIUM | MEDIUM | LOW | **P2** |
| R3: Privacy law changes | 2.8 | MEDIUM | MEDIUM | LOW | **P2** |
| M2: Market consolidation | 2.8 | MEDIUM | MEDIUM | LOW | **P2** |
| M3: Pricing pressure | 2.8 | MEDIUM | MEDIUM | LOW | **P3** |
| C4: StrateGene goes interactive | 2.4 | MEDIUM | LOW | LOW | **P3** |
| M1: Market decline | 2.4 | LOW | HIGH | LOW | **P3** |
| O3: Supplier dependency | 2.2 | LOW | MEDIUM | LOW | **P3** |
| T2: Data source changes | 2.2 | LOW | MEDIUM | LOW | **P3** |
| T3: Platform obsolescence | 2.0 | LOW | MEDIUM | LOW | **P3** |

---

## Monitoring Plan

### Quarterly Threat Review

| Activity | Frequency | Owner |
|----------|-----------|-------|
| Competitor monitoring | Weekly | Marketing |
| Regulatory scanning | Monthly | Legal |
| Security assessment | Quarterly | Engineering |
| Full threat review | Quarterly | Leadership |

### Early Warning Dashboard

| Threat | Indicator | Monitoring Method |
|--------|-----------|-------------------|
| C1 (SelfDecode) | Job postings, press | Google Alerts, LinkedIn |
| C2 (New entrants) | Funding announcements | Crunchbase, TechCrunch |
| R1 (FDA) | Warning letters | FDA database, industry news |
| O1 (Breach) | Security incidents | Monitoring tools, industry reports |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [23-COMPETITORS](../20-market/23-COMPETITORS.md) | Competitive analysis |
| [25-WHITE-SPACE](../20-market/25-WHITE-SPACE.md) | Competitive response |
| [90-RISK-OVERVIEW](./90-RISK-OVERVIEW.md) | Risk framework |
| [93-MITIGATIONS](./93-MITIGATIONS.md) | Response strategies |

---

## Open Questions

- [ ] Assess probability of specific regulatory scenarios
- [ ] Determine competitor monitoring budget
- [ ] Evaluate cyber insurance coverage needs
- [ ] Plan for acquisition defense or exit scenarios

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Strategy | Complete threat analysis |
