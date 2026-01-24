# KPIs and Metrics

**Document ID:** 84-METRICS
**Status:** Final
**Owner:** Product/Analytics
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

North Star Metric: Weekly Active Researchers (WAR) — users engaging with the platform weekly to research their health. Key metrics organized into four categories: User (activation, retention, engagement), Business (MRR, LTV, CAC), Product (feature adoption, AI usage), and Health Outcomes (user-reported improvements). Dashboard cadence: daily for operational, weekly for tactical, monthly for strategic.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| North Star Metric | Weekly Active Researchers | Captures value delivery, leads revenue | Jan 2026 |
| Metric framework | AARRR (Pirate Metrics) | Industry standard, comprehensive | Jan 2026 |
| Reporting cadence | Daily/Weekly/Monthly | Match decision frequency | Jan 2026 |
| Health outcomes | Self-reported | Avoid medical claims | Jan 2026 |

---

## North Star Metric

### Definition

**Weekly Active Researchers (WAR)**: Users who engage with the platform at least once per week to research their health data or explore the knowledge base.

### Why WAR?

| Alternative | Why Not |
|-------------|---------|
| MAU | Too broad, doesn't indicate value |
| Revenue/MRR | Lagging indicator |
| Daily Active | Too frequent for research use case |
| Reports generated | Doesn't capture ongoing value |

**WAR captures:**
- Active engagement (not just logged in)
- Value delivery (researching = getting value)
- Leading indicator for retention and revenue
- Aligned with mission (helping people research their health)

### WAR Calculation

```
WAR = COUNT(DISTINCT users WHERE
  actions.type IN ('view_report', 'explore_pathway', 'ask_ai', 'browse_knowledge')
  AND actions.timestamp > NOW() - 7 DAYS
)
```

### WAR Targets

| Phase | WAR Target | WAR/Paying User |
|-------|------------|-----------------|
| Beta (M1-2) | 60 | 60% |
| Launch (M3-6) | 300 | 50% |
| Growth (M7-12) | 1,500 | 45% |
| Scale (Y2) | 6,000 | 40% |

---

## AARRR Framework

### Metrics by Stage

```
AARRR FUNNEL
│
├── ACQUISITION: How do users find us?
│   ├── Website visitors
│   ├── Signup rate
│   └── Channel attribution
│
├── ACTIVATION: Do they have a great first experience?
│   ├── DNA upload completion
│   ├── First report viewed
│   └── Profile completion
│
├── RETENTION: Do they come back?
│   ├── Week 1 retention
│   ├── Month 1 retention
│   └── WAR
│
├── REVENUE: Do they pay?
│   ├── Free to paid conversion
│   ├── ARPU
│   └── LTV
│
└── REFERRAL: Do they tell others?
    ├── Referral rate
    ├── Invites sent
    └── NPS
```

---

## User Metrics

### Acquisition Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Website visitors | Unique visitors to marketing site | 50K/mo (Y1) | Daily |
| Signup rate | Visitors who create account | 5% | Weekly |
| Cost per signup | Marketing spend / signups | <$5 | Weekly |
| Channel mix | % from each acquisition channel | Diversified | Monthly |

### Activation Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| DNA upload rate | Signups who upload DNA | 70% | Daily |
| Upload success rate | Successful parses / uploads | 95% | Daily |
| First report viewed | Users who view SNP report | 80% | Daily |
| Profile completion | Users completing health profile | 60% | Weekly |
| Time to first value | Minutes from signup to "aha" | <10 min | Weekly |

### Activation Funnel

```
ACTIVATION FUNNEL
│
100% ─► Signup
│
│ 70% ─► DNA Upload
│
│ 56% ─► First Report
│
│ 45% ─► Explore Pathway
│
│ 35% ─► Ask AI Question
│
│ 25% ─► Return Visit (Day 7)
│
└──────────────────────────────────────►
```

### Retention Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Day 1 retention | Return within 24 hours | 50% | Daily |
| Week 1 retention | Return within 7 days | 40% | Weekly |
| Week 4 retention | Return within 28 days | 30% | Weekly |
| Month 3 retention | Active in month 3 | 25% | Monthly |
| WAR | Weekly active researchers | See targets | Daily |

### Retention Cohort Analysis

| Cohort | Week 1 | Week 4 | Week 12 | Week 24 |
|--------|--------|--------|---------|---------|
| Target | 40% | 30% | 25% | 20% |
| Minimum | 30% | 20% | 15% | 12% |

### Engagement Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Sessions per week | Avg sessions for active users | 2.5 | Weekly |
| Session duration | Avg time per session | 8 min | Weekly |
| Reports viewed | Reports per user per month | 3 | Monthly |
| AI queries | Questions asked per user | 5/mo | Weekly |
| Pathways explored | Unique pathways viewed | 4/mo | Monthly |

---

## Business Metrics

### Revenue Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| MRR | Monthly recurring revenue | See targets | Daily |
| ARR | Annual recurring revenue | MRR × 12 | Monthly |
| Revenue growth | MoM revenue change | 15%+ | Monthly |
| ARPU | Revenue / paying users | $95/yr | Monthly |
| Expansion revenue | Upgrades + add-ons | 20% of new | Monthly |

### MRR Targets

| Month | MRR Target | Paying Users |
|-------|------------|--------------|
| 3 | $5,000 | 100 |
| 6 | $25,000 | 500 |
| 12 | $100,000 | 2,000 |
| 18 | $200,000 | 5,000 |
| 24 | $400,000 | 12,000 |

### Conversion Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Free to paid | % free users who upgrade | 10% | Weekly |
| Trial to paid | % trial users who convert | 25% | Weekly |
| Upgrade rate | % users upgrading tier | 15% | Monthly |
| Time to conversion | Days from signup to payment | <30 | Weekly |

### Churn Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Monthly churn rate | % users canceling | <3% | Monthly |
| Annual churn rate | % users churning per year | <30% | Monthly |
| Revenue churn | Lost MRR / total MRR | <5% | Monthly |
| Net revenue retention | (MRR - churn + expansion) / MRR | >100% | Monthly |

### Unit Economics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| CAC | Customer acquisition cost | <$30 | Monthly |
| LTV | Lifetime value | >$175 | Monthly |
| LTV:CAC | Ratio | >5:1 | Monthly |
| Payback period | Months to recover CAC | <6 | Monthly |
| Gross margin | (Revenue - COGS) / Revenue | >75% | Monthly |

---

## Product Metrics

### Feature Adoption

| Feature | Definition | Target | Priority |
|---------|------------|--------|----------|
| DNA upload | % users with uploaded DNA | 70% | P0 |
| SNP reports | % users viewing reports | 80% | P0 |
| Pathway visualization | % users exploring pathways | 50% | P0 |
| AI Q&A | % users asking questions | 40% | P0 |
| TCM insights | % users viewing TCM content | 30% | P1 |
| Ayurveda insights | % users viewing Ayurveda | 25% | P1 |
| Community | % users joining community | 20% | P1 |
| Practitioner booking | % users booking consults | 10% | P2 |

### AI/LLM Metrics

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Queries per user | Avg AI questions per month | 5 | Weekly |
| Response satisfaction | Thumbs up / total responses | 80% | Daily |
| Response latency | Time to first token | <2s | Daily |
| Hallucination rate | Factually incorrect responses | <2% | Weekly |
| Citation rate | Responses with citations | >90% | Weekly |

### Search & Discovery

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Search success rate | Searches with clicked result | 70% | Daily |
| Zero result rate | Searches with no results | <5% | Daily |
| Knowledge coverage | % queries with good answers | 85% | Weekly |

### Technical Performance

| Metric | Definition | Target | Frequency |
|--------|------------|--------|-----------|
| Page load time | Time to interactive | <2s | Daily |
| API latency (p95) | 95th percentile response | <200ms | Daily |
| Error rate | % requests with errors | <1% | Daily |
| Uptime | % time available | 99.9% | Daily |

---

## Health Outcome Metrics

### Self-Reported Outcomes

*Note: These are self-reported and educational. We make no medical claims.*

| Metric | Definition | Measurement | Frequency |
|--------|------------|-------------|-----------|
| Insight gained | "I learned something useful" | Survey (1-5) | Monthly |
| Direction found | "I know what to explore next" | Survey (1-5) | Monthly |
| Symptom understanding | "I better understand my symptoms" | Survey (1-5) | Quarterly |
| Action taken | "I tried something based on insights" | Yes/No | Quarterly |
| Practitioner discussion | "I discussed findings with provider" | Yes/No | Quarterly |

### Outcome Survey Questions

**Monthly Pulse Survey (2 questions):**
1. "How useful was [Platform] in the past month?" (1-5)
2. "Would you recommend [Platform] to someone like you?" (0-10 NPS)

**Quarterly Outcome Survey (5 questions):**
1. "Have you gained new understanding of your health?" (1-5)
2. "Have you identified potential root causes to explore?" (1-5)
3. "Have you tried any new approaches based on your research?" (Yes/No)
4. "Have you discussed findings with a healthcare provider?" (Yes/No)
5. "Overall, how has [Platform] impacted your health journey?" (1-5)

### NPS Tracking

| Segment | Target NPS | Frequency |
|---------|------------|-----------|
| All users | >50 | Monthly |
| Paying users | >60 | Monthly |
| Power users (WAR >4x) | >70 | Monthly |
| Practitioners | >70 | Quarterly |

---

## Dashboard Structure

### Daily Dashboard (Operational)

| Metric | Yesterday | 7-Day Avg | Trend |
|--------|-----------|-----------|-------|
| Signups | | | |
| DNA uploads | | | |
| WAR | | | |
| AI queries | | | |
| Error rate | | | |
| Support tickets | | | |

### Weekly Dashboard (Tactical)

| Metric | This Week | Last Week | WoW | Target |
|--------|-----------|-----------|-----|--------|
| New users | | | | |
| Activation rate | | | | |
| Retention (W1) | | | | |
| MRR | | | | |
| Conversions | | | | |
| Referrals | | | | |

### Monthly Dashboard (Strategic)

| Metric | This Month | Last Month | MoM | YTD | Target |
|--------|------------|------------|-----|-----|--------|
| Total users | | | | | |
| Paying users | | | | | |
| MRR | | | | | |
| Churn | | | | | |
| LTV:CAC | | | | | |
| NPS | | | | | |
| EBITDA | | | | | |

---

## Metric Alert Thresholds

### Critical Alerts (Immediate Action)

| Metric | Threshold | Action |
|--------|-----------|--------|
| Error rate | >5% | Engineering on-call |
| API latency | >500ms | Performance investigation |
| Signup rate | <2% (vs 5% target) | Marketing review |
| Daily churn | >10 users | Customer success outreach |

### Warning Alerts (Review Required)

| Metric | Threshold | Review |
|--------|-----------|--------|
| WAR | <80% of target | Weekly product review |
| NPS | <40 | Monthly UX review |
| Activation | <60% | Weekly growth review |
| Response satisfaction | <70% | AI content review |

---

## Metric Governance

### Metric Definitions

All metrics must have:
- Clear, unambiguous definition
- Calculation methodology
- Data source
- Owner
- Target with rationale

### Metric Review Process

| Cadence | Meeting | Participants | Focus |
|---------|---------|--------------|-------|
| Daily | Standup | Ops team | Operational metrics |
| Weekly | Growth review | Product, Marketing | Funnel, conversion |
| Monthly | Business review | Leadership | Strategic metrics |
| Quarterly | Board report | Board | Overall performance |

### Data Quality

| Requirement | Implementation |
|-------------|----------------|
| Single source of truth | Data warehouse |
| Automated calculation | No manual data entry |
| Audit trail | Version history |
| Access control | Role-based permissions |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [80-FINANCIALS-OVERVIEW](./80-FINANCIALS-OVERVIEW.md) | Financial context |
| [60-GTM-OVERVIEW](../60-gtm/60-GTM-OVERVIEW.md) | Marketing metrics |
| [41-FEATURES](../40-product/41-features.md) | Feature adoption |

---

## Open Questions

- [ ] Define AI quality metrics beyond satisfaction
- [ ] Establish health outcome measurement ethics
- [ ] Determine cohort analysis granularity
- [ ] Assess need for experimentation platform

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Analytics | Complete metrics framework |
