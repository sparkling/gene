# Financials Overview

**Document ID:** 80-FINANCIALS-OVERVIEW
**Status:** Final
**Owner:** Finance
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Hybrid revenue model: recurring subscriptions (core) + marketplace transactions (growth). Year 1 target: $100K MRR. Year 3 target: $1M MRR. Break-even expected Month 18-24. Initial funding requirement: $150K for MVP + launch. Unit economics: LTV $175, CAC $25, LTV:CAC 7:1. Primary revenue drivers: Gene Annual ($79/year) and practitioner consulting fees (15-20% take rate).

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Revenue model | Hybrid (subscription + marketplace) | Multiple revenue streams, network effects | Jan 2026 |
| Pricing anchor | $79/year (Gene Annual) | Sweet spot for value/accessibility | Jan 2026 |
| Take rate | 15-20% consulting, 30% templates | Industry standard, sustainable | Jan 2026 |
| Funding approach | Bootstrap → Seed | Prove PMF before raising | Jan 2026 |
| Break-even target | Month 18-24 | Sustainable growth | Jan 2026 |

---

## Financial Highlights

### Key Metrics Summary

| Metric | Year 1 | Year 2 | Year 3 |
|--------|--------|--------|--------|
| Registered Users | 10,000 | 50,000 | 200,000 |
| Paying Users | 2,000 | 12,000 | 50,000 |
| MRR (End of Year) | $100K | $400K | $1M |
| ARR | $1.2M | $4.8M | $12M |
| Gross Margin | 75% | 80% | 85% |

### Revenue Mix (Year 3 Target)

```
REVENUE COMPOSITION
│
├── Subscriptions (60%)
│   ├── Consumer ($79-179/yr)    45%
│   └── Practitioner ($49-299/yr) 15%
│
├── Marketplace (35%)
│   ├── Consulting fees          25%
│   └── Template sales           10%
│
└── Other (5%)
    ├── Affiliate                 3%
    └── Data licensing            2%
```

---

## Key Assumptions

### Market Assumptions

| Assumption | Value | Confidence | Basis |
|------------|-------|------------|-------|
| TAM (US genetics consumers) | 40M | HIGH | 23andMe + Ancestry users |
| SAM (health-focused) | 4M | MEDIUM | ~10% seeking health insights |
| SOM (Year 3) | 50,000 | MEDIUM | 1.25% of SAM |
| Market growth rate | 20% CAGR | HIGH | Industry reports |

### Conversion Assumptions

| Assumption | Value | Confidence | Basis |
|------------|-------|------------|-------|
| Visitor to signup | 5% | MEDIUM | Industry benchmark |
| Signup to free user | 80% | MEDIUM | Similar platforms |
| Free to paid conversion | 10% | LOW | To be validated |
| Annual churn rate | 30% | MEDIUM | SaaS benchmark |
| Practitioner churn | 20% | MEDIUM | B2B is stickier |

### Pricing Assumptions

| Assumption | Value | Confidence | Basis |
|------------|-------|------------|-------|
| Average consumer ARPU | $95/year | MEDIUM | Mix of tiers |
| Average practitioner ARPU | $150/year | MEDIUM | Mix of tiers |
| Consulting avg transaction | $150 | MEDIUM | Market research |
| Template avg transaction | $25 | MEDIUM | Market research |

---

## Revenue Model Detail

### Revenue Streams

| Stream | Model | Year 1 | Year 2 | Year 3 |
|--------|-------|--------|--------|--------|
| Consumer Subscriptions | Recurring | $600K | $2.4M | $6M |
| Practitioner Subscriptions | Recurring | $120K | $480K | $1.2M |
| Consulting Marketplace | Transaction | $240K | $960K | $2.4M |
| Template Marketplace | Transaction | $60K | $240K | $600K |
| Affiliate | Commission | $24K | $120K | $360K |
| Data Licensing | B2B | $0 | $96K | $240K |
| **Total Revenue** | | **$1.04M** | **$4.3M** | **$10.8M** |

### Subscription Tier Economics

| Tier | Price | Users (Y3) | Revenue (Y3) |
|------|-------|------------|--------------|
| Free | $0 | 150,000 | $0 |
| Gene Report | $19 one-time | 20,000 | $380K |
| Gene Annual | $79/year | 25,000 | $1.98M |
| Gene Pro | $179/year | 5,000 | $895K |
| **Consumer Total** | | **200,000** | **$3.25M** |
| Practitioner Basic | $49/month | 500 | $294K |
| Practitioner Pro | $149/month | 200 | $357K |
| Practitioner Enterprise | $299/month | 50 | $179K |
| **Practitioner Total** | | **750** | **$830K** |

### Marketplace Transaction Projections

| Stream | Y1 Transactions | Y2 Transactions | Y3 Transactions |
|--------|-----------------|-----------------|-----------------|
| Consulting sessions | 2,000 | 8,000 | 20,000 |
| Avg transaction value | $150 | $150 | $150 |
| Gross marketplace volume | $300K | $1.2M | $3M |
| Take rate | 20% | 20% | 20% |
| **Net consulting revenue** | **$60K** | **$240K** | **$600K** |
| Template sales | 3,000 | 12,000 | 30,000 |
| Avg template price | $25 | $25 | $25 |
| Gross template volume | $75K | $300K | $750K |
| Take rate | 30% | 30% | 30% |
| **Net template revenue** | **$22.5K** | **$90K** | **$225K** |

---

## Cost Structure

### Fixed Costs (Monthly)

| Category | Month 1-6 | Month 7-12 | Year 2 | Year 3 |
|----------|-----------|------------|--------|--------|
| Engineering salaries | $15,000 | $25,000 | $50,000 | $100,000 |
| Marketing salaries | $5,000 | $10,000 | $25,000 | $50,000 |
| Operations | $2,000 | $5,000 | $10,000 | $20,000 |
| Legal/Compliance | $2,000 | $2,000 | $5,000 | $10,000 |
| Office/Admin | $1,000 | $2,000 | $5,000 | $10,000 |
| **Total Fixed** | **$25,000** | **$44,000** | **$95,000** | **$190,000** |

### Variable Costs (Per Unit)

| Category | Cost | Basis |
|----------|------|-------|
| Infrastructure (per user/mo) | $0.50 | Hosting, CDN |
| Data APIs (per user/mo) | $0.10 | NCBI, PubMed |
| AI/LLM (per query) | $0.02 | Claude API |
| Payment processing | 2.9% + $0.30 | Stripe |
| Support (per user/mo) | $0.25 | Intercom, staffing |
| **Total variable (per user)** | **~$1.00/mo** | |

### Gross Margin Analysis

| Component | Year 1 | Year 2 | Year 3 |
|-----------|--------|--------|--------|
| Revenue | $1.04M | $4.3M | $10.8M |
| Variable costs | $260K | $860K | $1.62M |
| **Gross profit** | **$780K** | **$3.44M** | **$9.18M** |
| **Gross margin** | **75%** | **80%** | **85%** |

---

## Funding Requirements

### Use of Funds (Seed Stage)

| Category | Amount | % | Purpose |
|----------|--------|---|---------|
| Engineering | $60,000 | 40% | MVP development, data integration |
| Marketing | $30,000 | 20% | Launch, content, community |
| Data/APIs | $15,000 | 10% | Database licenses, API costs |
| Infrastructure | $15,000 | 10% | Hosting, security, CDN |
| Legal | $15,000 | 10% | Compliance, terms, incorporation |
| Buffer | $15,000 | 10% | Contingency |
| **Total** | **$150,000** | 100% | |

### Funding Milestones

```
FUNDING ROADMAP
│
├── NOW: Bootstrap ($0-50K)
│   └── Founder savings, revenue
│   └── Goal: MVP + 100 beta users
│
├── MONTH 6: Pre-Seed ($150K)
│   └── Angels, friends & family
│   └── Goal: 1,000 paying users
│
├── MONTH 18: Seed ($500K-1M)
│   └── Seed funds, angels
│   └── Goal: $50K MRR, PMF proven
│
└── MONTH 36: Series A ($3-5M)
    └── Institutional VCs
    └── Goal: $500K+ MRR, scale
```

---

## Path to Profitability

### Monthly P&L Projection

| Month | Users | MRR | Costs | Profit/(Loss) | Cumulative |
|-------|-------|-----|-------|---------------|------------|
| 1 | 100 | $2K | $25K | ($23K) | ($23K) |
| 3 | 300 | $8K | $28K | ($20K) | ($63K) |
| 6 | 500 | $25K | $35K | ($10K) | ($108K) |
| 12 | 2,000 | $100K | $70K | $30K | ($48K) |
| 18 | 5,000 | $200K | $120K | $80K | $192K |
| 24 | 12,000 | $400K | $200K | $200K | $792K |

### Break-Even Analysis

```
BREAK-EVEN PROJECTION
│
│  Revenue ────────────────────────────────────╱
│                                            ╱
│                                          ╱
│  Costs ─────────────────────────────────╱───────
│                                       ╱
│                                     ╱
│                                   ╱
│                    ┌─────────────╱─────────┐
│                    │ BREAK-EVEN │          │
│                    │  Month 18  │          │
│                    └────────────┴──────────┘
├────────────────────────────────────────────────▶
0     6      12     18     24     30     36  Months
```

**Break-even requirements:**
- ~3,000 paying users at $50 ARPU
- ~$150K MRR
- Expected: Month 18-24

---

## Unit Economics

### Customer Acquisition Cost (CAC)

| Channel | CAC | % of Acquisition |
|---------|-----|------------------|
| Community/organic | $10 | 50% |
| Content/SEO | $15 | 25% |
| Referral | $12 | 15% |
| Paid | $50 | 10% |
| **Blended CAC** | **$18** | 100% |

### Lifetime Value (LTV)

| Segment | ARPU | Lifespan | LTV |
|---------|------|----------|-----|
| Free user | $0 | - | $0 |
| Gene Report | $19 | 1x | $19 |
| Gene Annual | $79 | 2.5 years | $198 |
| Gene Pro | $179 | 3 years | $537 |
| Practitioner | $120/mo | 4 years | $5,760 |
| **Blended consumer LTV** | | | **$175** |

### LTV:CAC Analysis

| Metric | Value | Target | Status |
|--------|-------|--------|--------|
| Blended LTV | $175 | >$150 | GOOD |
| Blended CAC | $25 | <$50 | GOOD |
| LTV:CAC | 7:1 | >3:1 | EXCELLENT |
| Payback period | 4 months | <12 mo | EXCELLENT |

---

## Financial Projections Summary

### 3-Year P&L

| | Year 1 | Year 2 | Year 3 |
|---|--------|--------|--------|
| **Revenue** | $1.04M | $4.3M | $10.8M |
| COGS | $260K | $860K | $1.62M |
| **Gross Profit** | $780K | $3.44M | $9.18M |
| S&M | $300K | $1M | $2.5M |
| R&D | $250K | $800K | $2M |
| G&A | $150K | $400K | $1M |
| **Operating Expenses** | $700K | $2.2M | $5.5M |
| **EBITDA** | $80K | $1.24M | $3.68M |
| **EBITDA Margin** | 8% | 29% | 34% |

### Key Ratios

| Ratio | Year 1 | Year 2 | Year 3 |
|-------|--------|--------|--------|
| Gross margin | 75% | 80% | 85% |
| S&M as % revenue | 29% | 23% | 23% |
| R&D as % revenue | 24% | 19% | 19% |
| G&A as % revenue | 14% | 9% | 9% |
| EBITDA margin | 8% | 29% | 34% |

---

## Sensitivity Analysis

### Revenue Sensitivity

| Scenario | Conversion Rate | Year 3 Revenue | Year 3 EBITDA |
|----------|-----------------|----------------|---------------|
| Bear | 5% | $5.4M | $1.2M |
| Base | 10% | $10.8M | $3.68M |
| Bull | 15% | $16.2M | $6M |

### Churn Sensitivity

| Scenario | Annual Churn | Year 3 LTV | Year 3 Revenue |
|----------|--------------|------------|----------------|
| High churn | 40% | $140 | $8.6M |
| Base | 30% | $175 | $10.8M |
| Low churn | 20% | $220 | $13.5M |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [51-REVENUE-STREAMS](../50-business/51-REVENUE-STREAMS.md) | Revenue detail |
| [52-PRICING](../50-business/52-PRICING.md) | Pricing strategy |
| [54-UNIT-ECONOMICS](../50-business/54-UNIT-ECONOMICS.md) | Unit economics detail |
| [84-METRICS](./84-METRICS.md) | KPI definitions |

---

## Open Questions

- [ ] Validate conversion assumptions with beta data
- [ ] Determine optimal marketing spend by channel
- [ ] Assess international expansion costs
- [ ] Model practitioner marketplace economics

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Finance | Complete financial overview |
