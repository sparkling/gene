# Product Roadmap

**Document ID:** 42-ROADMAP
**Status:** Final
**Owner:** Product
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Three-phase roadmap over 12 months: Phase 1 (MVP) delivers core DNA analysis, AI Q&A, and basic subscription in months 1-4. Phase 2 adds interactive pathway visualization, TCM/Ayurveda databases, and practitioner marketplace in months 5-8. Phase 3 scales with full marketplace, community features, and advanced AI in months 9-12. MVP prioritizes fastest path to value validation.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| MVP timeline | 4 months | Balance speed with quality | Jan 2026 |
| Phase 1 scope | Core analysis + AI + subscription | Fastest to value, validate demand | Jan 2026 |
| Traditional medicine timing | Phase 2 | Need user base first | Jan 2026 |
| Marketplace timing | Phase 2 (basic) / Phase 3 (full) | Chicken-egg problem | Jan 2026 |

---

## Roadmap Overview

```
PHASE 1: MVP                    PHASE 2: ENHANCEMENT           PHASE 3: SCALE
(Months 1-4)                    (Months 5-8)                   (Months 9-12)
━━━━━━━━━━━━━                   ━━━━━━━━━━━━━━━━               ━━━━━━━━━━━━━

┌─────────────┐                 ┌─────────────────┐            ┌─────────────────┐
│ DNA Upload  │                 │ Interactive     │            │ Full Marketplace│
│ SNP Analysis│                 │ Pathways        │            │ Community       │
│ Gene Reports│                 │ TCM Database    │            │ Advanced AI     │
│ AI Q&A     │─────────────────►│ Ayurveda DB    │────────────►│ Practitioner    │
│ Research   │                 │ Lab Integration │            │ Tools           │
│ Subscription│                 │ Basic Mktplace  │            │ Research Alerts │
└─────────────┘                 └─────────────────┘            └─────────────────┘

    VALIDATE                        DIFFERENTIATE                    SCALE
    DEMAND                          POSITION                         NETWORK
```

---

## Phase 1: MVP (Months 1-4)

### Objective
Validate demand with core functionality. Get users uploading DNA, receiving insights, and paying for value.

### Must Have Features

| Feature | Description | Acceptance Criteria | Priority |
|---------|-------------|---------------------|----------|
| **User Auth** | Secure login, register, 2FA | < 5s login, secure, OAuth support | Must |
| **DNA Upload** | 23andMe, AncestryDNA parsing | Parse in < 30s, handle errors | Must |
| **Core SNP Analysis** | Top 200 health-relevant SNPs | Display variant + meaning | Must |
| **Gene Reports** | Basic explanations per gene | Each gene has description | Must |
| **Research Links** | PubMed citations | ≥1 citation per SNP | Must |
| **AI Q&A** | Simple Claude-powered questions | Accurate, sourced responses | Must |
| **Dashboard** | Personalized home screen | Clear, actionable | Must |
| **Subscription** | Basic paywall (Stripe) | Process payments reliably | Must |

### Should Have (Phase 1.1)

| Feature | Description | Rationale |
|---------|-------------|-----------|
| Lab Entry | Manual input of lab results | Complete health picture |
| Symptom Tracking | Log symptoms | Correlation analysis |
| PDF Export | Downloadable report | Share with doctors |
| Basic Pathways | Static visualization | Visual differentiation |

### Timeline

| Month | Focus | Deliverables |
|-------|-------|--------------|
| **Month 1** | Foundation | Auth, database schema, DNA parser prototype |
| **Month 2** | Core Analysis | SNP analysis engine, gene reports, research links |
| **Month 3** | AI + UX | AI Q&A integration, dashboard, user flows |
| **Month 4** | Launch Prep | Subscription, testing, beta launch |

### Success Metrics (Phase 1)

| Metric | Target | Measurement |
|--------|--------|-------------|
| Beta Users | 500+ | Sign-ups |
| DNA Uploads | 200+ | Completed uploads |
| Conversion | 5%+ | Free to paid |
| AI Queries | 10+ per user | Usage tracking |
| NPS | 30+ | Survey |

---

## Phase 2: Enhancement (Months 5-8)

### Objective
Differentiate with interactive visualization, traditional medicine integration, and marketplace foundation.

### Features

| Feature | Description | Priority |
|---------|-------------|----------|
| **Interactive Pathways** | Cytoscape.js visualization | Must |
| **Lab Integration** | Manual entry, basic parsing | Should |
| **TCM Database** | 500+ herbs with gene connections | Must |
| **Ayurveda Database** | 200+ herbs, dosha typing | Should |
| **Gene-Herb Mapping** | Core differentiator feature | Must |
| **Basic Marketplace** | Practitioner directory, booking | Should |
| **Evidence Grading** | Clinical vs animal vs cell | Should |

### Timeline

| Month | Focus | Deliverables |
|-------|-------|--------------|
| **Month 5** | Visualization | Cytoscape.js pathways, interactive features |
| **Month 6** | Traditional Medicine | TCM database, herb-gene connections |
| **Month 7** | Ayurveda + Labs | Ayurveda DB, lab result entry |
| **Month 8** | Marketplace v1 | Practitioner profiles, basic booking |

### Success Metrics (Phase 2)

| Metric | Target | Measurement |
|--------|--------|-------------|
| Paid Users | 2,000+ | Subscriptions |
| Pathway Views | 5+ per user | Analytics |
| TCM Queries | 20% of users | Feature usage |
| Practitioner Sign-ups | 50+ | Onboarding |
| Retention (30-day) | 40%+ | Cohort analysis |

---

## Phase 3: Scale (Months 9-12)

### Objective
Scale with full marketplace, community features, and advanced analytics.

### Features

| Feature | Description | Priority |
|---------|-------------|----------|
| **Full Marketplace** | Reports, protocols, community | Could |
| **Community Features** | Condition groups, peer support | Could |
| **Advanced AI** | Multi-turn analysis, deep reasoning | Could |
| **Practitioner Tools** | Client management, white-label | Could |
| **Polygenic Risk Scores** | Multi-gene risk calculations | Could |
| **Research Alerts** | New studies notification | Could |
| **Custom SNP Requests** | Analyze specific rsIDs | Could |
| **Kampo Database** | Japanese herbal formulas | Could |

### Timeline

| Month | Focus | Deliverables |
|-------|-------|--------------|
| **Month 9** | Marketplace Full | Report templates, protocol library |
| **Month 10** | Community | Condition groups, peer features |
| **Month 11** | Advanced AI | Multi-turn, deeper analysis |
| **Month 12** | B2B + Scale | Practitioner tools, optimization |

### Success Metrics (Phase 3)

| Metric | Target | Measurement |
|--------|--------|-------------|
| Paid Users | 10,000+ | Subscriptions |
| Marketplace GMV | $50K+/month | Transactions |
| Community Members | 5,000+ | Active users |
| Practitioner Revenue | $10K+/month | Take rate |
| Retention (90-day) | 60%+ | Cohort analysis |

---

## Feature Dependencies

```
User Auth
    │
    ├──► DNA Upload ──► SNP Analysis ──► Gene Reports
    │                        │
    │                        └──► AI Q&A ──► Advanced AI
    │
    ├──► Lab Entry ──► Correlation Analysis
    │
    └──► Subscription ──► Premium Features
                              │
                              ├──► Pathway Viz (P2)
                              │
                              ├──► TCM Database ──► Gene-Herb Mapping
                              │
                              ├──► Ayurveda DB ──► Dosha Typing
                              │
                              └──► Marketplace ──► Practitioner Tools (P3)
                                      │
                                      └──► Community Features (P3)
```

---

## Risk Factors

| Risk | Impact | Mitigation |
|------|--------|------------|
| DNA parsing complexity | Delays Phase 1 | Start with 23andMe only |
| AI cost per query | Margins | Caching, query optimization |
| TCM data quality | Feature quality | Validate top 100 herbs first |
| Practitioner chicken-egg | Marketplace launch | Seed with initial practitioners |
| Regulatory scrutiny | Legal | Educational positioning, disclaimers |

---

## Resource Requirements

### Phase 1

| Role | FTE | Focus |
|------|-----|-------|
| Full-Stack Engineer | 2 | Core platform |
| Data Engineer | 1 | DNA parsing, database |
| AI/ML Engineer | 0.5 | RAG pipeline |
| Designer | 0.5 | UX/UI |

### Phase 2

| Role | FTE | Focus |
|------|-----|-------|
| Full-Stack Engineer | 2 | Features, marketplace |
| Data Engineer | 1 | TCM/Ayurveda data |
| AI/ML Engineer | 1 | Advanced AI |
| Designer | 1 | Visualization, marketplace |
| Content/Medical | 0.5 | Data curation |

### Phase 3

| Role | FTE | Focus |
|------|-----|-------|
| Full-Stack Engineer | 3 | Scale, optimization |
| Data Engineer | 1 | Data pipelines |
| AI/ML Engineer | 1 | Advanced features |
| Designer | 1 | Community, B2B |
| Content/Medical | 1 | Ongoing curation |
| Growth | 1 | Marketing, community |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [40-PRODUCT-OVERVIEW](./40-PRODUCT-OVERVIEW.md) | Product context |
| [41-FEATURES](./41-FEATURES.md) | Feature definitions |
| [43-DATA-SOURCES](./43-DATA-SOURCES.md) | Data availability |
| [50-BUSINESS-MODEL](../50-business/50-BUSINESS-MODEL.md) | Revenue timing |

---

## Open Questions

- [ ] Finalize top 200 SNP panel for MVP
- [ ] Validate Cytoscape.js for pathway visualization
- [ ] Determine practitioner onboarding sequence
- [ ] Assess Phase 2 TCM data quality

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Product | Complete roadmap |
