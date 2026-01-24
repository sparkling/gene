# Business Model

**Document ID:** 50-BUSINESS-MODEL
**Status:** Final
**Owner:** Business Strategy
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

The platform operates as a **marketplace-maker**, not just a SaaS tool. Four-sided marketplace connects Clients with Knowledge, Practitioners, Creators, and Each Other. Hybrid revenue model: subscriptions ($19-179/year) provide base recurring revenue, marketplace transactions (15-30% take rate) provide scalable upside. Network effects create defensible moat. Unit economics target LTV:CAC > 3:1.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Business model | Marketplace-maker (hybrid) | Network effects, multiple revenue streams | Jan 2026 |
| Primary revenue | Subscriptions + marketplace fees | Recurring + scalable | Jan 2026 |
| B2C vs B2B | B2C first, B2B Phase 3 | Easier acquisition, marketplace economics | Jan 2026 |
| Take rate | 15-30% depending on service | Competitive with similar platforms | Jan 2026 |

---

## Business Model Canvas

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              BUSINESS MODEL CANVAS                               │
├─────────────────┬─────────────────┬─────────────────┬─────────────────┬─────────┤
│ KEY PARTNERS    │ KEY ACTIVITIES  │ VALUE PROP      │ CUSTOMER REL    │CUSTOMER │
│                 │                 │                 │                 │SEGMENTS │
│ • Data providers│ • Data curation │ Three Worlds    │ • Self-service  │         │
│   (dbSNP,TCMSP) │ • AI development│ unified:        │ • AI-powered    │•Chronic │
│ • DNA testing   │ • Platform dev  │ genetics +      │ • Community     │ illness │
│   companies     │ • Community mgmt│ traditional     │ • Practitioner  │•Biohack │
│ • Practitioners │ • Practitioner  │ medicine +      │   matching      │•Parents │
│ • Research      │   onboarding    │ nutrition       │                 │•Prevent │
│   institutions  │                 │                 │                 │•Practnr │
├─────────────────┼─────────────────┼─────────────────┴─────────────────┼─────────┤
│ KEY RESOURCES   │                 │ CHANNELS                          │         │
│                 │                 │                                   │         │
│ • Knowledge     │                 │ • Web platform                    │         │
│   graph (40+ DB)│                 │ • Community forums                │         │
│ • AI/LLM        │                 │ • Practitioner referrals          │         │
│ • Engineering   │                 │ • Content marketing               │         │
│ • Community     │                 │ • Word of mouth                   │         │
│                 │                 │                                   │         │
├─────────────────┴─────────────────┼───────────────────────────────────┴─────────┤
│ COST STRUCTURE                    │ REVENUE STREAMS                             │
│                                   │                                             │
│ • Engineering (60%)               │ • Subscriptions: $19-179/year               │
│ • AI/LLM API costs (15%)          │ • Marketplace: 15-30% take rate             │
│ • Data licensing (5%)             │ • Affiliate: 5-15% commission               │
│ • Infrastructure (10%)            │ • Data licensing (future): negotiated       │
│ • Marketing (10%)                 │                                             │
└───────────────────────────────────┴─────────────────────────────────────────────┘
```

---

## Four-Sided Marketplace

### Marketplace Architecture

```
DEMAND SIDE                              SUPPLY SIDE
━━━━━━━━━━━                              ━━━━━━━━━━━

┌─────────────┐                          ┌─────────────┐
│   CLIENTS   │                          │ KNOWLEDGE   │
│             │◄────────────────────────►│  (40+ DBs)  │
│ • Upload    │         SIDE 1           │ • Genetics  │
│   data      │    Knowledge Access      │ • TCM       │
│ • Ask       │                          │ • Ayurveda  │
│   questions │                          │ • Nutrition │
└─────────────┘                          └─────────────┘
       │                                        │
       │          ┌──────────────┐              │
       └─────────►│   PLATFORM   │◄─────────────┘
                  │              │
       ┌─────────►│  AI Engine   │◄─────────────┐
       │          │  Marketplace │              │
       │          └──────────────┘              │
       │                                        │
┌─────────────┐          SIDE 2          ┌─────────────┐
│PRACTITIONERS│    Consulting Services   │  CREATORS   │
│             │◄────────────────────────►│             │
│ • Consult   │                          │ • Reports   │
│ • Interpret │         SIDE 3           │ • Templates │
│ • Guide     │    Content Marketplace   │ • Protocols │
└─────────────┘                          └─────────────┘
       │                                        │
       │          ┌──────────────┐              │
       └─────────►│  COMMUNITY   │◄─────────────┘
                  │              │
                  │    SIDE 4    │
                  │ Peer Support │
                  │ • Protocols  │
                  │ • Shared SNPs│
                  └──────────────┘
```

### Side Descriptions

| Side | Participants | Value Exchange | Revenue Model |
|------|--------------|----------------|---------------|
| **1. Knowledge** | Clients ↔ Database | AI-powered insights from 40+ databases | Subscription |
| **2. Consulting** | Clients ↔ Practitioners | One-on-one sessions, interpretation | 15-20% take rate |
| **3. Content** | Clients ↔ Creators | Report templates, protocols | 30% take rate |
| **4. Community** | Clients ↔ Clients | Peer support, shared experiences | Subscription (access) |

---

## Why Marketplace > Pure SaaS

| Factor | Pure SaaS | Marketplace Model |
|--------|-----------|-------------------|
| **Revenue ceiling** | Limited by subscribers | Unlimited (% of all transactions) |
| **Network effects** | None | Strong (more providers → more clients → more providers) |
| **Moat** | Features (copyable) | Network (hard to replicate) |
| **Unit economics** | CAC per user | CAC amortized across transactions |
| **Exit multiple** | 5-10x revenue | 10-20x+ revenue |
| **Defensibility** | Low (features copied) | High (network locked in) |

---

## Network Effects Flywheel

```
                    ┌─────────────────────┐
                    │    MORE CLIENTS     │
                    │  Upload DNA, labs   │
                    └──────────┬──────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────┐
│                                                      │
│   More Practitioners ──► Better Matching ──► More   │
│   join platform           quality service    Clients│
│                                                      │
└──────────────────────────────────────────────────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────┐
│                                                      │
│   More Creators ──► Better Templates ──► More       │
│   make reports      and protocols       Purchases   │
│                                                      │
└──────────────────────────────────────────────────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────┐
│                                                      │
│   More Community ──► More Shared ──► More Value     │
│   members            knowledge       for everyone   │
│                                                      │
└──────────────────────────────────────────────────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │  DEFENSIBLE MOAT    │
                    │  Network effects    │
                    │  compound over time │
                    └─────────────────────┘
```

---

## Revenue Model (Hybrid)

### Revenue Streams Overview

| Stream | Model | Rate | Phase | Year 1 Target |
|--------|-------|------|-------|---------------|
| **Client Subscriptions** | Recurring | $19-179/year | Phase 1 | $200K |
| **Consulting Marketplace** | Transaction % | 15-20% | Phase 2 | $50K |
| **Report Template Sales** | Transaction % | 30% | Phase 2 | $30K |
| **Practitioner Subscriptions** | Recurring | $49-299/year | Phase 3 | $20K |
| **Affiliate** | Commission | 5-15% | Phase 2 | $20K |
| **Data Licensing** | B2B | Negotiated | Phase 4+ | — |

### Revenue Mix Evolution

```
Year 1                  Year 2                  Year 3
━━━━━━                  ━━━━━━                  ━━━━━━

┌─────────────┐        ┌─────────────┐        ┌─────────────┐
│ Subscriptions│        │ Subscriptions│        │ Subscriptions│
│     80%     │        │     50%     │        │     35%     │
├─────────────┤        ├─────────────┤        ├─────────────┤
│ Marketplace │        │ Marketplace │        │ Marketplace │
│     15%     │        │     35%     │        │     45%     │
├─────────────┤        ├─────────────┤        ├─────────────┤
│   Other     │        │   Other     │        │   Other     │
│     5%      │        │     15%     │        │     20%     │
└─────────────┘        └─────────────┘        └─────────────┘

Target: $300K          Target: $1.5M          Target: $5M
```

---

## Value Creation Logic

### For Clients (Demand Side)

| Input | Platform Process | Output Value |
|-------|------------------|--------------|
| Genetic data | Parse, analyze, correlate | Personalized insights |
| Lab results | Map to pathways | Root cause understanding |
| Symptoms | AI pattern matching | Treatment suggestions |
| Questions | RAG + Claude | Evidence-based answers |
| Need for guidance | Practitioner matching | Expert support |
| Need for community | Group matching | Peer connection |

### For Practitioners (Supply Side)

| Input | Platform Process | Output Value |
|-------|------------------|--------------|
| Expertise | Profile, visibility | Client acquisition |
| Time | Booking, payments | Revenue stream |
| Knowledge | Report creation | Passive income |
| Client data | Analytics | Better outcomes |

### For Creators

| Input | Platform Process | Output Value |
|-------|------------------|--------------|
| Domain expertise | Template creation | Scalable income |
| Protocols | Systemization | Reach more people |
| SNP collections | Curation | Community value |

---

## Competitive Moat Strategy

### Moat Components

| Moat Type | Description | Defensibility |
|-----------|-------------|---------------|
| **Network Effects** | More users → more value → more users | HIGH |
| **Data Network** | Unique TCM/Ayurveda gene connections | HIGH |
| **Switching Costs** | Health history, practitioner relationships | MEDIUM |
| **Brand** | Trust in health space | MEDIUM |
| **Content** | Community protocols, reports | MEDIUM |

### Moat Building Timeline

| Phase | Moat Focus | Actions |
|-------|------------|---------|
| Phase 1 | Data + Content | Build unique THREE WORLDS connections |
| Phase 2 | Network Effects | Seed practitioners, build community |
| Phase 3 | Switching Costs | Deepen user data, relationships |
| Phase 4 | Brand | Establish authority, research partnerships |

---

## Key Partnerships

| Partner Type | Examples | Value Exchange |
|--------------|----------|----------------|
| **Data Providers** | dbSNP, TCMSP, USDA | Data access for attribution |
| **DNA Testing** | 23andMe, AncestryDNA | User referrals both ways |
| **Practitioners** | Functional medicine networks | Clients for visibility |
| **Supplement Brands** | Quality brands (no conflicts) | Affiliate revenue |
| **Research Institutions** | Universities | Data for research credit |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [51-REVENUE-STREAMS](./51-REVENUE-STREAMS.md) | Detailed revenue model |
| [52-PRICING](./52-PRICING.md) | Pricing strategy |
| [53-MARKETPLACE](./53-MARKETPLACE.md) | Marketplace dynamics |
| [54-UNIT-ECONOMICS](./54-UNIT-ECONOMICS.md) | Unit economics |

---

## Open Questions

- [ ] Optimal take rate for practitioner consultations
- [ ] Practitioner onboarding incentives
- [ ] Community monetization approach
- [ ] B2B pricing validation

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Business Strategy | Complete business model |
