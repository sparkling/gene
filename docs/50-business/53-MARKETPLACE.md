# Marketplace Dynamics

**Document ID:** 53-MARKETPLACE
**Status:** Final
**Owner:** Business
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Four-sided marketplace connecting clients with knowledge, practitioners, creators, and each other. Network effects create defensible moat: more clients → more practitioners → better matching → more clients. Chicken-and-egg solution: lead with self-service knowledge access (free/low-cost), then layer marketplace. Platform take rate: 15-20% consulting, 30% report templates. Long-term revenue shift from subscriptions to transaction fees as marketplace matures.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Marketplace model | Four-sided platform | Maximum network effects | Jan 2026 |
| Launch strategy | Knowledge-first, marketplace-later | Solve chicken-egg | Jan 2026 |
| Consulting take rate | 15-20% | Competitive with Upwork/Clarity | Jan 2026 |
| Report template take rate | 30% | Higher margin, lower volume | Jan 2026 |

---

## Marketplace Architecture

### Four-Sided Platform Model

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    FOUR-SIDED MARKETPLACE                                │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│                         ┌───────────────┐                               │
│                         │   KNOWLEDGE   │                               │
│                         │   (AI + DB)   │                               │
│                         └───────┬───────┘                               │
│                                 │                                        │
│                    ┌────────────┼────────────┐                          │
│                    │            │            │                          │
│                    ▼            ▼            ▼                          │
│            ┌───────────┐ ┌───────────┐ ┌───────────┐                   │
│            │  CLIENTS  │ │PRACTITIONERS│ │ CREATORS │                   │
│            │           │ │           │ │           │                   │
│            │ • Upload  │ │ • Consult │ │ • Reports │                   │
│            │ • Query   │ │ • Manage  │ │ • Protocols│                  │
│            │ • Learn   │ │ • Earn    │ │ • Sell    │                   │
│            └─────┬─────┘ └─────┬─────┘ └─────┬─────┘                   │
│                  │             │             │                          │
│                  └─────────────┼─────────────┘                          │
│                                │                                        │
│                                ▼                                        │
│                         ┌───────────┐                                   │
│                         │ COMMUNITY │                                   │
│                         │           │                                   │
│                         │ • Peers   │                                   │
│                         │ • Forums  │                                   │
│                         │ • Support │                                   │
│                         └───────────┘                                   │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### Side Definitions

| Side | Description | Value Received | Value Provided |
|------|-------------|----------------|----------------|
| **Clients** | End users seeking health insights | Knowledge, recommendations, help | Data, fees, content |
| **Practitioners** | Licensed health professionals | Clients, tools, income | Services, expertise |
| **Creators** | Report/protocol developers | Distribution, income | Templates, protocols |
| **Community** | Peer groups, condition networks | Support, shared knowledge | Engagement, content |

---

## Network Effects Analysis

### Network Effect Types

```
                    NETWORK EFFECTS
                          │
    ┌─────────────────────┼─────────────────────┐
    │                     │                     │
  SAME-SIDE          CROSS-SIDE            DATA
  EFFECTS            EFFECTS              EFFECTS
    │                     │                     │
┌───┴───┐             ┌───┴───┐            ┌───┴───┐
│More   │             │More   │            │More   │
│clients│             │practi-│            │data   │
│= more │             │tioners│            │= bet- │
│commun-│             │= bet- │            │ter AI │
│ity    │             │ter    │            │       │
│value  │             │match  │            │       │
└───────┘             └───────┘            └───────┘
```

### Network Effect Strength

| Effect Type | Participants | Strength | Defense |
|-------------|--------------|----------|---------|
| Same-side (community) | Clients ↔ Clients | HIGH | Hard to replicate |
| Cross-side (consulting) | Clients ↔ Practitioners | HIGH | Switching costs |
| Cross-side (templates) | Clients ↔ Creators | MEDIUM | Content moat |
| Data network | All → AI | HIGH | Proprietary training |

### Flywheel Dynamics

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    MARKETPLACE FLYWHEEL                                  │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│         ┌──────────────────────────────────────────────────┐            │
│         │                                                  │            │
│         ▼                                                  │            │
│    ┌─────────┐                                        ┌─────────┐       │
│    │ More    │                                        │ Better  │       │
│    │ Clients │ ─────────────────────────────────────► │ AI/Data │       │
│    └────┬────┘                                        └────┬────┘       │
│         │                                                  │            │
│         │                                                  │            │
│         ▼                                                  ▼            │
│    ┌─────────┐                                        ┌─────────┐       │
│    │ More    │                                        │ Better  │       │
│    │Practit- │ ◄───────────────────────────────────── │ Match   │       │
│    │ ioners  │                                        │ Quality │       │
│    └────┬────┘                                        └─────────┘       │
│         │                                                  ▲            │
│         │                                                  │            │
│         └──────────────────────────────────────────────────┘            │
│                                                                          │
│                    VIRTUOUS CYCLE                                       │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Marketplace Transactions

### Transaction Types

| Transaction | Buyer | Seller | Take Rate | Est. Volume |
|-------------|-------|--------|-----------|-------------|
| Consulting Session | Client | Practitioner | 15-20% | 1,000/mo (Y2) |
| Report Template | Client | Creator | 30% | 500/mo (Y2) |
| Protocol Access | Client | Creator | 30% | 300/mo (Y2) |
| SNP Collection | Client | Creator | 30% | 200/mo (Y2) |

### Transaction Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    CONSULTING TRANSACTION FLOW                           │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  1. DISCOVERY                                                            │
│  ─────────────                                                          │
│  Client searches for practitioner by:                                    │
│  • Specialty (functional medicine, TCM)                                  │
│  • Condition (ME/CFS, autoimmune)                                        │
│  • Modality (Ayurveda, Kampo)                                           │
│  • Location / Availability                                               │
│                                                                          │
│  2. BOOKING                                                              │
│  ──────────                                                             │
│  • Client selects time slot                                              │
│  • Pre-fills consultation form                                           │
│  • Optionally shares genetic report                                      │
│  • Pays session fee (held in escrow)                                     │
│                                                                          │
│  3. SESSION                                                              │
│  ─────────                                                              │
│  • Video/audio consultation                                              │
│  • Practitioner views shared health data                                 │
│  • Notes captured in platform                                            │
│                                                                          │
│  4. COMPLETION                                                           │
│  ─────────────                                                          │
│  • Session ends                                                          │
│  • Client prompted to review                                             │
│  • Funds released to practitioner (minus take rate)                      │
│  • Follow-up recommendations sent                                        │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### Take Rate Benchmarks

| Platform | Category | Take Rate |
|----------|----------|-----------|
| Upwork | Freelance | 10-20% |
| Clarity.fm | Expert calls | 15% |
| Gumroad | Digital products | 10% + fees |
| Etsy | Handmade goods | 6.5% + fees |
| App Store | Digital | 15-30% |
| **Gene Platform** | **Consulting** | **15-20%** |
| **Gene Platform** | **Templates** | **30%** |

---

## Chicken-and-Egg Strategy

### The Classic Problem

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    CHICKEN-AND-EGG PROBLEM                               │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌─────────────┐           ┌─────────────┐                              │
│  │             │           │             │                              │
│  │   CLIENTS   │           │PRACTITIONERS│                              │
│  │             │           │             │                              │
│  │  "Where are │           │ "Where are  │                              │
│  │   the       │ ◄──────► │  the        │                              │
│  │  experts?"  │           │  clients?"  │                              │
│  │             │           │             │                              │
│  └─────────────┘           └─────────────┘                              │
│                                                                          │
│                    DEADLOCK                                              │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### Our Solution: Value-First Strategy

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    CHICKEN-AND-EGG SOLUTION                              │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  PHASE 1: KNOWLEDGE-FIRST (Pre-marketplace)                              │
│  ──────────────────────────────────────────                             │
│  • Free/low-cost self-service genetic analysis                          │
│  • AI-powered Q&A without practitioners                                  │
│  • TCM/Ayurveda knowledge base browsing                                  │
│  • Build user base on knowledge value alone                              │
│                                                                          │
│              │                                                           │
│              ▼                                                           │
│                                                                          │
│  PHASE 2: SUPPLY SEEDING (Early marketplace)                             │
│  ─────────────────────────────────────────────                          │
│  • Recruit 50-100 practitioners manually                                 │
│  • Offer reduced take rate (10%) for early adopters                      │
│  • Provide free client management tools                                  │
│  • Guarantee minimum referrals                                           │
│                                                                          │
│              │                                                           │
│              ▼                                                           │
│                                                                          │
│  PHASE 3: DEMAND ACTIVATION (Growing marketplace)                        │
│  ────────────────────────────────────────────────                       │
│  • Promote practitioner profiles to user base                            │
│  • "Recommended for your genes" matching                                 │
│  • Incentivize first consultation (discounts)                            │
│  • Collect and display reviews                                           │
│                                                                          │
│              │                                                           │
│              ▼                                                           │
│                                                                          │
│  PHASE 4: FLYWHEEL (Mature marketplace)                                  │
│  ──────────────────────────────────────                                 │
│  • Organic supply growth (practitioners apply)                           │
│  • Organic demand growth (word-of-mouth)                                 │
│  • Increase take rate to standard 15-20%                                 │
│  • Launch creator marketplace (templates)                                │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### Practitioner Acquisition Tactics

| Tactic | Description | Target |
|--------|-------------|--------|
| Manual outreach | Direct contact to functional medicine practitioners | 50 practitioners |
| Conference presence | Booth at integrative medicine conferences | Awareness |
| Early adopter program | Reduced fees, guaranteed referrals | 100 practitioners |
| Tool-first | Free client management tools even without marketplace | Stickiness |
| Credential verification | Professional verification process | Trust |

---

## Quality Control

### Trust & Safety Framework

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    TRUST & SAFETY                                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  PRACTITIONERS                                                           │
│  ─────────────                                                          │
│  • License verification (manual review)                                  │
│  • Credential display (MD, ND, LAc, etc.)                               │
│  • Specialization claims validated                                       │
│  • Background check (optional, displayed)                                │
│  • Review system with moderation                                         │
│  • Response time tracking                                                │
│  • Session completion rate                                               │
│                                                                          │
│  CREATORS                                                                │
│  ────────                                                               │
│  • Content review before listing                                         │
│  • Medical accuracy check                                                │
│  • Evidence requirement (citations)                                      │
│  • User ratings and reviews                                              │
│  • Refund policy compliance                                              │
│                                                                          │
│  CLIENTS                                                                 │
│  ───────                                                                │
│  • Review guidelines enforcement                                         │
│  • No-show tracking                                                      │
│  • Payment verification                                                  │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### Review System

| Dimension | Measured | Display |
|-----------|----------|---------|
| Overall rating | 1-5 stars | Star average |
| Communication | 1-5 stars | Included in average |
| Knowledge | 1-5 stars | Included in average |
| Value | 1-5 stars | Included in average |
| Written review | Free text | Moderated display |
| Response rate | % within 24h | Badge |
| Session count | Total completed | Number shown |

---

## Marketplace Economics

### Unit Economics (Consulting)

| Metric | Value |
|--------|-------|
| Average session price | $150 |
| Platform take rate | 18% ($27) |
| Practitioner payout | $123 |
| Payment processing | 2.9% + $0.30 ($4.65) |
| Net revenue | $22.35 |
| Margin | ~83% |

### Unit Economics (Templates)

| Metric | Value |
|--------|-------|
| Average template price | $29 |
| Platform take rate | 30% ($8.70) |
| Creator payout | $20.30 |
| Payment processing | 2.9% + $0.30 ($1.14) |
| Net revenue | $7.56 |
| Margin | ~87% |

### Revenue Mix Projection

| Year | Subscriptions | Marketplace | Mix |
|------|---------------|-------------|-----|
| Y1 | 90% | 10% | Subscription-heavy |
| Y2 | 70% | 30% | Growing marketplace |
| Y3 | 50% | 50% | Balanced |
| Y5 | 30% | 70% | Marketplace-dominant |

---

## Competitive Moats

### Why Marketplace is Defensible

| Moat Type | Description | Replicability |
|-----------|-------------|---------------|
| Network effects | More users → more practitioners → more users | Very hard |
| Data effects | More usage → better AI → better service | Very hard |
| Content moat | User-generated templates, protocols | Hard |
| Trust/reviews | Established reputation system | Hard |
| Switching costs | Health data, practitioner relationships | Medium |

### Marketplace vs Pure SaaS

| Factor | Pure SaaS | Marketplace |
|--------|-----------|-------------|
| Revenue ceiling | Limited by subscribers | % of all transactions |
| Network effects | None | Strong |
| Moat | Features (copyable) | Network (hard to replicate) |
| Unit economics | CAC per user | CAC amortized |
| Exit multiple | 5-10x revenue | 10-20x+ revenue |

---

## Risk Factors

### Marketplace Risks

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Disintermediation | MEDIUM | HIGH | Stickiness via tools, data lock-in |
| Quality issues | MEDIUM | HIGH | Verification, reviews, moderation |
| Regulatory | MEDIUM | MEDIUM | Legal review, disclaimers |
| Chicken-and-egg failure | LOW | HIGH | Value-first strategy |
| Take rate pressure | MEDIUM | MEDIUM | Differentiated value, data moat |

### Disintermediation Prevention

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    DISINTERMEDIATION PREVENTION                          │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  WHY STAY ON PLATFORM:                                                   │
│                                                                          │
│  For Practitioners:                                                      │
│  • Client genetic data integration                                       │
│  • Client health history access                                          │
│  • Professional dashboard/CRM                                            │
│  • Payment processing/invoicing                                          │
│  • Review accumulation                                                   │
│  • New client discovery                                                  │
│                                                                          │
│  For Clients:                                                            │
│  • Payment protection (escrow)                                           │
│  • Review visibility                                                     │
│  • Health data integration                                               │
│  • Session records                                                       │
│  • Dispute resolution                                                    │
│  • Quality assurance                                                     │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [50-BUSINESS-MODEL](./50-BUSINESS-MODEL.md) | Business model context |
| [51-REVENUE-STREAMS](./51-REVENUE-STREAMS.md) | Revenue details |
| [52-PRICING](./52-PRICING.md) | Pricing strategy |
| [54-UNIT-ECONOMICS](./54-UNIT-ECONOMICS.md) | Economics detail |

---

## Open Questions

- [ ] Define practitioner vetting process in detail
- [ ] Determine minimum viable liquidity threshold
- [ ] Establish refund and dispute resolution policies
- [ ] Plan for international practitioner licensing

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Business | Complete marketplace specification |
