# Feature Inventory

**Document ID:** 41-FEATURES
**Status:** Final
**Owner:** Product
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Platform features span 6 categories: Core Platform, Knowledge Access (THREE WORLDS), Marketplace Services, Community Features, Practitioner Tools, and Advanced Analysis. MVP focuses on DNA upload, core SNP analysis, basic AI, and subscription. Unique features no competitor offers: TCM/Ayurveda/Kampo gene-treatment mapping, interactive pathway visualization, four-sided marketplace.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| MVP scope | Core analysis + AI + subscription | Fastest to value, validate demand | Jan 2026 |
| Visualization library | Cytoscape.js | Best graph visualization, React integration | Jan 2026 |
| Traditional medicine priority | TCM first, Ayurveda second | Larger data sources available | Jan 2026 |
| Marketplace timing | Phase 2 | Need user base first | Jan 2026 |

---

## Feature Categories

| Category | Description | MVP Scope |
|----------|-------------|-----------|
| **Core Platform** | Base functionality everyone uses | Yes |
| **Knowledge Access** | THREE WORLDS databases | Partial |
| **Marketplace** | Transactions between users | No |
| **Community** | Peer connections and sharing | No |
| **Practitioner Tools** | B2B features | No |
| **Advanced Analysis** | Premium analytics | Partial |

---

## Category 1: Core Platform Features

### What Users Can DO

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **User Authentication** | Secure login, 2FA | MVP | Must |
| **DNA File Upload** | 23andMe, AncestryDNA, VCF | MVP | Must |
| **Profile Management** | Demographics, goals, preferences | MVP | Must |
| **Dashboard** | Personalized home screen | MVP | Must |
| **Subscription Management** | Tiers, billing, upgrade | MVP | Must |
| **Lab Result Entry** | Manual input, file parsing | Phase 2 | Should |
| **Symptom Tracking** | Log and correlate | Phase 2 | Should |
| **Treatment History** | What worked, what didn't | Phase 2 | Should |
| **Progress Tracking** | Changes over time | Phase 2 | Could |
| **Data Export** | Download your data | Phase 2 | Should |
| **Sharing Controls** | Privacy, permissions | Phase 2 | Should |

---

## Category 2: Knowledge Access (THREE WORLDS)

### World 1: Modern Genetics

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **SNP Analysis** | Parse and interpret variants | MVP | Must |
| **Gene Reports** | What each gene does | MVP | Must |
| **Variant Impacts** | rsID explanations | MVP | Must |
| **Research Links** | PubMed citations | MVP | Must |
| **Pharmacogenomics** | Drug interactions | MVP | Should |
| **Polygenic Risk Scores** | Multi-gene risk | Phase 2 | Could |
| **SNP Imputation** | Expand to more variants | Phase 3 | Could |

### World 2: Traditional Medicine

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **TCM Herb Database** | 500+ herbs, gene connections | Phase 2 | Should |
| **TCM Pattern Matching** | Constitutional analysis | Phase 2 | Could |
| **Ayurveda Herbs** | 200+ herbs, gene connections | Phase 2 | Should |
| **Dosha Typing** | Vata/Pitta/Kapha | Phase 2 | Should |
| **Kampo Formulas** | Japanese herbal | Phase 2 | Could |
| **Western Herbal** | European/American herbs | Phase 3 | Could |
| **Gene-Herb-Condition Mapping** | Core differentiator | Phase 2 | Must |

### World 3: Nutritional Science

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **Nutrient Database** | USDA integration | MVP | Should |
| **Food-Gene Connections** | Nutrigenomics | Phase 2 | Should |
| **Supplement Analysis** | Form recommendations | Phase 2 | Could |
| **Diet Recommendations** | Personalized nutrition | Phase 3 | Could |

---

## Category 3: Marketplace Services

### What Users Can BUY/SELL

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **Practitioner Directory** | Find providers | Phase 2 | Should |
| **Consultation Booking** | Schedule sessions | Phase 2 | Should |
| **Report Templates** | Purchase condition reports | Phase 2 | Could |
| **Protocol Library** | Buy/sell protocols | Phase 3 | Could |
| **SNP Collection Sharing** | Curated panels | Phase 3 | Could |
| **White-Label Reports** | For practitioners | Phase 3 | Could |

### Transaction Models

| Model | Take Rate | Applies To |
|-------|-----------|------------|
| Consultation fees | 15-20% | Practitioner sessions |
| Report sales | 30% | Template purchases |
| Protocol sales | 30% | Protocol marketplace |

---

## Category 4: Community Features

### Who Users Can CONNECT With

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **Condition Groups** | ME/CFS, ADHD, autoimmune | Phase 2 | Should |
| **Geographic Groups** | Local meetups | Phase 3 | Could |
| **Interest Groups** | Methylation, nootropics | Phase 3 | Could |
| **Peer Support** | Share experiences | Phase 2 | Should |
| **Anonymous Sharing** | Privacy-preserving protocols | Phase 2 | Should |
| **Discussion Forums** | Q&A, advice | Phase 3 | Could |

---

## Category 5: Practitioner Tools (B2B)

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **Client Dashboard** | Manage patients | Phase 3 | Could |
| **White-Label Reports** | Branded outputs | Phase 3 | Could |
| **Protocol Templates** | Reusable protocols | Phase 3 | Could |
| **Client Invitations** | Onboard clients | Phase 3 | Could |
| **Analytics** | Track outcomes | Phase 3 | Could |

---

## Category 6: Advanced Analysis

| Feature | Description | Phase | Priority |
|---------|-------------|-------|----------|
| **AI-Powered Q&A** | Natural language queries | MVP | Must |
| **Pathway Visualization** | Interactive diagrams | Phase 2 | Must |
| **Cross-Modality Analysis** | TCM + genetics + nutrition | Phase 2 | Must |
| **Evidence Grading** | Clinical vs animal vs cell | Phase 2 | Should |
| **Research Alerts** | New studies notification | Phase 3 | Could |
| **Custom SNP Requests** | Analyze specific rsIDs | Phase 3 | Could |

---

## MVP Features (Phase 1)

### Must Have

| Feature | Description | Acceptance Criteria |
|---------|-------------|---------------------|
| User Auth | Login, register, 2FA | < 5s login, secure |
| DNA Upload | 23andMe, AncestryDNA | Parse in < 30s |
| Core SNP Analysis | Top 200 SNPs | Display variant + meaning |
| Gene Reports | Basic explanations | Each gene has description |
| Research Links | PubMed citations | ≥1 citation per SNP |
| AI Q&A | Simple questions | Accurate, sourced responses |
| Dashboard | User home | Personalized, clear |
| Subscription | Paywall | Stripe integration |

### Should Have (Phase 1.1)

| Feature | Description | Rationale |
|---------|-------------|-----------|
| Lab Entry | Manual input | Complete health picture |
| Symptom Tracking | Log symptoms | Correlation analysis |
| PDF Export | Downloadable report | Share with doctors |
| Basic Pathways | Static visualization | Visual differentiation |

---

## Feature Prioritization Matrix

| Feature | User Value (40%) | Differentiation (30%) | Effort (30%) | **Score** |
|---------|------------------|----------------------|--------------|-----------|
| AI Q&A | 10 | 8 | 7 | **8.5** |
| Pathway Visualization | 9 | 10 | 5 | **8.0** |
| TCM Database | 8 | 10 | 6 | **7.9** |
| Gene-Herb Mapping | 9 | 10 | 5 | **8.0** |
| Marketplace | 8 | 9 | 4 | **7.0** |
| Community | 7 | 7 | 6 | **6.7** |
| Polygenic Risk | 6 | 5 | 5 | **5.4** |

---

## Build vs Buy Decisions

| Component | Decision | Rationale |
|-----------|----------|-----------|
| **Auth** | Buy (Auth0/Clerk) | Security, compliance |
| **Payments** | Buy (Stripe) | Industry standard |
| **DNA Parsing** | Build | Core competency |
| **Graph DB** | Buy (Neo4j) | Mature, proven |
| **Vector Search** | Buy (Pinecone/pgvector) | Specialized |
| **LLM** | API (Claude) | Best quality |
| **Pathway Viz** | Build (Cytoscape.js) | Custom needs |
| **Data Ingestion** | Build | Unique sources |

---

## Feature Dependencies

```
User Auth
    │
    ├──► DNA Upload ──► SNP Analysis ──► Gene Reports
    │                        │
    │                        └──► AI Q&A
    │
    ├──► Lab Entry ──► Correlation Analysis
    │
    └──► Subscription ──► Premium Features
                              │
                              ├──► Pathway Viz
                              │
                              ├──► TCM Database ──► Gene-Herb Mapping
                              │
                              └──► Marketplace ──► Practitioner Tools
```

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [40-PRODUCT-OVERVIEW](./40-PRODUCT-OVERVIEW.md) | Context |
| [42-ROADMAP](./42-ROADMAP.md) | Timeline |
| [43-DATA-SOURCES](./43-DATA-SOURCES.md) | Data availability |
| [44-ARCHITECTURE](./44-ARCHITECTURE.md) | Technical feasibility |

---

## Open Questions

- [ ] Finalize top 200 SNP panel for MVP
- [ ] Assess AI cost per query
- [ ] Validate TCM data quality

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Product | Complete feature inventory |
