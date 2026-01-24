# Product Overview

**Document ID:** 40-PRODUCT-OVERVIEW
**Status:** Final
**Owner:** Product
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

The platform is a personalized health intelligence system bridging THREE WORLDS: modern genetics, traditional medicine, and nutritional science. It operates as both a self-service research tool (AI-powered knowledge exploration) and a four-sided marketplace (practitioners, creators, community). Core differentiation: no competitor connects genetics to TCM/Ayurveda/Kampo TREATMENTS with interactive pathway visualization.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Platform architecture | Graph DB + Vector DB + LLM | Relationships, semantic search, natural language | Jan 2026 |
| Visualization approach | Interactive Cytoscape.js | Build on StrateGene concept, exceed it | Jan 2026 |
| AI interface | Claude RAG | Best reasoning, privacy-focused | Jan 2026 |
| Data strategy | Open sources first | CC0/CC BY licenses, no vendor lock-in | Jan 2026 |

---

## Product Vision

**Enable anyone to understand their unique biology and find what works for their body by connecting their health data to the world's collective healing wisdom.**

### Vision Pillars

| Pillar | Description |
|--------|-------------|
| **Personalized** | Your data, your biology, your recommendations |
| **Comprehensive** | THREE WORLDS unified (genetics, traditional medicine, nutrition) |
| **Accessible** | Natural language interface, affordable pricing |
| **Evidence-linked** | Every insight backed by research citations |
| **Community-powered** | Peer protocols, shared discoveries |

---

## Core Value Delivery

### Dual-Mode Operation

| Mode | Description | Pricing | User Need |
|------|-------------|---------|-----------|
| **Self-Service Research** | Explore unified knowledge library with your own health data | Free / $79-179/year | "I want to understand my own biology" |
| **Marketplace** | Connect with practitioners, buy reports, join communities | Transaction fees (15-30%) | "I want expert help or community support" |

### Value Delivery Framework

```
USER INPUTS                    PLATFORM                      OUTPUTS
━━━━━━━━━━━                    ━━━━━━━━                      ━━━━━━━

┌─────────────┐              ┌──────────────┐              ┌─────────────┐
│ Genetic Data│──────────────►│              │──────────────►│Personalized │
│ (23andMe,   │              │   UNIFIED    │              │  Reports    │
│  Ancestry)  │              │  KNOWLEDGE   │              │             │
├─────────────┤              │    GRAPH     │              ├─────────────┤
│ Lab Results │──────────────►│              │──────────────►│AI-Powered  │
│ (blood,     │              │ 40+ databases│              │  Insights   │
│  hormones)  │              │              │              │             │
├─────────────┤              │  ┌────────┐  │              ├─────────────┤
│ Symptoms    │──────────────►│  │   AI   │  │──────────────►│Interactive │
│             │              │  │ Engine │  │              │  Pathways   │
├─────────────┤              │  └────────┘  │              │             │
│ Lifestyle   │──────────────►│              │──────────────►├─────────────┤
│ (diet,sleep)│              │  ┌────────┐  │              │Practitioner │
├─────────────┤              │  │Market- │  │──────────────►│  Matching   │
│ Treatment   │──────────────►│  │ place  │  │              │             │
│ History     │              │  └────────┘  │              ├─────────────┤
└─────────────┘              └──────────────┘              │ Community   │
                                                           │  Protocols  │
                                                           └─────────────┘
```

---

## Platform Architecture Overview

### Technical Stack

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         Frontend (Next.js/React)                        │
│  Dashboard │ Data Upload │ Reports │ Pathways │ Marketplace │ Community │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
┌─────────────────────────────────────────────────────────────────────────┐
│                          API Layer (Node/Python)                        │
│  Auth │ Data Parsers │ Analysis Engine │ AI/LLM │ Search │ Marketplace │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
┌──────────────────┬────────────────────────┬─────────────────────────────┐
│   PostgreSQL     │   Neo4j (Graph DB)     │   Vector DB (Pinecone)      │
│   Users, orders  │   Genes, SNPs,         │   Research embeddings,      │
│   subscriptions  │   pathways, nutrients, │   semantic search,          │
│   health data    │   herbs, relationships │   similarity matching       │
└──────────────────┴────────────────────────┴─────────────────────────────┘
                                    │
┌─────────────────────────────────────────────────────────────────────────┐
│                           AI/ML Layer                                    │
│  LLM (Claude) │ RAG Pipeline │ Recommendation Engine │ Embeddings       │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
┌─────────────────────────────────────────────────────────────────────────┐
│                      External Data Sources (THREE WORLDS)                │
│  dbSNP │ ClinVar │ PharmGKB │ SNPedia │ Reactome │ WikiPathways         │
│  TCMSP │ IMPPAT │ KampoDB │ HERB │ USDA FoodData │ FooDB               │
└─────────────────────────────────────────────────────────────────────────┘
```

### Database Architecture

| Database | Purpose | Technology |
|----------|---------|------------|
| **Relational** | Users, subscriptions, orders | PostgreSQL |
| **Graph** | Genes, pathways, relationships | Neo4j |
| **Vector** | Semantic search, embeddings | Pinecone/pgvector |
| **Cache** | Session, frequently accessed | Redis |
| **Object** | File uploads, reports | S3/Cloudflare R2 |

---

## Feature Prioritization Approach

### MoSCoW Framework

| Priority | Definition | Phase |
|----------|------------|-------|
| **Must Have** | Core functionality for MVP | Phase 1 |
| **Should Have** | Important but not critical | Phase 1.1 |
| **Could Have** | Nice to have if time | Phase 2 |
| **Won't Have Yet** | Deferred to future | Phase 3+ |

### Value vs Effort Matrix

```
                    LOW ←─── Development Effort ───→ HIGH

    HIGH    ┌─────────────────────────────────────────────────────┐
            │                                                     │
    User    │  Quick Wins        │         Big Bets              │
    Value   │  • SNP Analysis    │         • AI Interface        │
            │  • Basic Reports   │         • Pathways Viz        │
            │                    │         • Marketplace         │
            ├────────────────────┼─────────────────────────────────┤
            │                    │                                │
            │  Fill-Ins          │         Money Pits            │
            │  • PDF Export      │         • Mobile Apps         │
            │  • Sharing         │         • Wearables           │
    LOW     │                    │                                │
            └─────────────────────────────────────────────────────┘
```

---

## Roadmap Summary

### Phase 1: MVP (Months 1-4)

| Feature | Description | Priority |
|---------|-------------|----------|
| DNA Upload | 23andMe, AncestryDNA parsing | Must |
| Core SNP Analysis | Top 200 health-relevant SNPs | Must |
| Basic Reports | Gene/pathway explanations | Must |
| Research Links | PubMed citations per SNP | Must |
| AI Insights | Simple Claude-powered Q&A | Must |
| Subscription | Basic paywall | Must |

### Phase 2: Enhancement (Months 5-8)

| Feature | Description | Priority |
|---------|-------------|----------|
| Interactive Pathways | Cytoscape.js visualization | Should |
| Lab Integration | Manual entry, parsing | Should |
| TCM Database | Herb-gene connections | Should |
| Ayurveda Database | Treatments + doshas | Should |
| Practitioner Marketplace | Basic matching | Should |

### Phase 3: Scale (Months 9-12)

| Feature | Description | Priority |
|---------|-------------|----------|
| Full Marketplace | Reports, protocols, community | Could |
| Community Features | Peer support, sharing | Could |
| Advanced AI | Multi-turn analysis | Could |
| Practitioner Tools | Client management | Could |

---

## Competitive Differentiation

| Feature | Competitors | Us |
|---------|-------------|-----|
| SNP Coverage | SelfDecode: 200M+, StrateGene: 200 | 1M+ (quality over quantity) |
| Pathway Viz | StrateGene: static PDF | Interactive Cytoscape.js |
| TCM Integration | NONE | Full (500+ herbs) |
| Ayurveda | ADNTRO: doshas only | Full treatments |
| AI Interface | SelfDecode: DecodyGPT | Claude RAG |
| Marketplace | NONE | Four-sided |
| Evidence Links | Promethease: basic | Every insight cited |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [41-FEATURES](./41-features.md) | Feature inventory |
| [42-ROADMAP](./42-roadmap.md) | Detailed roadmap |
| [43-DATA-SOURCES](./43-data-sources.md) | THREE WORLDS data |
| [44-ARCHITECTURE](./44-architecture.md) | Technical details |

---

## Open Questions

- [ ] Finalize pathway visualization library (Cytoscape vs D3)
- [ ] AI cost optimization strategy
- [ ] SNP imputation approach

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Product | Complete product overview |
