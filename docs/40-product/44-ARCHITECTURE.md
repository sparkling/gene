# Technical Architecture

**Document ID:** 44-ARCHITECTURE
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Three-tier architecture with Next.js frontend, Node/Python API layer, and hybrid database backend (PostgreSQL for relational, Neo4j for graph relationships, Pinecone/pgvector for semantic search). AI layer powered by Claude RAG for natural language Q&A. Core innovation: graph database connecting genes → pathways → nutrients → herbs → conditions across THREE WORLDS. Security: SOC 2 Type II target, HIPAA-ready, privacy-first design.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Frontend framework | Next.js 14 + React | SSR, performance, ecosystem | Jan 2026 |
| Primary database | PostgreSQL | Mature, ACID, pgvector support | Jan 2026 |
| Graph database | Neo4j | Best for relationship queries | Jan 2026 |
| Vector database | Pinecone (cloud) / pgvector (self-hosted) | Semantic search at scale | Jan 2026 |
| AI provider | Anthropic Claude | Best reasoning, RAG performance | Jan 2026 |
| Visualization | Cytoscape.js | Interactive pathways, touch-friendly | Jan 2026 |

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         Frontend (Next.js 14/React)                     │
│  Dashboard │ Data Upload │ Reports │ Pathways │ Marketplace │ Learning │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                          API Layer (Node.js/Python)                     │
│  Auth │ Data Parsers │ Analysis Engine │ AI/LLM │ Search │ Marketplace │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        ▼                           ▼                           ▼
┌──────────────────┐   ┌────────────────────────┐   ┌─────────────────────┐
│   PostgreSQL     │   │   Neo4j (Graph DB)     │   │   Vector DB         │
│   Users, orders  │   │   Genes, SNPs,         │   │   Research embed-   │
│   subscriptions  │   │   pathways, nutrients, │   │   dings, semantic   │
│   health data    │   │   herbs, relationships │   │   search, similar-  │
│   (pgvector)     │   │                        │   │   ity matching      │
└──────────────────┘   └────────────────────────┘   └─────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                           AI/ML Layer                                    │
│  LLM (Claude) │ RAG Pipeline │ Recommendation Engine │ Embeddings       │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                      External Data Sources                               │
│  dbSNP │ ClinVar │ PubMed │ KEGG │ PharmGKB │ SNPedia │ DrugBank        │
│  TCMSP │ IMPPAT │ KampoDB │ HERB │ USDA FoodData │ FooDB               │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Technology Stack

### Frontend

| Layer | Technology | Purpose |
|-------|------------|---------|
| Framework | Next.js 14 (App Router) | SSR, routing, API routes |
| UI Library | React 18 | Component architecture |
| State | Zustand + React Query | Client state + server cache |
| Styling | Tailwind CSS + Radix UI | Utility-first + accessible components |
| Visualization | Cytoscape.js | Interactive pathway diagrams |
| Charts | Recharts / Visx | Data visualization |
| Forms | React Hook Form + Zod | Form handling + validation |

### Backend

| Layer | Technology | Purpose |
|-------|------------|---------|
| Runtime | Node.js 20 LTS | Primary API server |
| Framework | Fastify / Express | HTTP handling |
| Python Services | FastAPI | ML/Analysis microservices |
| Queue | BullMQ + Redis | Background job processing |
| Cache | Redis | Session, rate limiting, cache |
| File Storage | S3 / Cloudflare R2 | User uploads, exports |

### Databases

| Database | Use Case | Scale Target |
|----------|----------|--------------|
| PostgreSQL 15 | Users, orders, health data | 10M+ records |
| pgvector | Embedding search (self-hosted) | 1M+ vectors |
| Neo4j 5 | Knowledge graph relationships | 100M+ nodes |
| Pinecone | Research embedding search (cloud) | 10M+ vectors |
| Redis | Cache, sessions, queues | High throughput |

### AI/ML

| Component | Technology | Purpose |
|-----------|------------|---------|
| LLM | Claude 3.5 Sonnet | Q&A, report generation |
| Embeddings | OpenAI text-embedding-3-small | Document embeddings |
| RAG | LangChain / LlamaIndex | Retrieval pipeline |
| ML | scikit-learn, PyTorch | Custom models |

### Infrastructure

| Layer | Technology | Purpose |
|-------|------------|---------|
| Hosting | Vercel (frontend) + Railway/Render (backend) | PaaS deployment |
| CDN | Cloudflare | Edge caching, DDoS protection |
| Monitoring | Datadog / Sentry | APM, error tracking |
| Logging | Axiom / Papertrail | Log aggregation |
| CI/CD | GitHub Actions | Automated deployment |

---

## System Components

### 1. Data Ingestion Layer

```
┌─────────────────────────────────────────────────────────────────┐
│                    DATA INGESTION PIPELINE                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐     │
│  │ 23andMe  │   │Ancestry  │   │  VCF     │   │  Manual  │     │
│  │  Parser  │   │ Parser   │   │ Parser   │   │  Entry   │     │
│  └────┬─────┘   └────┬─────┘   └────┬─────┘   └────┬─────┘     │
│       │              │              │              │            │
│       └──────────────┴──────────────┴──────────────┘            │
│                              │                                   │
│                              ▼                                   │
│                    ┌──────────────────┐                         │
│                    │   Normalizer     │                         │
│                    │   & Validator    │                         │
│                    └────────┬─────────┘                         │
│                             │                                    │
│              ┌──────────────┼──────────────┐                    │
│              ▼              ▼              ▼                    │
│       ┌──────────┐   ┌──────────┐   ┌──────────┐               │
│       │ dbSNP    │   │ SNP      │   │ User     │               │
│       │ Lookup   │   │ Store    │   │ Profile  │               │
│       └──────────┘   └──────────┘   └──────────┘               │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**Supported Formats:**
- 23andMe v5 (txt)
- AncestryDNA (txt)
- VCF/gVCF
- PLINK BED/BIM/FAM

### 2. Knowledge Graph Layer

```
┌─────────────────────────────────────────────────────────────────┐
│                    KNOWLEDGE GRAPH (Neo4j)                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│    [Gene]──has_variant──▶[SNP]                                  │
│       │                    │                                     │
│       │                    │                                     │
│    involved_in         affects                                   │
│       │                    │                                     │
│       ▼                    ▼                                     │
│   [Pathway]──────────▶[Nutrient]                                │
│       │                    │                                     │
│       │                    │                                     │
│    treats              derived_from                              │
│       │                    │                                     │
│       ▼                    ▼                                     │
│   [Condition]◀────────[Herb/Formula]                            │
│       │                    │                                     │
│       │                    │                                     │
│    has_research        in_tradition                              │
│       │                    │                                     │
│       ▼                    ▼                                     │
│   [Publication]       [Modality]                                │
│                     (TCM/Ayur/Kampo)                            │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**Node Types:**
| Node | Count (Target) | Sources |
|------|----------------|---------|
| Gene | ~20,000 | HGNC, Ensembl |
| SNP | ~1.2M (curated) | dbSNP, ClinVar |
| Pathway | ~2,000 | Reactome, KEGG |
| Nutrient | ~1,000 | USDA, FooDB |
| Herb | ~2,000 | TCMSP, IMPPAT, KampoDB |
| Condition | ~5,000 | MeSH, ICD-10 |
| Publication | ~500,000 | PubMed |

### 3. AI/RAG Layer

```
┌─────────────────────────────────────────────────────────────────┐
│                       RAG PIPELINE                               │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────┐                                               │
│  │ User Query   │                                               │
│  └──────┬───────┘                                               │
│         │                                                        │
│         ▼                                                        │
│  ┌──────────────┐   ┌──────────────┐   ┌──────────────┐        │
│  │ Query        │   │ User Context │   │ Knowledge    │        │
│  │ Embedding    │──▶│ Retrieval    │──▶│ Retrieval    │        │
│  └──────────────┘   └──────────────┘   └──────────────┘        │
│                                               │                  │
│                                               ▼                  │
│                                        ┌──────────────┐         │
│                                        │ Context      │         │
│                                        │ Assembly     │         │
│                                        └──────┬───────┘         │
│                                               │                  │
│                                               ▼                  │
│                                        ┌──────────────┐         │
│                                        │ Claude LLM   │         │
│                                        │ Generation   │         │
│                                        └──────┬───────┘         │
│                                               │                  │
│                                               ▼                  │
│                                        ┌──────────────┐         │
│                                        │ Response +   │         │
│                                        │ Citations    │         │
│                                        └──────────────┘         │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**RAG Context Sources:**
1. User's genetic profile (SNPs, interpretations)
2. User's health data (symptoms, labs, conditions)
3. Knowledge graph (pathways, nutrients, herbs)
4. Research literature (PubMed embeddings)

### 4. Visualization Layer

```
┌─────────────────────────────────────────────────────────────────┐
│                    CYTOSCAPE.JS PATHWAYS                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Features:                                                       │
│  • Interactive zoom/pan                                          │
│  • Personal SNP overlay (red/yellow/green)                       │
│  • Click-to-expand nodes                                         │
│  • Filterable by category                                        │
│  • Touch-friendly (mobile)                                       │
│  • Publication-quality PNG/SVG export                            │
│                                                                  │
│  Pathway Templates:                                              │
│  • Methylation cycle                                             │
│  • BH4/Biopterin cycle                                           │
│  • Transsulfuration                                              │
│  • Neurotransmitter synthesis                                    │
│  • Detoxification (Phase I/II)                                   │
│  • Histamine metabolism                                          │
│  • Folate cycle                                                  │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Data Flow

### Upload Flow

```
User Upload → Parser → Validator → Normalizer → dbSNP Lookup
    → SNP Classification → Report Generation → Notification
```

### Query Flow

```
User Question → Embedding → Vector Search → Context Assembly
    → Claude RAG → Response Generation → Citation Linking
```

### Recommendation Flow

```
User Profile + SNPs → Pathway Analysis → Knowledge Graph Query
    → Recommendation Engine → Evidence Scoring → Ranked Results
```

---

## Security Architecture

### Defense in Depth

```
┌─────────────────────────────────────────────────────────────────┐
│                    SECURITY LAYERS                               │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Layer 1: Edge                                                   │
│  ├── Cloudflare WAF                                             │
│  ├── DDoS protection                                            │
│  └── Rate limiting                                              │
│                                                                  │
│  Layer 2: Application                                            │
│  ├── JWT authentication                                         │
│  ├── RBAC authorization                                         │
│  ├── Input validation (Zod)                                     │
│  └── CSRF protection                                            │
│                                                                  │
│  Layer 3: Data                                                   │
│  ├── Encryption at rest (AES-256)                               │
│  ├── Encryption in transit (TLS 1.3)                            │
│  ├── Field-level encryption (genetic data)                      │
│  └── Audit logging                                              │
│                                                                  │
│  Layer 4: Infrastructure                                         │
│  ├── VPC isolation                                              │
│  ├── Private subnets                                            │
│  └── Security groups                                            │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### Compliance Requirements

| Standard | Requirement | Status |
|----------|-------------|--------|
| HIPAA | PHI protection, BAA | Target |
| GDPR | EU data rights, consent | Target |
| CCPA | California privacy | Target |
| GINA | Genetic non-discrimination | Compliant |
| SOC 2 Type II | Security controls | Target |

---

## Scalability

### Horizontal Scaling

| Component | Scaling Strategy |
|-----------|------------------|
| API Servers | Auto-scaling behind load balancer |
| PostgreSQL | Read replicas, connection pooling |
| Neo4j | Read replicas, query routing |
| Redis | Cluster mode |
| Background Jobs | Horizontal worker scaling |

### Performance Targets

| Metric | Target |
|--------|--------|
| Page load time | <2s |
| API response time (p95) | <200ms |
| Report generation | <30s |
| AI Q&A response | <5s |
| Concurrent users | 10,000+ |

---

## Integration Points

### External APIs

| Service | Purpose | Rate Limit |
|---------|---------|------------|
| dbSNP (NCBI) | SNP reference data | 3 req/s |
| PubMed (NCBI) | Research literature | 3 req/s |
| ClinVar | Clinical significance | 3 req/s |
| Anthropic Claude | LLM generation | 60 req/min |
| OpenAI | Embeddings | 3,000 req/min |
| Stripe | Payments | Standard |
| SendGrid | Email | Standard |

### Webhooks

| Event | Destination |
|-------|-------------|
| Upload complete | User notification |
| Report ready | User notification |
| Subscription change | Billing system |
| New practitioner signup | Admin notification |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [45-DATA-MODEL](./45-DATA-MODEL.md) | Entity definitions |
| [43-DATA-SOURCES](./43-DATA-SOURCES.md) | External data sources |
| [73-TECHNOLOGY](../70-operations/73-TECHNOLOGY.md) | Tech operations |
| [74-COMPLIANCE](../70-operations/74-COMPLIANCE.md) | Security compliance |

---

## Open Questions

- [ ] Evaluate Neo4j vs Amazon Neptune for graph database
- [ ] Assess pgvector vs dedicated Pinecone for embeddings
- [ ] Determine self-hosted vs managed database strategy
- [ ] Plan for mobile app architecture (React Native vs native)

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Complete architecture specification |
