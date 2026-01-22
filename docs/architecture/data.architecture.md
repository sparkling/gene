# Data Architecture Decisions Summary (ADR 001-004)

This document summarizes the core data architecture decisions for the Gene platform.

## Overview

The Gene platform uses a three-database architecture, each optimized for its specific purpose:

| Database | Purpose | Key Technology |
|----------|---------|----------------|
| QLever | Knowledge graph (SPARQL/RDF) | Genes, SNPs, pathways, herbs, conditions |
| RuVector | Vector search + AI memory | Embeddings, similarity search, SONA learning |
| PostgreSQL | User data (HIPAA) | Profiles, genetic uploads, lab results |

---

## ADR-001: Three-Database Architecture

**Status**: Accepted

**Decision**: Adopt a three-database architecture where each database handles its intended use case.

**Key Points**:
- QLever for structured knowledge graph with native SPARQL and ontology support
- RuVector for vector similarity search (150x-12,500x speedup) and AI agent memory
- PostgreSQL for HIPAA-compliant user data storage
- All databases are FREE (Apache 2.0, MIT, PostgreSQL licenses)

**Trade-offs**:
- (+) Each database optimized for its use case
- (+) Native ontology support (Gene Ontology, Disease Ontology)
- (-) Cross-database queries require application-level joins
- (-) Three systems to deploy and maintain

---

## ADR-002: Knowledge Graph Engine Selection

**Status**: Accepted

**Decision**: Use QLever as the knowledge graph engine.

**Dataset Scale**: ~15-35M triples (genes, SNPs, pathways, herbs, conditions)

**Resource Requirements**:
| Metric | Requirement |
|--------|-------------|
| Index size | ~1-3 GB |
| Query RAM | 4-8 GB |
| Indexing RAM | 4-8 GB |
| Indexing time | Seconds-minutes |

**Why QLever**:
- Native SPARQL 1.1 support with federation
- Handles Gene Ontology, Disease Ontology natively
- Proven at Wikidata scale (20B+ triples)
- Disk-backed with OS page cache
- Apache 2.0 license (FREE)

**Alternatives Rejected**:
- RuVector Graph: No ontology support, not designed for structured knowledge graphs
- Apache AGE: No native ontology support, no federation
- Neo4j: Enterprise features require paid license

---

## ADR-003: Embedding Model Selection

**Status**: Accepted

**Decision**: Use all-mpnet-base-v2 (768 dimensions) as the primary embedding model. Start simple, add specialized models only if needed.

**Phased Approach**:
| Phase | Model | Use For |
|-------|-------|---------|
| Phase 1 (Now) | all-mpnet-base-v2 | Everything - solid for biomedical text |
| Phase 2 (If needed) | PubMedBERT | Better gene/SNP description similarity |
| Phase 3 (If needed) | Gene2Vec | Gene-gene functional relationships |
| Phase 4 (If sequences) | ESM-2 | Protein sequence similarity |

**Key Principle**: Always embed with context, not bare identifiers.

```typescript
// BAD - bare identifier
await embed("MTHFR");

// GOOD - context-enriched
await embed("MTHFR: Methylenetetrahydrofolate reductase, key enzyme in methylation cycle");
```

---

## ADR-004: SPARQL Federation Strategy

**Status**: Accepted

**Decision**: Local import preferred, with external SPARQL federation only for data too large to import.

**Strategy (In Order of Preference)**:
1. **Local Import** (Preferred) - Import all needed data into local QLever
2. **Multiple QLever Instances** - Internal federation if data exceeds single instance
3. **External SPARQL Federation** - Only for data too large to import (e.g., full UniProt 500M+ triples)

**Data Import Plan**:
| Source | Triples | Strategy |
|--------|---------|----------|
| Human proteins (UniProt subset) | ~500K | Local import |
| dbSNP/ClinVar variants | ~5M | Local import |
| Reactome/KEGG pathways | ~500K | Local import |
| Gene Ontology | ~2M | Local import |
| Disease Ontology | ~200K | Local import |
| Herb databases | ~1M | Local import |
| Full UniProt (all species) | ~500M | External federation |

**Total local import**: ~10M triples - well within QLever's capacity.

---

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                    Gene Platform API                             │
└─────────────────────────────────────────────────────────────────┘
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│     QLever      │     │    RuVector     │     │   PostgreSQL    │
│   (SPARQL/RDF)  │     │  (Vectors/AI)   │     │   (User Data)   │
├─────────────────┤     ├─────────────────┤     ├─────────────────┤
│ Knowledge Graph │     │ Gene embeddings │     │ User profiles   │
│ - Genes         │     │ SNP embeddings  │     │ Genetic data    │
│ - SNPs          │     │ Similarity      │     │ Lab results     │
│ - Pathways      │     │ AI agent memory │     │ Symptoms        │
│ - Herbs         │     │ SONA learning   │     │ Treatments      │
│ - Conditions    │     │ Trajectories    │     │                 │
└─────────────────┘     └─────────────────┘     └─────────────────┘
```

---

## Related Documents

- [ADR-001](../adr/adr-001-three-database-architecture.md): Full three-database architecture decision
- [ADR-002](../adr/adr-002-knowledge-graph-engine.md): Knowledge graph engine selection details
- [ADR-003](../adr/adr-003-embedding-model-selection.md): Embedding model selection rationale
- [ADR-004](../adr/adr-004-federation-strategy.md): Federation strategy details
