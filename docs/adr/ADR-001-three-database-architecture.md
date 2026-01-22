# ADR-001: Three-Database Architecture

## Status
Accepted

## Context

The Gene platform requires a data architecture that supports:
1. **Biomedical knowledge graph** - Genes, SNPs, pathways, herbs, conditions with complex relationships
2. **Vector similarity search** - Find similar genes/SNPs, AI agent memory
3. **User data storage** - HIPAA-compliant storage for genetic data, lab results, PII

Previous architecture attempted to use RuVector Graph for the knowledge graph, but research revealed:
- RuVector is designed for **AI agent memory**, not general-purpose graph databases
- Knowledge graphs benefit from SPARQL and native ontology support
- Each database should be used for its intended purpose

## Decision

Adopt a **three-database architecture** where each database handles its intended use case:

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
│ - SNPs          │     │ Similarity search│    │ Lab results     │
│ - Pathways      │     │ AI agent memory │     │ Symptoms        │
│ - Herbs         │     │ SONA learning   │     │ Treatments      │
│ - Conditions    │     │ Trajectories    │     │                 │
└─────────────────┘     └─────────────────┘     └─────────────────┘
   Apache 2.0              MIT License           PostgreSQL
```

### Database Responsibilities

| Database | Purpose | Data Types |
|----------|---------|------------|
| **QLever** | Structured knowledge graph | Genes, SNPs, pathways, herbs, conditions, relationships, ontologies |
| **RuVector** | Vector search + AI memory | Gene embeddings, SNP embeddings, agent patterns, trajectories, SONA learning |
| **PostgreSQL** | User data (HIPAA) | Profiles, genetic uploads, lab results, PII |

## Consequences

### Positive
- Each database used for its intended purpose
- QLever provides native SPARQL and ontology support (Gene Ontology, Disease Ontology)
- RuVector excels at vector search (150x-12,500x speedup) and AI memory
- PostgreSQL provides mature HIPAA-compliant user data storage
- All databases are FREE (Apache 2.0, MIT, PostgreSQL licenses)
- Clear separation of concerns simplifies debugging and scaling

### Negative
- Cross-database queries require application-level joins
- Three systems to deploy, monitor, and maintain
- Learning curve for SPARQL (QLever)

### Neutral
- Migration from RuVector Graph to QLever for knowledge graph
- Need to implement cross-database query layer in API

## Options Considered

### Option 1: Single Database (PostgreSQL + pgvector)
- **Pros**: Simple, single system
- **Cons**: Poor graph traversal, limited ontology support, vector search not optimized

### Option 2: RuVector for Everything
- **Pros**: Unified stack
- **Cons**: RuVector not designed for structured knowledge graphs; using it beyond intended purpose

### Option 3: Three-Database Architecture (Chosen)
- **Pros**: Each database optimized for its use case, native ontology support, scalable
- **Cons**: Cross-database complexity

### Option 4: Neo4j + Pinecone + PostgreSQL
- **Pros**: Industry standard
- **Cons**: Neo4j Enterprise is NOT FREE ($$$), vendor lock-in

## Related Decisions
- [ADR-002](./ADR-002-knowledge-graph-engine.md): Knowledge Graph Engine Selection
- [ADR-003](./ADR-003-embedding-model-selection.md): Embedding Model Selection
- [ADR-005](./ADR-005-claude-flow-memory.md): Claude Flow Memory Configuration

## References
- [QLever GitHub](https://github.com/ad-freiburg/qlever) - SPARQL engine
- [RuVector Documentation](https://github.com/ruvnet/ruvector) - Vector database
- Research analysis: `docs/research/graphs/`
