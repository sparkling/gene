# Three-Database Architecture Plan

**Document ID:** 46-THREE-DATABASE-ARCHITECTURE
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Three-database architecture using each database for its intended purpose: **QLever** (SPARQL/RDF) for the biomedical knowledge graph, **RuVector** for vector search and AI agent memory, and **PostgreSQL** for HIPAA-compliant user data. This replaces the previous RuVector Graph approach after research revealed RuVector is designed for AI agent memory, not general-purpose graphs. All databases are FREE (Apache 2.0, MIT, PostgreSQL licenses).

---

## Key Decisions

| Decision | Choice | Rationale | Date | ADR |
|----------|--------|-----------|------|-----|
| Database architecture | Three-database | Each DB for intended purpose | Jan 2026 | [ADR-001](../adr/ADR-001-three-database-architecture.md) |
| Knowledge graph engine | QLever | SPARQL, ontology support, 20B+ scale | Jan 2026 | [ADR-002](../adr/ADR-002-knowledge-graph-engine.md) |
| Embedding model | all-mpnet-base-v2 | Simple, proven, FREE | Jan 2026 | [ADR-003](../adr/ADR-003-embedding-model-selection.md) |
| Federation strategy | Local preferred | Reliability over freshness | Jan 2026 | [ADR-004](../adr/ADR-004-federation-strategy.md) |
| Claude Flow memory | RuVector + HNSW | 150x faster pattern search | Jan 2026 | [ADR-005](../adr/ADR-005-claude-flow-memory.md) |

---

## Architecture Overview

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
│                 │     │ Patterns        │     │                 │
│ Local preferred │     │ HNSW indexing   │     │ HIPAA compliant │
│ (federation for │     │                 │     │ Field encryption│
│ large sources)  │     │                 │     │                 │
└─────────────────┘     └─────────────────┘     └─────────────────┘
   Apache 2.0              MIT License           PostgreSQL
```

---

## Database Responsibilities

| Database | Purpose | Data Types |
|----------|---------|------------|
| **QLever** | Structured knowledge graph | Genes, SNPs, pathways, herbs, conditions, relationships, ontologies |
| **RuVector** | Vector search + AI memory | Gene embeddings, SNP embeddings, agent patterns, trajectories, SONA learning |
| **PostgreSQL** | User data (HIPAA) | Profiles, genetic uploads, lab results, PII |

---

## Why This Architecture

### QLever Provides
- **Local import preferred** for reliability and speed
- **Disk-backed storage** with OS page cache (not purely in-memory)
- Native ontology support (Gene Ontology, Disease Ontology)
- Optimized multi-hop traversals (gene→pathway→condition→herb)
- Scales to 20B+ triples (Wikidata scale)
- Apache 2.0 license (FREE)
- **Fallback strategies**:
  1. Multiple QLever instances with internal federation
  2. External SPARQL federation for data too large to import

### QLever Resource Requirements (Our Scale)

Our dataset (~1-2M SNPs, ~20K genes, ~2K pathways) ≈ **15-35M triples**

| Metric | Requirement | Notes |
|--------|-------------|-------|
| Index size | ~1-3 GB | Stored on disk (SSD preferred) |
| Query RAM | 4-8 GB | OS caches hot data automatically |
| Indexing RAM | 4-8 GB | One-time build process |
| Indexing time | Seconds-minutes | QLever indexes 5B triples/hour |
| Disk space | 5-10 GB | Index + source RDF files |

**This is tiny by QLever standards** - Wikidata (20B triples) needs only ~20GB query RAM.

### RuVector Provides
- HNSW vector search (150x-12,500x speedup)
- Gene/SNP embeddings for similarity search
- AI agent memory (its intended use case)
- SONA learning system
- Pattern storage with confidence scores
- MIT license (FREE)

### PostgreSQL Provides
- HIPAA-compliant user data storage
- Field-level encryption
- Mature, battle-tested
- Audit logging

---

## Query Flow Examples

### Example 1: "Find genes similar to MTHFR"
```
1. App → RuVector: similarity_search(mthfr_embedding, k=10)
2. RuVector returns: [MTRR, MTR, SHMT1, CBS, ...]
3. App → QLever: SPARQL query for gene details
```

### Example 2: "What pathways involve MTHFR and affect folate?"
```
1. App → QLever: SPARQL query (all data local)
   - gene-pathway relationships
   - protein annotations (imported from UniProt)
2. Returns: Methylation cycle, one-carbon metabolism, ...
```

### Example 3: "Find herbs for user's SNP profile"
```
1. App → PostgreSQL: Get user's SNP list
2. App → RuVector: Get similar SNP patterns from learned data
3. App → QLever: SPARQL query herbs→conditions→pathways→SNPs
4. Combine results with AI-learned recommendations
```

---

## Implementation Phases

### Phase 1: Infrastructure Setup
- [ ] Deploy QLever Docker container (adfreiburg/qlever)
- [ ] Configure RuVector with redb + HNSW
- [ ] Ensure PostgreSQL is running (existing)
- [ ] Create cross-database query layer in API

### Phase 2: QLever Knowledge Graph
- [ ] Download and convert biomedical data to RDF/N-Triples:
  - UniProt human proteins (~20K)
  - dbSNP/ClinVar variants (~1.2M clinically relevant)
  - Reactome/KEGG pathways (~2K)
  - Gene Ontology, Disease Ontology
- [ ] Use rdflib (Python) or Apache Jena for RDF conversion
- [ ] Index genes (~20K), SNPs (~1.2M), pathways (~2K), herbs (~2K)
- [ ] Test ontology queries

### Phase 3: RuVector Setup
- [ ] Initialize RuVector with appropriate storage backend
- [ ] Generate gene embeddings (using all-mpnet-base-v2)
- [ ] Generate SNP embeddings
- [ ] Configure HNSW index parameters (M=32, ef_construction=200, ef_search=100)
- [ ] Setup AI agent memory namespaces

### Phase 3b: Claude Flow Memory Configuration
- [ ] Audit current `.claude/settings.json` hooks configuration
- [ ] Fix memory backend configuration in `claude-flow.config.json`
- [ ] Test memory store/search/retrieve commands
- [ ] Verify SONA learning trajectories persist across sessions
- [ ] Configure namespaces: `patterns`, `trajectories`, `tasks`, `solutions`

### Phase 4: Integration Layer
- [ ] Create unified API that routes queries to appropriate DB
- [ ] Implement cross-database join patterns
- [ ] Build caching layer for frequent queries
- [ ] Setup monitoring for all three databases

---

## Configuration

### QLever Docker
```yaml
services:
  qlever:
    image: adfreiburg/qlever:latest
    container_name: gene-qlever
    ports:
      - "7001:7001"
    volumes:
      - qlever_data:/data
      - ./data/rdf:/input
    environment:
      MEMORY_FOR_QUERIES: 8G
      CACHE_MAX_SIZE: 4G
```

### RuVector Config (Hardware-Optimized)

**Target Hardware:** Ryzen 9 7950X3D + 187GB RAM + NVMe RAID

```typescript
const ruvector = new RuVector({
  storagePath: './data/claude-flow-memory',
  dimensions: 768,
  metric: 'cosine',
  hnswM: 48,           // Higher - RAM abundant
  hnswEfConstruction: 400,  // Higher - NVMe fast
  hnswEfSearch: 200,   // Higher - V-Cache helps
  quantization: 'none', // Full f32 - no need to compress
  cacheSize: 50000     // Large cache
});
```

### Claude Flow Memory (Hardware-Optimized)
```json
{
  "memory": {
    "backend": "hybrid",
    "storagePath": "./data/claude-flow-memory",
    "hnsw": {
      "enabled": true,
      "dimensions": 768,
      "M": 48,
      "efConstruction": 400,
      "efSearch": 200
    },
    "quantization": "none",
    "cacheSize": 50000,
    "maxEntries": 10000000
  },
  "learning": {
    "sona": { "enabled": true },
    "trajectories": { "persist": true },
    "patterns": { "autoStore": true }
  },
  "ruvector": {
    "moe": { "enabled": true, "experts": 8 },
    "flash": { "enabled": true, "blockSize": 128 }
  }
}
```

---

## Embedding Strategy

### Key Principle: Always Embed with Context

Gene symbols and SNP IDs alone don't embed well. Always enrich with context:

```typescript
// ❌ POOR - bare identifier has no semantic meaning
await embed("MTHFR");

// ✅ GOOD - context-enriched embedding
await embed("MTHFR: Methylenetetrahydrofolate reductase, key enzyme in methylation cycle");

// ✅ GOOD - SNP with clinical significance
await embed("rs1801133 (MTHFR C677T): Reduces enzyme activity, impacts folate requirements");
```

### Model Selection: Start Simple

| Phase | Model | Dimensions | Use For |
|-------|-------|------------|---------|
| **Phase 1 (Now)** | all-mpnet-base-v2 | 768 | Everything |
| **Phase 2 (If needed)** | Add PubMedBERT | 768 | Biomedical text |
| **Phase 3 (If needed)** | Add Gene2Vec | 200 | Gene-gene relationships |
| **Phase 4 (If sequences)** | Add ESM-2 | 320-480 | Protein sequences |

### Namespace Strategy in RuVector

| Namespace | Content | Model |
|-----------|---------|-------|
| `genes` | Gene embeddings (symbol + description) | all-mpnet-base-v2 |
| `snps` | SNP embeddings (rsid + clinical) | all-mpnet-base-v2 |
| `patterns` | AI agent learned patterns | all-mpnet-base-v2 |
| `trajectories` | SONA learning trajectories | all-mpnet-base-v2 |

---

## Verification

### Test QLever
```bash
# Query local knowledge graph
curl "http://localhost:7001/api?query=SELECT * WHERE { ?s ?p ?o } LIMIT 10"

# Query gene-pathway relationships
curl "http://localhost:7001/api?query=SELECT ?gene ?pathway WHERE { ?gene :participatesIn ?pathway } LIMIT 10"
```

### Test RuVector
```typescript
// Test vector similarity search
const similar = await ruvector.search(geneEmbedding, { k: 10 });

// Test AI memory
const patterns = await ruvector.memory.search("methylation patterns");
```

### Test Claude Flow Memory
```bash
npx @claude-flow/cli@latest memory store --key "test" --value "hello" --namespace patterns
npx @claude-flow/cli@latest memory search --query "hello" --namespace patterns
npx @claude-flow/cli@latest hooks intelligence stats
```

### Test Embedding Quality
```typescript
// Verify similar genes cluster together
const mthfr_mtr = cosineSimilarity(mthfrEmbed, mtrEmbed);
const mthfr_brca1 = cosineSimilarity(mthfrEmbed, brca1Embed);
assert(mthfr_mtr > mthfr_brca1, 'Methylation genes should cluster together');
```

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| RuVector 0.1.x stability | Gene embeddings are regenerable; not critical user data |
| QLever learning curve (SPARQL) | Team training; start with simple queries |
| Cross-DB complexity | Clear separation of concerns; documented query patterns |
| Data staleness | Periodic re-import from upstream sources (quarterly) |

---

## Backup & Recovery

| Database | Strategy | Recovery |
|----------|----------|----------|
| **QLever** | RDF source files in git/S3 | Re-run indexer (~minutes for our scale) |
| **RuVector** | Embeddings are regenerable | Re-embed from gene/SNP descriptions |
| **PostgreSQL** | Standard pg_dump; PITR | Restore from backup |

**Note**: RuVector embeddings are derived data, not source of truth. If lost, regenerate from QLever gene/SNP descriptions.

---

## Files to Modify

| File | Changes |
|------|---------|
| `docker-compose.yml` | Add QLever service |
| `src/lib/qlever-client.ts` | New - SPARQL query client |
| `src/lib/ruvector-client.ts` | New - Vector search + AI memory client |
| `src/lib/knowledge-graph.ts` | Update - Route to QLever |
| `src/lib/cross-db-queries.ts` | Update - Add QLever integration |
| `claude-flow.config.json` | Fix memory backend config |
| `.claude/settings.json` | Audit/fix hooks configuration |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [44-ARCHITECTURE](./44-architecture.md) | Parent architecture (to be updated) |
| [43-DATA-SOURCES](../../data/source/INDEX.md) | Data source inventory |
| [ADR Index](../adr/README.md) | Architecture decisions |

---

## Open Questions

- [x] **Embedding model**: all-mpnet-base-v2 first
- [x] **RuVector use case**: AI memory + gene/SNP embeddings
- [x] **Federation strategy**: Local import preferred
- [x] **QLever storage**: Disk-backed, 4-8GB RAM sufficient
- [x] **Documentation**: ADRs + McKinsey framework

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial three-database architecture plan |
