# ADR-002: Knowledge Graph Engine Selection

## Status
Accepted

## Context

The Gene platform requires a knowledge graph engine to store and query:
- ~20K genes with annotations
- ~1.2M SNPs with clinical significance
- ~2K pathways (Reactome, KEGG)
- ~2K herbs with traditional medicine mappings
- ~5K conditions (MeSH, ICD-10)
- Complex relationships between all entities
- Ontology support (Gene Ontology, Disease Ontology)

Requirements:
1. Multi-hop graph traversals (gene→pathway→condition→herb)
2. SPARQL or compatible query language
3. Native ontology support
4. Scales to millions of entities
5. FREE and open source
6. Federation capability for external data sources

## Decision

Use **QLever** as the knowledge graph engine.

### Why QLever

| Feature | QLever | RuVector Graph | Apache AGE |
|---------|--------|----------------|------------|
| Query Language | SPARQL 1.1 | Cypher-like | Cypher + SQL |
| Ontology Support | Native | None | Limited |
| Scale | 20B+ triples | Unknown | Millions |
| Federation | Native SPARQL | None | None |
| License | Apache 2.0 | MIT | Apache 2.0 |
| RAM (our scale) | 4-8 GB | Unknown | 8-16 GB |
| Storage | Disk-backed | Disk-backed | Disk-backed |

### Resource Requirements (Our Scale)

Dataset: ~15-35M triples (genes, SNPs, pathways, herbs, conditions)

| Metric | Requirement | Notes |
|--------|-------------|-------|
| Index size | ~1-3 GB | Stored on disk (SSD preferred) |
| Query RAM | 4-8 GB | OS caches hot data automatically |
| Indexing RAM | 4-8 GB | One-time build process |
| Indexing time | Seconds-minutes | QLever indexes 5B triples/hour |
| Disk space | 5-10 GB | Index + source RDF files |

**This is tiny by QLever standards** - Wikidata (20B triples) needs only ~20GB query RAM.

## Consequences

### Positive
- Native SPARQL 1.1 support with federation
- Handles Gene Ontology, Disease Ontology natively
- Proven at Wikidata scale (20B+ triples)
- Disk-backed with OS page cache (not purely in-memory)
- Apache 2.0 license (FREE)
- Optimized for multi-hop traversals

### Negative
- SPARQL learning curve for team
- Requires RDF/N-Triples data conversion
- Less common than Neo4j in industry

### Neutral
- Need to convert biomedical data to RDF format
- Docker deployment (adfreiburg/qlever)

## Options Considered

### Option 1: RuVector Graph (@ruvector/graph-node)
- **Pros**: Unified with vector database, Cypher-like queries
- **Cons**: No ontology support, not designed for structured knowledge graphs, limited documentation

### Option 2: Apache AGE (PostgreSQL extension)
- **Pros**: PostgreSQL ecosystem, SQL + Cypher
- **Cons**: No native ontology support, no federation, less optimized for graph traversals

### Option 3: QLever (Chosen)
- **Pros**: SPARQL, ontology support, federation, proven scale, FREE
- **Cons**: SPARQL learning curve

### Option 4: Neo4j
- **Pros**: Industry standard, excellent tooling
- **Cons**: Enterprise features require paid license, NOT FREE for commercial use

### Option 5: TypeDB
- **Pros**: Powerful type system
- **Cons**: Smaller community, less biomedical ecosystem support

## Related Decisions
- [ADR-001](./ADR-001-three-database-architecture.md): Three-Database Architecture
- [ADR-004](./ADR-004-federation-strategy.md): SPARQL Federation Strategy

## References
- [QLever GitHub](https://github.com/ad-freiburg/qlever)
- [QLever Performance Evaluation](https://github.com/ad-freiburg/qlever/wiki/QLever-performance-evaluation-and-comparison-to-other-SPARQL-engines)
- [SPARQL 1.1 Specification](https://www.w3.org/TR/sparql11-query/)
- Research: `docs/40-product/data-sources/research/hive-ruvector-analysis/07-QLEVER-ANALYSIS.md`
