---
id: reference-architecture
title: Architecture Documentation
description: System architecture, schema designs, and data model specifications
category: reference
parent: ../_index.md
weight: 20
last_updated: 2026-01-23
status: active
tags: [architecture, schema, data-model, ruvector, design]
---

# Architecture Documentation

System architecture documentation including schema designs, data models, and integration patterns for the knowledge base.

## Contents

### [Three Worlds Schema](./three-worlds-schema.md)
RuVector database schema design for the THREE WORLDS data model:
- **World 1: Genetics** - Genes, variants, haplotypes (1.2B+ variants)
- **World 2: Traditional Medicine** - Compounds, organisms, formulas (800K+ compounds)
- **World 3: Nutrition** - Foods, nutrients, phytochemicals (380K+ foods)
- Cross-domain entities: Pathways, diseases, phenotypes
- Hyperedge patterns for N-ary relationships
- Vector embedding integration points
- Sharding strategy for distributed storage

### [Unified Schema Analysis](./unified-schema-analysis.md)
Cross-database entity mapping and relationship patterns:
- Common entity type analysis (Gene, Variant, Compound, Disease, Pathway)
- Identifier cross-reference matrix
- Relationship pattern analysis
- Proposed unified data model
- Evidence/provenance tracking schema
- Integration architecture recommendations

### [Sample Data](./sample-data.md)
Example API responses and data structures:
- Reactome pathway and entity examples
- WikiPathways search results
- DisGeNET gene-disease associations
- Cross-database relationship mapping
- Database statistics summary

---

## Architecture Overview

```
+------------------+     +------------------+     +------------------+
|   Source Layer   |     |  Transform Layer |     |   Target Layer   |
+------------------+     +------------------+     +------------------+
                              |
  250+ Databases  ------>  Unified Schema  ------>  RuVector
                              |
                    ID Resolution + Evidence Tracking
```

## Key Design Principles

1. **Hub Identifiers**: Use canonical IDs (NCBI Gene, InChIKey, MONDO, Reactome) for cross-database mapping
2. **Evidence Tracking**: Full provenance for all assertions with source, method, and confidence
3. **Hyperedge Support**: Native N-ary relationships for formulas, reactions, interactions
4. **Vector Integration**: Semantic embeddings for similarity search alongside graph queries
5. **Hot/Cold Tiering**: Frequently accessed data at full precision, rare data compressed

---

## Related Documentation

- [Formats](../formats/_index.md) - Data format specifications
- [Schemas Index](../../operations/schemas/_index.md) - Database-specific schemas
- [Integration Guide](../../operations/integration/_index.md) - Integration procedures
