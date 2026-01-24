---
id: reference-three-worlds-schema
title: RuVector THREE WORLDS Schema Design
category: reference
parent: _index.md
last_updated: 2026-01-23
status: active
migrated_from: operations/schemas/ruvector-three-worlds-schema.md
tags: [schema, database, architecture, ruvector]
---

# RuVector THREE WORLDS Schema Design

**Document ID:** RUVECTOR-THREE-WORLDS-SCHEMA
**Status:** Architecture Decision
**Owner:** System Architecture
**Date:** January 21, 2026
**Version:** 1.0
**Parent:** [Architecture Documentation](./_index.md)

---

## Executive Summary

This document defines the optimal schema for storing THREE WORLDS data (Genetics, Traditional Medicine, Nutritional Science) using RuVector's native capabilities.

### Scale Requirements

| Domain | Entity Type | Estimated Count | Storage Strategy |
|--------|-------------|-----------------|------------------|
| **WORLD 1: Genetics** | SNPs/Variants | ~1.2B | Domain-sharded, quantized vectors |
| **WORLD 1: Genetics** | Genes | ~25K | Hot cache, full precision |
| **WORLD 2: Traditional Medicine** | Compounds (NP) | ~800K | HNSW-indexed, structure vectors |
| **WORLD 2: Traditional Medicine** | Formulas | ~60K | Hyperedge patterns |
| **WORLD 3: Nutrition** | Foods | ~380K | Multi-vector (composition + embedding) |
| **Cross-Domain** | Pathways | ~10K | Graph topology, hyperedge-native |
| **Cross-Domain** | Diseases | ~80K | Ontology hierarchy in graph |

---

## Part 1: RuVector Capability Mapping

### Core Components Used

| RuVector Component | THREE WORLDS Use Case | Why This Component |
|--------------------|----------------------|-------------------|
| VectorDB | Gene/compound embeddings, semantic search | 150x faster HNSW, metadata filtering |
| CodeGraph | Pathways, gene networks, compound-target | Native hyperedges for N-ary relations |
| RuvectorCluster | Cross-node queries, domain sharding | Raft consensus, auto-rebalancing |
| SONA Engine | Learning from query patterns | Micro-LoRA adaptation (<0.1ms) |

---

## Part 2: Node Type Definitions

### WORLD 1: Genetics Node Types

**Gene Node** - Primary identifiers: NCBI Gene ID, HGNC, Ensembl, UniProt
**Variant Node** - Primary identifiers: rsID, ClinVar VCV, gnomAD ID
**Haplotype Node** - Star alleles for pharmacogenomics

### WORLD 2: Traditional Medicine Node Types

**Compound Node** - Primary identifier: InChIKey
**Organism Node** - Source species with NCBI Taxonomy ID
**Formula Node** - Traditional preparations from KampoDB, BATMAN-TCM, IMPPAT

### WORLD 3: Nutritional Science Node Types

**Food Node** - Sources: USDA, FooDB, OpenFoodFacts
**Nutrient Node** - USDA nutrient IDs with RDA values

### Cross-Domain Node Types

**Disease Node** - Primary identifier: MONDO
**Pathway Node** - Primary identifier: Reactome
**Phenotype Node** - HPO terms

---

## Part 3: Edge Type Definitions

### Full Edge Type Catalog

| Edge Type | From Node | To Node | Properties |
|-----------|-----------|---------|------------|
| HAS_VARIANT | Gene | Variant | consequence, impact, exon |
| ASSOCIATED_WITH | Variant | Disease | significance, inheritance |
| INHIBITS | Compound | Gene/Protein | activity_type, activity_value |
| PRODUCES | Organism | Compound | plant_part, concentration |
| PARTICIPATES_IN | Gene | Pathway | role, evidence |
| HAS_PHENOTYPE | Disease | Phenotype | frequency, onset |
| METABOLIZED_BY | Compound | Gene | pharmgkb_level, guideline |

---

## Part 4: Hyperedge Patterns (N-ary Relationships)

### Why Hyperedges?

Traditional binary edges cannot capture relationships involving 3+ entities:
- Formula composition (Formula + ingredients + dosages)
- Pathway reactions (Substrates + enzymes + products)
- Pharmacogenomic interactions (Drug + gene + variant + phenotype)

### Hyperedge Types

- FORMULA_COMPOSITION
- PATHWAY_REACTION
- PHARMACOGENOMIC_INTERACTION
- NUTRIENT_PROFILE
- DRUG_DRUG_GENE_INTERACTION

---

## Part 5: Vector Embedding Integration Points

### Embedding Storage Strategy

| Entity Type | VectorDB Key Format | Dimensions | Use Case |
|-------------|---------------------|------------|----------|
| Gene semantic | GENE_SEM:{ncbi_id} | 384 | Functional similarity |
| Variant population | VAR_POP:{rsid} | 8 | Population similarity |
| Compound structure | COMP_FP:{inchikey} | 2048 | Structure similarity |
| Compound semantic | COMP_SEM:{inchikey} | 384 | Bioactivity similarity |
| Food nutrient | FOOD_NUT:{id} | 50 | Nutrient similarity |
| Disease semantic | DIS_SEM:{mondo_id} | 384 | Disease similarity |

---

## Part 6: Sharding Strategy

### Domain-Based Sharding (16 Shards)

- Shards 0-3: GENETICS_HOT (~10M clinically relevant variants)
- Shards 4-7: GENETICS_COLD (~1.19B rare variants, quantized)
- Shards 8-9: GENES (~25K genes)
- Shards 10-11: COMPOUNDS (~800K compounds)
- Shards 12-13: FOODS (~380K foods)
- Shards 14-15: ONTOLOGIES (diseases, phenotypes, pathways)

---

## Part 7: Performance Considerations

### Expected Query Performance

| Query Type | Expected Latency |
|------------|------------------|
| Single node lookup | <5ms |
| Vector similarity (k=20) | <50ms |
| Graph traversal (2 hops) | <100ms |
| Hyperedge lookup | <20ms |
| Cross-shard query (10 shards) | <500ms |

### Storage Estimates

| Data Tier | Storage Estimate |
|-----------|------------------|
| Hot vectors | ~20 GB |
| Cold vectors (1.2B variants) | ~150 GB |
| Graph storage | ~50 GB |
| Hyperedge storage | ~20 GB |
| Indices | ~100 GB |
| **Total per node** | **~340 GB** |

---

## Glossary

| Term | Definition |
|------|------------|
| node | Graph database entity |
| edge | Binary relationship between nodes |
| hyperedge | N-ary relationship connecting 3+ nodes |
| shard | Partition of distributed database |
| HNSW | Hierarchical Navigable Small World |
| SONA | Self-Optimizing Neural Architecture |

---

*Full content preserved from original source at operations/schemas/ruvector-three-worlds-schema.md*
