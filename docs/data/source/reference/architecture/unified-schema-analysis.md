---
id: reference-unified-schema-analysis
title: "Unified Schema Analysis"
category: reference
parent: _index.md
last_updated: 2026-01-23
status: active
migrated_from: operations/schemas/unified-schema-analysis.md
tags: [schema, unified-model, integration, data-model, entity-relationship, cross-database]
---

**Parent:** [Architecture Documentation](./_index.md)

# Unified Schema Analysis

**Document ID:** UNIFIED-SCHEMA-ANALYSIS
**Status:** Final
**Owner:** Architecture Team
**Last Updated:** January 2026
**Version:** 1.0

---

## Executive Summary

This document synthesizes schema patterns across 250+ databases to propose a unified data model for the knowledge base. The analysis identifies five core entity types (Gene, Variant, Compound, Disease, Pathway), maps 15+ identifier systems, and proposes a normalized schema with provenance tracking.

### Key Findings

1. **Hub Identifiers**: UniProt (proteins), NCBI Gene ID (genes), InChIKey (compounds), MONDO (diseases), and Reactome (pathways) serve as optimal canonical identifiers
2. **Common Patterns**: All databases follow similar entity-relationship model with Gene -> Variant -> Disease and Compound -> Target -> Pathway chains
3. **Evidence Complexity**: Evidence quality varies dramatically requiring robust provenance tracking
4. **Integration Challenges**: Identifier fragmentation, naming inconsistencies, and schema heterogeneity require careful normalization

---

## Part 1: Common Entity Type Analysis

### 1.1 Gene Entity Comparison

| Database | Gene ID Type | Gene Count | Key Attributes |
|----------|--------------|------------|----------------|
| dbSNP | NCBI Gene ID | ~30K human | Symbol, location, variants |
| PharmGKB | PA ID + HGNC | 1,761 | Haplotypes, PGx annotations |
| Reactome | UniProt ID | 11,020 | Pathway membership |
| DisGeNET | NCBI Gene ID | 21,671 | Gene-disease associations |
| Gene Ontology | UniProt/NCBI | 20K+ human | Functions, processes |

**Hub Recommendation: NCBI Gene ID (Entrez) with UniProt for protein-level**

### 1.2 Variant Entity Comparison

| Database | Variant ID | Variant Count | Key Attributes |
|----------|------------|---------------|----------------|
| dbSNP | rsID | 1B+ | SPDI format, frequencies |
| ClinVar | VCV/RCV/SCV | 2M+ | Clinical significance |
| gnomAD | chrom-pos-ref-alt | 807M+ | Population frequencies |
| GWAS Catalog | rsID + trait | 400K+ | P-values, effect sizes |
| PharmGKB | rsID + haplotype | 10K+ | Drug response |

**Hub Recommendation: rsID for SNVs, SPDI for normalization**

### 1.3 Compound Entity Comparison

| Database | Compound ID | Count | Key Attributes |
|----------|-------------|-------|----------------|
| ChEMBL | CHEMBL ID | 2.4M | Bioactivity, assays |
| COCONUT | COCONUT ID | 695K | NP-likeness |
| FooDB | FooDB ID | 26K+ | Food sources |
| PubChem | CID | 116M+ | Structures |
| DrugBank | DB ID | 15K | Drug properties |

**Hub Recommendation: InChIKey (structure-based, universal)**

### 1.4 Disease Entity Comparison

| Database | Disease ID | Count | Key Attributes |
|----------|------------|-------|----------------|
| MONDO | MONDO ID | 24K+ | Unified ontology |
| HPO | HP ID | 18K | Phenotypes |
| Orphanet | ORPHA ID | 10K+ | Rare diseases |
| DisGeNET | UMLS CUI | 30K+ | Gene-disease scores |
| OMIM | MIM ID | 27K+ | Mendelian genetics |

**Hub Recommendation: MONDO (best cross-reference coverage)**

---

## Part 2: Identifier Cross-Reference Matrix

### Master Cross-Reference Sources

- UniProt ID Mapping (286 databases)
- NCBI ELink (38 Entrez databases)
- PubChem ID Exchange
- OLS (Ontology Lookup Service)
- Wikidata SPARQL

---

## Part 3: Relationship Pattern Analysis

### Core Relationship Chains

**Chain 1: Genetic Causation**
Gene -> Variant -> Disease -> Phenotype

**Chain 2: Pharmacological Action**
Compound -> Target -> Pathway -> Disease

**Chain 3: Natural Product Discovery**
Organism -> Compound -> Target -> Condition

---

## Part 4: Proposed Unified Schema

### Entity-Relationship Diagram

\`\`\`
+------------+       +------------+       +------------+
|    Gene    |-------|  Variant   |-------|  Disease   |
+------------+       +------------+       +------------+
      |                    |                    |
      v                    v                    v
+------------+       +------------+       +------------+
|  Pathway   |       |  Evidence  |       | Phenotype  |
+------------+       +------------+       +------------+
      |                    ^                    ^
      v                    |                    |
+------------+       +------------+       +------------+
|   Target   |<------|  Compound  |-------|  Organism  |
+------------+       +------------+       +------------+
\`\`\`

---

## Part 5: Integration Architecture

### Storage Recommendations

| Collection | Estimated Records | Storage (indexed) |
|------------|-------------------|-------------------|
| Gene | ~25K | 1 GB |
| Variant | ~10M (clinical) | 100 GB |
| Compound | ~1M | 20 GB |
| Disease | ~50K | 1 GB |
| Pathway | ~10K | 500 MB |
| Evidence | ~100M | 400 GB |
| Associations | ~500M | 200 GB |
| **Total** | - | **~720 GB** |

---

## Glossary

| Term | Definition |
|------|------------|
| hub_identifier | Canonical identifier for cross-database mapping |
| entity | Core data object type in unified schema |
| association_table | Junction table linking entities with evidence |
| provenance | Origin and evidence trail for data assertions |
| SPDI | Sequence-Position-Deletion-Insertion notation |
| InChIKey | 27-character hash of structure identifier |

---

*Full content preserved from original source at operations/schemas/unified-schema-analysis.md*
