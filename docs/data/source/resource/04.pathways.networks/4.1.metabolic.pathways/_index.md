---
id: metabolic.pathways
title: "Metabolic Pathways"
type: subcategory
parent: ../_index.md
last_updated: 2026-01-23
status: active
tags: [metabolism, biochemistry, kegg, reactome, enzymes]
---

# Metabolic Pathways

**Parent:** [Pathways & Networks](../_index.md)

## Overview

Metabolic pathway databases catalog the biochemical reactions and enzymes involved in cellular metabolism. These resources provide pathway maps, reaction stoichiometry, and enzyme annotations essential for metabolic modeling and analysis.

Key resources include KEGG (comprehensive pathway database), Reactome (curated pathway knowledgebase), and WikiPathways (community-curated pathways). Together they provide multiple perspectives on metabolic organization.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [KEGG](./kegg/_index.md) | 1 | Kyoto Encyclopedia of Genes and Genomes |
| [Reactome](./reactome/_index.md) | 1 | Curated pathway knowledgebase |
| [WikiPathways](./wikipathways/_index.md) | 2 | Community pathway database |

## Integration Notes

KEGG provides the most comprehensive pathway coverage with KEGG IDs. Reactome offers detailed reaction mechanisms with evidence. WikiPathways provides community curation and GPML format. Use pathway enrichment tools to map gene lists to pathways. Cross-reference enzymes with UniProt.

---

## Database Selection Guide

### TL;DR

Primary pathway databases provide curated biological pathway data essential for understanding gene-compound-disease relationships. Reactome (CC BY 4.0, 2,712 human pathways) and WikiPathways (CC0, 3,100+ pathways) are recommended as open-access priorities, with KEGG, MetaCyc, SMPDB, PharmGKB Pathways, PANTHER, Pathway Commons, and NDEx providing specialized coverage for metabolism, pharmacogenomics, and network analysis.

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Primary pathway source | Reactome | CC BY 4.0 license allows commercial use; highest curation quality; Neo4j + REST API |
| Secondary pathway source | WikiPathways | CC0 public domain; community contributions capture emerging biology |
| Drug metabolism pathways | SMPDB | Best HMDB/DrugBank integration; comprehensive PK pathways |
| Pharmacogenomics pathways | PharmGKB Pathways | Clinical variant-drug relationships; CPIC guideline integration |
| KEGG approach | Reference only | Bulk FTP requires subscription; use Reactome/WikiPathways for open alternatives |
| Integration layer | Pathway Commons | Pre-integrated data from multiple sources; BioPAX standard format |
| Network repository | NDEx | Access specialized networks; NCI-PID archive preserved |
| Primary compound identifier | ChEBI | Used by Reactome, WikiPathways, Pathway Commons; open ontology |

### Comparative Statistics

| Database | Pathways | Genes | Key Feature | License |
|----------|----------|-------|-------------|---------|
| Reactome | 2,712 | ~12K | High curation quality | CC BY 4.0 |
| WikiPathways | 3,100+ | ~8K | Community-curated | CC0 |
| KEGG | ~500 human | ~8K | Reference standard | Subscription |
| SMPDB | ~50K | - | Drug metabolism | Open |
| Pathway Commons | Aggregated | - | Multi-source integration | Mixed |
| PANTHER | ~177 | ~2K | Gene family context | Open |

### Tier Priorities

| Tier | Databases | Priority Reason |
|------|-----------|-----------------|
| **1 (Critical)** | Reactome, WikiPathways, SMPDB | Open licenses; comprehensive coverage; drug metabolism |
| **2 (High Value)** | Pathway Commons, PharmGKB Pathways, PANTHER | Integration layer; pharmacogenomics; gene families |
| **3 (Reference)** | KEGG, MetaCyc, NDEx | Reference standard; specialized metabolism; network sharing |

### Licensing Summary

| License Type | Databases | Commercial Use |
|--------------|-----------|----------------|
| CC BY 4.0 | Reactome | Yes (with attribution) |
| CC0 (Public Domain) | WikiPathways | Yes |
| Open Access | SMPDB, PANTHER, PathBank | Yes |
| Mixed Sources | Pathway Commons | Varies by source |
| Subscription | KEGG (bulk FTP) | Academic free; commercial requires license |
