---
id: pharmaceuticals
title: "Pharmaceuticals"
type: subcategory
parent: ../_index.md
last_updated: 2026-01-23
status: active
tags: [drugs, pharmaceuticals, approved, clinical, fda]
---

# Pharmaceuticals

**Parent:** [Compounds & Molecules](../_index.md)

## Overview

Pharmaceutical databases catalog approved drugs, drug candidates, and their associated bioactivity, targets, and clinical information. These resources support drug repurposing, pharmacovigilance, and understanding drug mechanisms of action.

Key resources include DrugBank (comprehensive drug database), ChEMBL (bioactivity data), RxNorm (clinical drug nomenclature), DailyMed (FDA label information), and Orange Book (approved drug products). Together they provide complete drug information from discovery to market.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [DrugBank](./drugbank/_index.md) | 1 | Comprehensive drug database |
| [ChEMBL](./chembl/_index.md) | 1 | Bioactivity database |
| [RxNorm](./rxnorm/_index.md) | 1 | Clinical drug nomenclature |
| [DailyMed](./dailymed/_index.md) | 2 | FDA drug labeling |
| [Orange Book](./orange.book/_index.md) | 2 | FDA approved products |

## Integration Notes

DrugBank provides the most complete drug annotations including targets and pathways. ChEMBL offers extensive bioactivity assay data. RxNorm standardizes clinical drug names for EHR integration. DailyMed provides official FDA labeling text. Use DrugBank IDs for cross-referencing.

---

## Database Selection Guide

### TL;DR

Pharmaceutical and pharmacogenomics databases provide the foundation for SNP-drug relationship mapping, clinical dosing recommendations, and drug-target interactions. PharmGKB/CPIC serve as the gold standard for pharmacogenomics with clinically validated dosing guidelines, while DrugBank and ChEMBL provide comprehensive drug information and bioactivity data. Integration of 22+ databases enables complete coverage from variant annotation through regulatory compliance.

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Primary PGx source | PharmGKB + CPIC + DPWG | Gold standard with clinical validation |
| Drug information | DrugBank (academic) + ChEMBL | Comprehensive targets + open bioactivity |
| Drug-gene aggregator | DGIdb + Open Targets | Free aggregation from 40+ sources |
| Regulatory data | OpenFDA + DailyMed | Authoritative FDA labels and adverse events |
| Binding affinity | BindingDB + ChEMBL | Quantitative IC50/Ki/Kd data |
| Allele nomenclature | PharmVar | Official star allele definitions |

### Comparative Statistics

| Database | Drugs | Targets | Interactions | SNP-Drug | License |
|----------|-------|---------|--------------|----------|---------|
| PharmGKB | 1,000+ | 700+ genes | 20,000+ annotations | Yes | CC BY-SA 4.0 |
| CPIC | 164 | 34 genes | - | Yes | CC0 |
| DrugBank | 15,000+ | 5,000+ | 2.5M DDI | Yes | Academic/Commercial |
| ChEMBL | 2.8M compounds | 17,803 | 21M bioactivities | No | Open |
| DGIdb | 20,000+ | 10,000+ | 70,000+ | No | Open |
| BindingDB | 1.3M | Thousands | 2.9M affinity | No | Open |

### Tier Priorities

| Tier | Databases | Priority Reason |
|------|-----------|-----------------|
| **1 (Critical)** | PharmGKB, CPIC, DrugBank, OpenFDA, DailyMed, UniProt | PGx gold standard; regulatory compliance; target mapping |
| **2 (High Value)** | ChEMBL, DGIdb, Open Targets, BindingDB, GtoPdb, CIViC | Bioactivity; drug-gene interactions; oncology |
| **3 (Enrichment)** | PubChem, TTD, SuperDRUG2, KEGG Drug | Chemistry reference; druggability; pathway links |

### Licensing Summary

| License Type | Databases | Commercial Use |
|--------------|-----------|----------------|
| CC0 (Public Domain) | CPIC, CIViC, Open Targets | Yes |
| Public Domain | OpenFDA, DailyMed, PubChem | Yes |
| CC BY-SA 4.0 | PharmGKB | Yes (with attribution) |
| Open Access | ChEMBL, DGIdb, BindingDB, GtoPdb | Yes |
| Academic Only | DrugBank (full data), KEGG Drug | License required |
