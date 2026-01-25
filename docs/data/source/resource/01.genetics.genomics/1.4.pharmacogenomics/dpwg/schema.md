---
id: schema-dpwg
title: "DPWG Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, database, pharmacogenomics, guidelines, dutch]
---

# DPWG Schema Documentation

**Document ID:** SCHEMA-DPWG
**Version:** 1.0
**Source Version:** Current (via PharmGKB)

---

## TL;DR

DPWG provides Dutch pharmacogenetics guidelines integrated into the G-Standaard national drug database. Guidelines are accessible internationally via PharmGKB with structured gene-drug recommendations and urgency classifications.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Gene-Drug Pairs | 100+ | Covered |
| Active Guidelines | 90+ | Current |
| Genes Covered | 15+ | Primary PGx genes |
| Evidence Levels | Strong, Moderate, Limited | Classification |
| Implementation | G-Standaard | Dutch national |

---

## Entity Relationship Overview

```mermaid
---
config:
  layout: elk
---
flowchart LR
    accTitle: DPWG Pharmacogenomics Guideline Structure
    accDescr: Shows relationships between genes, guidelines, and drugs in Dutch pharmacogenomics guidelines

    %% Style definitions using Cagle palette
    classDef gene fill:#E8F5E9,stroke:#2E7D32,stroke-width:2px,color:#1B5E20
    classDef guideline fill:#E3F2FD,stroke:#1565C0,stroke-width:2px,color:#0D47A1
    classDef drug fill:#FFF8E1,stroke:#F57F17,stroke-width:2px,color:#E65100

    %% Entities with fields
    GENE["Gene<br/>────────────<br/>symbol<br/>phenotypes<br/>alleles"]:::gene
    GUIDE["Guideline<br/>────────────<br/>dpwg_id<br/>evidence<br/>urgency"]:::guideline
    DRUG["Drug<br/>────────────<br/>generic_name<br/>atc_code<br/>dose_adjust"]:::drug

    %% Relationships
    GENE --> GUIDE --> DRUG
```

---

## Core Tables/Entities

### Guideline

**Description:** DPWG recommendation

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| gene | string | Yes | HGNC symbol |
| drug | string | Yes | Generic name |
| atc_code | string | Yes | WHO ATC code |
| phenotype | string | Yes | Metabolizer status |
| recommendation | string | Yes | Action to take |
| urgency | string | Yes | Urgent/Non-urgent |
| evidence_level | string | Yes | Strong/Moderate/Limited |

### Phenotype Mapping

**Description:** Genotype to phenotype translation

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| gene | string | Yes | Gene symbol |
| diplotype | string | Yes | Allele combination |
| phenotype | string | Yes | Predicted phenotype |
| activity_score | float | No | When applicable |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | G-Standaard (Dutch) |
| International | PharmGKB JSON/TSV |
| Encoding | UTF-8 |
| Access | Web + API |

---

## Sample Record

```json
{
  "gene": "CYP2C19",
  "drug": "clopidogrel",
  "atc_code": "B01AC04",
  "phenotype": "CYP2C19 Poor Metabolizer",
  "recommendation": "Avoid clopidogrel. Use alternative antiplatelet (prasugrel, ticagrelor).",
  "urgency": "Urgent",
  "evidence_level": "Strong"
}
```

---

## Urgency Classification

| Level | Description | Action |
|-------|-------------|--------|
| Urgent | Therapeutic failure or adverse event likely | Immediate intervention |
| Non-urgent | Suboptimal therapy possible | Consider at next prescription |

---

## Glossary

| Term | Definition |
|------|------------|
| G-Standaard | Dutch national drug database |
| KNMP | Royal Dutch Pharmacists Association |
| ATC | Anatomical Therapeutic Chemical classification |

---

## References

1. https://www.pharmgkb.org/page/dpwg
2. https://www.knmp.nl/
3. Swen et al. (2011) Clin Pharmacol Ther. DOI: 10.1038/clpt.2011.34
