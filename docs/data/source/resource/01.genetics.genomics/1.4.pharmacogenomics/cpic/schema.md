---
id: schema-cpic
title: "CPIC Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, database, pharmacogenomics, guidelines, clinical]
---

# CPIC Schema Documentation

**Document ID:** SCHEMA-CPIC
**Version:** 1.0
**Source Version:** Current (continuously updated)

---

## TL;DR

CPIC provides structured clinical pharmacogenomics guidelines linking gene-drug pairs to therapeutic recommendations. Data includes allele function tables, diplotype-phenotype translations, and dosing recommendations accessed via PharmGKB.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Published Guidelines | 33 | CPIC website |
| Gene-Drug Pairs | 75+ | Covered |
| Priority Genes | 25 | Level A evidence |
| Evidence Levels | A, B, C, D | Classification |
| Update Cycle | Continuous | As evidence emerges |

---

## Entity Relationship Overview

```mermaid
---
config:
  layout: elk
---
flowchart TD
    accTitle: CPIC Pharmacogenomics Guideline Structure
    accDescr: Shows relationships between genes, guidelines, drugs, diplotypes, and therapeutic recommendations

    %% Style definitions using Cagle palette
    classDef gene fill:#E8F5E9,stroke:#2E7D32,stroke-width:2px,color:#1B5E20
    classDef guideline fill:#E3F2FD,stroke:#1565C0,stroke-width:2px,color:#0D47A1
    classDef drug fill:#FFF8E1,stroke:#F57F17,stroke-width:2px,color:#E65100
    classDef diplotype fill:#F3E5F5,stroke:#7B1FA2,stroke-width:2px,color:#4A148C
    classDef recommendation fill:#E1F5FE,stroke:#0277BD,stroke-width:2px,color:#01579B

    %% Entities with fields
    GENE["Gene<br/>────────────<br/>symbol<br/>alleles<br/>functions"]:::gene
    GUIDE["Guideline<br/>────────────<br/>DOI<br/>publication<br/>level"]:::guideline
    DRUG["Drug<br/>────────────<br/>generic_name<br/>ATC_code<br/>indications"]:::drug
    DIPLO["Diplotype<br/>────────────<br/>allele1<br/>allele2<br/>phenotype"]:::diplotype
    REC["Recommendation<br/>────────────<br/>phenotype<br/>dosing<br/>alternatives"]:::recommendation

    %% Relationships
    GENE --> GUIDE
    GUIDE --> DRUG
    GENE --> DIPLO
    DIPLO --> REC
    GUIDE --> REC
```

---

## Core Tables/Entities

### Guideline

**Description:** CPIC guideline document

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| guideline_id | string | Yes | Unique identifier |
| gene | string | Yes | HGNC symbol |
| drug | string | Yes | Generic drug name |
| publication_doi | string | Yes | DOI reference |
| cpic_level | string | Yes | A, B, C, or D |
| last_updated | date | Yes | Update date |

### Allele Function Table

**Description:** Allele-to-function mapping

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| gene | string | Yes | Gene symbol |
| allele | string | Yes | Star allele (e.g., *1) |
| function | string | Yes | Normal/Decreased/No function |
| clinical_function | string | No | Clinical interpretation |
| activity_score | float | No | Numeric activity |

### Recommendation

**Description:** Therapeutic recommendation

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| phenotype | string | Yes | Metabolizer status |
| drug | string | Yes | Drug name |
| recommendation | string | Yes | Dosing guidance |
| strength | string | Yes | Strong/Moderate/Optional |
| implications | string | No | Clinical implications |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | Excel supplements |
| Alternative | PharmGKB JSON/TSV |
| Encoding | UTF-8 |
| Access | Web + API |

---

## Sample Record

```json
{
  "gene": "CYP2D6",
  "drug": "codeine",
  "phenotype": "CYP2D6 Poor Metabolizer",
  "recommendation": "Avoid codeine use due to lack of efficacy. Consider alternative analgesics.",
  "strength": "Strong",
  "level": "A"
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| Level A | Genetic information should be used to change prescribing |
| Level B | Genetic information could be used to change prescribing |
| Activity Score | Numeric representation of enzyme activity |
| Diplotype | Combination of two alleles |

---

## References

1. https://cpicpgx.org/
2. https://www.pharmgkb.org/cpic
3. Relling & Klein (2011) Clin Pharmacol Ther. DOI: 10.1038/clpt.2010.279
