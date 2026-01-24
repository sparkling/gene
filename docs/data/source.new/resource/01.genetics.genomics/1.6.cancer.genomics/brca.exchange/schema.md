---
id: schema-brca-exchange
title: "BRCA Exchange Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, database, cancer, brca1, brca2, hereditary]
---

# BRCA Exchange Schema Documentation

**Document ID:** SCHEMA-BRCA-EXCHANGE
**Version:** 1.0
**Source Version:** Current (continuously updated)

---

## TL;DR

BRCA Exchange aggregates BRCA1/BRCA2 variant classifications from multiple sources including ClinVar, LOVD, and ENIGMA consortium. Data includes expert-reviewed classifications, functional evidence, and source attributions for hereditary cancer risk assessment.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| BRCA1 Variants | 22,000+ | Aggregated |
| BRCA2 Variants | 30,000+ | Aggregated |
| Expert Reviewed | 15,000+ | ENIGMA classified |
| Data Sources | 10+ | Contributing |
| ENIGMA Classified | 5,000+ | Gold standard |

---

## Entity Relationship Overview

```
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│   Variant     │────▶│ Classification│────▶│    Source     │
├───────────────┤     ├───────────────┤     ├───────────────┤
│ gene          │     │ significance  │     │ name          │
│ hgvs_cdna     │     │ evidence      │     │ submitter     │
│ hgvs_protein  │     │ date          │     │ date          │
└───────────────┘     └───────────────┘     └───────────────┘
```

---

## Core Tables/Entities

### Variant

**Description:** BRCA1 or BRCA2 variant

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Gene_Symbol | string | Yes | BRCA1 or BRCA2 |
| Genomic_Coordinate_hg38 | string | Yes | chr:pos:ref:alt |
| HGVS_cDNA | string | Yes | Coding DNA change |
| HGVS_Protein | string | No | Protein change |
| Reference_Sequence | string | Yes | Transcript reference |
| Pathogenicity_expert | string | No | ENIGMA classification |
| Pathogenicity_all | string | Yes | Aggregated classification |

### Classification

**Description:** Clinical significance assessment

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Clinical_significance_ENIGMA | string | No | ENIGMA expert classification |
| Clinical_significance_ClinVar | string | No | ClinVar classification |
| Date_last_evaluated_ENIGMA | date | No | ENIGMA review date |
| Assertion_method_ENIGMA | string | No | Classification criteria |

### Source Data

**Description:** Contributing database information

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Source | string | Yes | Database name |
| Source_URL | string | No | Source link |
| Submitter_ClinVar | string | No | ClinVar submitter |
| SCV_ClinVar | string | No | ClinVar accession |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | TSV download |
| Alternative | JSON API |
| Encoding | UTF-8 |
| Access | Web + API |

---

## Sample Record

```json
{
  "Gene_Symbol": "BRCA1",
  "HGVS_cDNA": "c.68_69delAG",
  "HGVS_Protein": "p.Glu23ValfsTer17",
  "Genomic_Coordinate_hg38": "chr17:43124027:AG:A",
  "Pathogenicity_expert": "Pathogenic",
  "Clinical_significance_ENIGMA": "Pathogenic",
  "Clinical_significance_ClinVar": "Pathogenic",
  "Source": "ClinVar,ENIGMA,LOVD"
}
```

---

## Classification Categories

| Classification | Description |
|----------------|-------------|
| Pathogenic | Disease-causing |
| Likely pathogenic | >90% probability pathogenic |
| Uncertain significance | Insufficient evidence |
| Likely benign | >90% probability benign |
| Benign | Not disease-causing |
| Not yet reviewed | Awaiting ENIGMA review |

---

## Glossary

| Term | Definition |
|------|------------|
| ENIGMA | Evidence-based Network for Interpretation of Germline Mutant Alleles |
| LOVD | Leiden Open Variation Database |
| GA4GH | Global Alliance for Genomics and Health |

---

## References

1. https://brcaexchange.org/
2. Spurdle et al. (2012) Hum Mutat. DOI: 10.1002/humu.21628
