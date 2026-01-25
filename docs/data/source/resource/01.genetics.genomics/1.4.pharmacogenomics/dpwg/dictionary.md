# DPWG - Data Dictionary

## Overview

This data dictionary documents the schema for DPWG (Dutch Pharmacogenetics Working Group) guidelines.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dpwg |
| **Name** | DPWG |
| **Parent** | 1.4.pharmacogenomics |
| **Total Fields** | 12 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Guideline

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| guideline_id | string | 1:1 | Yes | Unique guideline identifier | PA166104945 |
| gene | string | 1:1 | Yes | Gene symbol | CYP2D6 |
| drug | string | 1:1 | Yes | Drug name | codeine |
| phenotype | string | 1:1 | Yes | Metabolizer phenotype | Poor Metabolizer |
| recommendation | string | 1:1 | Yes | Clinical recommendation | Avoid codeine |
| recommendation_strength | string | 1:1 | No | Evidence strength | Strong |

### Gene-Drug Pair

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_symbol | string | 1:1 | Yes | HGNC gene symbol | CYP2C19 |
| drug_name | string | 1:1 | Yes | Generic drug name | clopidogrel |
| atc_code | string | 1:1 | No | ATC classification | B01AC04 |
| phenotypes | array | 1:N | Yes | Relevant phenotypes | ["PM", "IM", "UM"] |

### Phenotype Definition

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| phenotype_code | string | 1:1 | Yes | Phenotype abbreviation | PM, IM, NM, UM |
| phenotype_name | string | 1:1 | Yes | Full phenotype name | Poor Metabolizer |
| activity_score | float | 1:1 | No | Activity score range | 0, 0.5, 1.0-2.0 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Guideline ID | PA + digits | PA166104945 | PharmGKB accession |
| Gene Symbol | HGNC symbol | CYP2D6 | Gene identifier |
| Drug ID | PA + digits | PA449053 | PharmGKB drug ID |
| ATC Code | Letter + digits | B01AC04 | WHO drug classification |

---

## Enumerations

### Phenotype Categories

| Code | Name | Description |
|------|------|-------------|
| PM | Poor Metabolizer | No/minimal enzyme activity |
| IM | Intermediate Metabolizer | Reduced enzyme activity |
| NM | Normal Metabolizer | Normal enzyme activity |
| UM | Ultrarapid Metabolizer | Increased enzyme activity |
| RM | Rapid Metabolizer | Moderately increased activity |

### Recommendation Strength

| Value | Meaning |
|-------|---------|
| Strong | High confidence, should be followed |
| Moderate | Moderate confidence |
| Optional | Consider if additional factors |
| None | No action required |

### Action Types

| Action | Description |
|--------|-------------|
| Avoid | Do not prescribe |
| Reduce dose | Use lower dose |
| Increase dose | Use higher dose |
| Alternative | Use different drug |
| Monitor | Enhanced monitoring |
| Standard | No change needed |

---

## Entity Relationships

### Guideline to Gene
- **Cardinality:** N:1
- **Description:** Multiple guidelines per gene
- **Key Fields:** guideline_id, gene

### Guideline to Drug
- **Cardinality:** N:1
- **Description:** Multiple guidelines per drug
- **Key Fields:** guideline_id, drug

### Gene to Phenotype
- **Cardinality:** 1:N
- **Description:** One gene maps to multiple phenotypes
- **Key Fields:** gene, phenotype

### Phenotype to Recommendation
- **Cardinality:** 1:1
- **Description:** One recommendation per gene-drug-phenotype
- **Key Fields:** gene, drug, phenotype

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DPWG | Dutch Pharmacogenetics Working Group | Guideline source |
| KNMP | Royal Dutch Pharmacists Association | Parent organization |
| PGx | Pharmacogenomics | Field abbreviation |
| PM | Poor Metabolizer | Phenotype |
| IM | Intermediate Metabolizer | Phenotype |
| NM | Normal Metabolizer | Phenotype |
| UM | Ultrarapid Metabolizer | Phenotype |
| CYP | Cytochrome P450 | Enzyme family |
| ATC | Anatomical Therapeutic Chemical | Drug classification |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PharmGKB | PA accession | Guideline repository |
| CPIC | Guideline ID | Parallel guidelines |
| DrugBank | DB ID | Drug information |
| HGNC | Symbol | Gene naming |
| PharmVar | Allele | Star allele definitions |

---

## Key Gene-Drug Pairs

| Gene | Example Drugs | Clinical Impact |
|------|---------------|-----------------|
| CYP2D6 | codeine, tramadol, tamoxifen | Efficacy/toxicity |
| CYP2C19 | clopidogrel, PPIs, voriconazole | Efficacy/toxicity |
| CYP2C9 | warfarin, NSAIDs | Toxicity/dosing |
| VKORC1 | warfarin | Dosing |
| DPYD | fluoropyrimidines | Severe toxicity |
| TPMT | azathioprine, mercaptopurine | Severe toxicity |
| UGT1A1 | irinotecan | Toxicity |

---

## Data Quality Notes

1. **Cardinality:** One recommendation per gene-drug-phenotype combination
2. **Updates:** Guidelines revised as evidence evolves
3. **Harmonization:** Largely consistent with CPIC guidelines
4. **Implementation:** Designed for clinical decision support
5. **Evidence Basis:** Systematic literature review
6. **Language:** Originally Dutch, translated to English
