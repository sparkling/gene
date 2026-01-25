# CIViC - Data Dictionary

## Overview

This data dictionary documents the schema for CIViC (Clinical Interpretation of Variants in Cancer) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | civic |
| **Name** | CIViC |
| **Parent** | 1.6.cancer.genomics |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | CIViC variant ID | 12 |
| name | string | 1:1 | Yes | Variant name | V600E |
| gene_id | integer | 1:1 | Yes | CIViC gene ID | 5 |
| gene_symbol | string | 1:1 | Yes | Gene symbol | BRAF |
| variant_types | array | 1:N | Yes | Sequence ontology terms | missense_variant |
| hgvs_expressions | array | 1:N | No | HGVS notation | ENST00000288602.10:c.1799T>A |
| coordinates | object | 1:1 | No | Genomic coordinates | {chr: 7, start: 140753336} |

### Evidence Item

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | Evidence ID | 1234 |
| variant_id | integer | 1:1 | Yes | Associated variant | 12 |
| disease | string | 1:1 | Yes | Cancer type | Melanoma |
| doid | string | 1:1 | No | Disease Ontology ID | DOID:1909 |
| drugs | array | 1:N | No | Associated drugs | ["vemurafenib", "dabrafenib"] |
| evidence_type | string | 1:1 | Yes | Evidence category | Predictive |
| evidence_direction | string | 1:1 | Yes | Direction | Supports |
| evidence_level | string | 1:1 | Yes | Evidence strength | A |
| clinical_significance | string | 1:1 | Yes | Clinical meaning | Sensitivity |
| description | string | 1:1 | Yes | Evidence summary | V600E predicts response... |
| source_id | integer | 1:1 | Yes | Publication source | 12345 |

### Assertion

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | Assertion ID | 100 |
| variant_id | integer | 1:1 | Yes | Associated variant | 12 |
| disease | string | 1:1 | Yes | Cancer type | Melanoma |
| amp_level | string | 1:1 | Yes | AMP/ASCO tier | Tier I Level A |
| acmg_codes | array | 1:N | No | ACMG criteria | ["PS1", "PM2"] |
| summary | string | 1:1 | Yes | Assertion summary | BRAF V600E is a... |
| status | string | 1:1 | Yes | Curation status | accepted |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| CIViC Variant ID | Integer | 12 | Variant identifier |
| CIViC Evidence ID | Integer | 1234 | Evidence item |
| CIViC Gene ID | Integer | 5 | Gene identifier |
| DOID | DOID: + digits | DOID:1909 | Disease ontology |
| PubMed ID | Integer | 12345678 | Publication reference |

---

## Enumerations

### Evidence Type

| Type | Description |
|------|-------------|
| Predictive | Drug response prediction |
| Diagnostic | Cancer diagnosis/subtype |
| Prognostic | Outcome prediction |
| Predisposing | Cancer risk |
| Oncogenic | Driver status |
| Functional | Biological function |

### Evidence Level

| Level | Description |
|-------|-------------|
| A | Validated association |
| B | Clinical evidence |
| C | Case study |
| D | Preclinical evidence |
| E | Inferential association |

### Evidence Direction

| Direction | Meaning |
|-----------|---------|
| Supports | Evidence supports association |
| Does Not Support | Evidence refutes association |

### Clinical Significance

| Significance | Description |
|--------------|-------------|
| Sensitivity | Predicts drug response |
| Resistance | Predicts drug resistance |
| Adverse Response | Predicts toxicity |
| Reduced Sensitivity | Partial resistance |
| Positive | Favorable prognosis |
| Negative | Poor prognosis |
| Poor Outcome | Reduced survival |
| Better Outcome | Improved survival |
| Pathogenic | Disease-causing |
| Likely Pathogenic | Probably pathogenic |
| Benign | Not pathogenic |
| Oncogenic | Driver mutation |
| Gain of Function | Activating |
| Loss of Function | Inactivating |

### AMP/ASCO/CAP Tiers

| Tier | Level | Description |
|------|-------|-------------|
| Tier I | Level A | FDA-approved therapy |
| Tier I | Level B | Well-powered studies |
| Tier II | Level C | FDA-approved other tumor |
| Tier II | Level D | Preclinical/case reports |
| Tier III | - | Uncertain significance |
| Tier IV | - | Benign/likely benign |

---

## Entity Relationships

### Gene to Variant
- **Cardinality:** 1:N
- **Description:** One gene has multiple variants
- **Key Fields:** gene_id, variant_id

### Variant to Evidence
- **Cardinality:** 1:N
- **Description:** One variant has multiple evidence items
- **Key Fields:** variant_id, evidence_id

### Variant to Assertion
- **Cardinality:** 1:N
- **Description:** One variant has multiple assertions
- **Key Fields:** variant_id, assertion_id

### Evidence to Source
- **Cardinality:** N:1
- **Description:** Multiple evidence items from one publication
- **Key Fields:** evidence_id, source_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CIViC | Clinical Interpretation of Variants in Cancer | Database name |
| AMP | Association for Molecular Pathology | Guideline source |
| ASCO | American Society of Clinical Oncology | Guideline source |
| CAP | College of American Pathologists | Guideline source |
| ACMG | American College of Medical Genetics | Classification criteria |
| DOID | Disease Ontology Identifier | Disease ontology |
| SO | Sequence Ontology | Variant classification |
| HGVS | Human Genome Variation Society | Nomenclature |
| FDA | Food and Drug Administration | Regulatory approval |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| DOID | Disease Ontology | Cancer type |
| PubMed | PMID | Publication |
| ClinVar | VCV | Clinical significance |
| dbSNP | rsID | Variant identifier |
| Ensembl | ENST | Transcript |
| DrugBank | DB ID | Drug information |
| NCIt | NCIt ID | Cancer terminology |

---

## Data Quality Notes

1. **Cardinality:** One evidence item per variant-disease-drug-publication
2. **Curation:** Community-curated with expert review
3. **Evidence Levels:** A-E hierarchy based on evidence strength
4. **Assertions:** Expert-reviewed summaries of evidence
5. **Open Access:** Freely available under CC0
6. **Updates:** Continuous community curation
