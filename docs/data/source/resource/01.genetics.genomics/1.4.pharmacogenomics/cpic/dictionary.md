# CPIC - Data Dictionary

## Overview

This data dictionary documents the schema for CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | cpic |
| **Name** | CPIC |
| **Parent** | 1.4.pharmacogenomics |
| **Total Fields** | 15 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Guideline

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| guideline_id | string | 1:1 | Yes | CPIC guideline identifier | cpic:PA166104945 |
| gene | string | 1:1 | Yes | Gene symbol | CYP2D6 |
| drug | string | 1:1 | Yes | Drug name | codeine |
| publication_date | date | 1:1 | Yes | Publication date | 2014-01-01 |
| update_date | date | 1:1 | No | Last update date | 2019-03-15 |
| pmid | string | 1:1 | Yes | PubMed ID | 24458010 |
| url | string | 1:1 | Yes | Guideline URL | https://cpicpgx.org/... |

### Recommendation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| phenotype | string | 1:1 | Yes | Metabolizer phenotype | Poor Metabolizer |
| activity_score | string | 1:1 | No | Activity score range | 0 |
| therapeutic_recommendation | string | 1:1 | Yes | Clinical action | Avoid codeine |
| classification | string | 1:1 | Yes | Recommendation strength | Strong |
| implications | string | 1:1 | No | Clinical implications | Reduced morphine formation |

### Allele Function

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| allele | string | 1:1 | Yes | Star allele | *4 |
| function | string | 1:1 | Yes | Functional status | No function |
| activity_value | float | 1:1 | No | Numeric activity | 0 |
| clinical_function | string | 1:1 | No | Clinical interpretation | Non-functional |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Guideline ID | cpic:PA + digits | cpic:PA166104945 | CPIC guideline |
| Gene Symbol | HGNC symbol | CYP2D6 | Gene identifier |
| Star Allele | Gene*number | CYP2D6*4 | Haplotype |
| PMID | digits | 24458010 | PubMed reference |
| PharmGKB ID | PA + digits | PA166104945 | Cross-reference |

---

## Enumerations

### Phenotype Categories

| Phenotype | Activity Score | Description |
|-----------|----------------|-------------|
| Poor Metabolizer (PM) | 0 | No/minimal enzyme activity |
| Intermediate Metabolizer (IM) | 0.5-1.0 | Reduced enzyme activity |
| Normal Metabolizer (NM) | 1.5-2.0 | Normal enzyme activity |
| Rapid Metabolizer (RM) | 2.5 | Increased enzyme activity |
| Ultrarapid Metabolizer (UM) | >3.0 | Highly increased activity |

### Allele Functional Status

| Status | Description | Activity Value |
|--------|-------------|----------------|
| No function | Complete loss | 0 |
| Decreased function | Reduced activity | 0.5 |
| Normal function | Normal activity | 1 |
| Increased function | Enhanced activity | 1.5-2 |

### Recommendation Strength

| Classification | Meaning |
|----------------|---------|
| Strong | High confidence, should follow |
| Moderate | Moderate confidence |
| Optional | Consider with other factors |

### CPIC Levels

| Level | Meaning |
|-------|---------|
| A | Prescribing action recommended |
| B | Prescribing action potentially recommended |
| C | No prescribing action recommended |
| D | Evidence insufficient for action |

---

## Entity Relationships

### Guideline to Gene
- **Cardinality:** 1:N
- **Description:** One guideline may cover multiple genes
- **Key Fields:** guideline_id, gene

### Guideline to Drug
- **Cardinality:** 1:N
- **Description:** One guideline may cover multiple drugs
- **Key Fields:** guideline_id, drug

### Gene to Allele
- **Cardinality:** 1:N
- **Description:** Multiple alleles per gene
- **Key Fields:** gene, allele

### Phenotype to Recommendation
- **Cardinality:** 1:1
- **Description:** One recommendation per gene-drug-phenotype
- **Key Fields:** gene, drug, phenotype

### Allele to Diplotype
- **Cardinality:** N:N
- **Description:** Allele pairs form diplotypes
- **Key Fields:** allele1, allele2

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CPIC | Clinical Pharmacogenetics Implementation Consortium | Guideline source |
| PGx | Pharmacogenomics | Field abbreviation |
| PM | Poor Metabolizer | Phenotype |
| IM | Intermediate Metabolizer | Phenotype |
| NM | Normal Metabolizer | Phenotype |
| RM | Rapid Metabolizer | Phenotype |
| UM | Ultrarapid Metabolizer | Phenotype |
| AS | Activity Score | Phenotype calculation |
| CYP | Cytochrome P450 | Enzyme family |
| PharmVar | Pharmacogene Variation | Allele database |
| PMID | PubMed Identifier | Reference ID |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PharmGKB | PA accession | Guideline repository |
| DPWG | Guideline | Parallel guidelines |
| PharmVar | Allele | Star allele definitions |
| PubMed | PMID | Publication reference |
| ClinGen | Gene | Gene validity |

---

## Guideline Categories

| Category | Example Drugs | Gene(s) |
|----------|---------------|---------|
| Opioids | codeine, tramadol | CYP2D6 |
| Antiplatelet | clopidogrel | CYP2C19 |
| Anticoagulants | warfarin | CYP2C9, VKORC1, CYP4F2 |
| Antidepressants | SSRIs, TCAs | CYP2D6, CYP2C19 |
| Immunosuppressants | tacrolimus | CYP3A5 |
| Oncology | fluoropyrimidines | DPYD |
| Thiopurines | azathioprine | TPMT, NUDT15 |
| Statins | simvastatin | SLCO1B1 |

---

## Data Quality Notes

1. **Cardinality:** One recommendation per gene-drug-phenotype combination
2. **Evidence Basis:** Systematic literature review + expert consensus
3. **Updates:** Guidelines revised as evidence evolves
4. **Implementation:** Designed for EHR clinical decision support
5. **Harmonization:** Coordinated with DPWG where possible
6. **Access:** Guidelines freely available, CC BY-SA license
