# UK Biobank - Data Dictionary

## Overview

This data dictionary documents the schema for UK Biobank genetic data.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | uk.biobank |
| **Name** | UK Biobank |
| **Parent** | 1.3.population.genetics |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| CHROM | string | 1:1 | Yes | Chromosome | 1, X |
| POS | integer | 1:1 | Yes | Genomic position (GRCh38) | 10177 |
| ID | string | 1:1 | No | Variant identifier | rs367896724 |
| REF | string | 1:1 | Yes | Reference allele | A |
| ALT | string | 1:N | Yes | Alternate allele(s) | C |
| FILTER | string | 1:1 | Yes | Filter status | PASS |

### Allele Frequencies

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| AF | float | 1:1 | Yes | Allele frequency | 0.0001 |
| AC | integer | 1:1 | Yes | Allele count | 100 |
| AN | integer | 1:1 | Yes | Allele number | 974000 |
| n_het | integer | 1:1 | No | Heterozygote count | 98 |
| n_hom_alt | integer | 1:1 | No | Homozygote alt count | 1 |

### Quality Metrics

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| call_rate | float | 1:1 | No | Genotype call rate | 0.99 |
| HWE_p | float | 1:1 | No | Hardy-Weinberg p-value | 0.45 |
| info | float | 1:1 | No | Imputation quality | 0.98 |

### Phenotype Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| field_id | integer | 1:1 | Yes | UK Biobank field ID | 21001 |
| beta | float | 1:1 | No | Effect size | 0.05 |
| se | float | 1:1 | No | Standard error | 0.01 |
| pval | float | 1:1 | No | P-value | 5e-8 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| rsID | rs + digits | rs367896724 | dbSNP identifier |
| Participant ID | 7 digits | 1234567 | UK Biobank participant |
| Field ID | digits | 21001 | Phenotype field code |
| Data-Field | digits-digits | 21001-0.0 | Instance-array notation |
| Application ID | digits | 12345 | Research application |

---

## Enumerations

### Genotyping Array

| Value | Description | Samples |
|-------|-------------|---------|
| UK Biobank Axiom | Primary array | 438,000 |
| UK BiLEVE Axiom | Lung study subset | 50,000 |

### Data Type

| Value | Description |
|-------|-------------|
| Directly genotyped | From array |
| Imputed | Statistical imputation |
| WES | Whole exome sequencing |
| WGS | Whole genome sequencing |

### Ancestry (Self-reported)

| Code | Description | Approximate N |
|------|-------------|---------------|
| 1 | White | 430,000 |
| 2 | Mixed | 3,000 |
| 3 | Asian/Asian British | 10,000 |
| 4 | Black/Black British | 8,000 |
| 5 | Chinese | 1,500 |
| 6 | Other | 4,500 |

### Data Categories

| ID | Category |
|----|----------|
| 100313 | Genomics |
| 100314 | Exome Sequencing |
| 100315 | Whole Genome Sequencing |
| 100316 | Imputation |

---

## Entity Relationships

### Variant to Participant
- **Cardinality:** 1:N
- **Description:** Each variant genotyped across participants
- **Key Fields:** CHROM, POS, participant_id

### Participant to Phenotype
- **Cardinality:** 1:N
- **Description:** Each participant has many phenotypes
- **Key Fields:** participant_id, field_id

### Variant to Association
- **Cardinality:** 1:N
- **Description:** Variants tested against multiple phenotypes
- **Key Fields:** CHROM, POS, field_id

### Participant to Sample
- **Cardinality:** 1:N
- **Description:** Participants may have multiple samples
- **Key Fields:** participant_id, sample_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| UKB | UK Biobank | Database abbreviation |
| WES | Whole Exome Sequencing | 450K participants |
| WGS | Whole Genome Sequencing | 500K planned |
| GWAS | Genome-Wide Association Study | Analysis type |
| PRS | Polygenic Risk Score | Genetic prediction |
| HWE | Hardy-Weinberg Equilibrium | QC metric |
| MAF | Minor Allele Frequency | Population metric |
| QC | Quality Control | Data filtering |
| EID | Encrypted Identifier | Participant ID |
| RAP | Research Analysis Platform | Cloud platform |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| gnomAD | Variant ID | Frequency comparison |
| GWAS Catalog | Study ID | Association studies |
| Open Targets | Association | Drug targets |
| PheWAS Catalog | Phenotype | Phenome-wide associations |

---

## Phenotype Field Examples

| Field ID | Name | Type |
|----------|------|------|
| 21001 | Body mass index (BMI) | Continuous |
| 21002 | Weight | Continuous |
| 30000 | White blood cell count | Continuous |
| 41202 | Diagnoses (ICD-10) | Categorical |
| 20116 | Smoking status | Categorical |
| 2443 | Diabetes diagnosed | Binary |

---

## Data Quality Notes

1. **Cardinality:** One genotype per variant-participant pair
2. **Sample Size:** 500,000 participants (largest biobank)
3. **Access Control:** Requires approved research application
4. **Ancestry:** Predominantly British/European (~94%)
5. **Age Range:** 40-69 at recruitment (2006-2010)
6. **Imputation:** Reference panel: Haplotype Reference Consortium + UK10K
7. **Sequencing:** WES for 450K, WGS for 500K participants
