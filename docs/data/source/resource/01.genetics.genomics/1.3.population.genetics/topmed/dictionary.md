# TOPMed - Data Dictionary

## Overview

This data dictionary documents the schema for TOPMed (Trans-Omics for Precision Medicine) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | topmed |
| **Name** | TOPMed |
| **Parent** | 1.3.population.genetics |
| **Total Fields** | 16 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| CHROM | string | 1:1 | Yes | Chromosome | chr1, chrX |
| POS | integer | 1:1 | Yes | Genomic position | 10177 |
| REF | string | 1:1 | Yes | Reference allele | A |
| ALT | string | 1:N | Yes | Alternate allele(s) | C |
| FILTER | string | 1:1 | Yes | Filter status | PASS |
| QUAL | float | 1:1 | No | Quality score | 100 |

### Allele Frequencies

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| AF | float | 1:1 | Yes | Overall allele frequency | 0.0001 |
| AC | integer | 1:1 | Yes | Allele count | 50 |
| AN | integer | 1:1 | Yes | Allele number (chromosomes) | 364000 |
| Het | integer | 1:1 | No | Heterozygote count | 48 |
| Hom | integer | 1:1 | No | Homozygote count | 1 |

### Quality Metrics

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| DP | integer | 1:1 | No | Total read depth | 500000 |
| MQ | float | 1:1 | No | Mapping quality | 60 |
| VRT | integer | 1:1 | No | Variant read threshold | 5 |
| ABHet | float | 1:1 | No | Allele balance heterozygotes | 0.5 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| TOPMed ID | chr-pos-ref-alt | chr1-10177-A-C | Variant identifier |
| Freeze | Freeze + number | Freeze 8 | Data release version |
| Study | Study name | MESA, FHS | Contributing study |
| Sample | NWD + 10 digits | NWD123456789 | Sample identifier |

---

## Enumerations

### Filter Status

| Value | Meaning |
|-------|---------|
| PASS | Passed all quality filters |
| SVM | Failed SVM filter |
| DISC | Mendelian/duplicate discordance |
| EXHET | Excess heterozygosity |

### Ancestry Groups

| Code | Name | Approximate % |
|------|------|---------------|
| EUR | European | 45% |
| AFR | African | 30% |
| HIS | Hispanic/Latino | 15% |
| EAS | East Asian | 5% |
| SAS | South Asian | 3% |
| Other | Other/Mixed | 2% |

### Contributing Studies

| Study | Focus | Samples |
|-------|-------|---------|
| MESA | Multi-Ethnic Atherosclerosis | 6,800 |
| FHS | Framingham Heart Study | 4,100 |
| JHS | Jackson Heart Study | 3,400 |
| ARIC | Atherosclerosis Risk | 15,800 |
| WHI | Women's Health Initiative | 11,000 |
| COPDGene | COPD Genetics | 10,200 |
| HVH | Heart and Vascular Health | 3,500 |

---

## Entity Relationships

### Variant to Sample
- **Cardinality:** 1:N
- **Description:** Each variant genotyped across all samples
- **Key Fields:** CHROM, POS, sample_id

### Sample to Study
- **Cardinality:** N:1
- **Description:** Multiple samples per contributing study
- **Key Fields:** sample_id, study_name

### Variant to Annotation
- **Cardinality:** 1:N
- **Description:** Variants annotated with functional predictions
- **Key Fields:** CHROM, POS, annotation_type

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| TOPMed | Trans-Omics for Precision Medicine | Program name |
| WGS | Whole Genome Sequencing | Data type |
| AF | Allele Frequency | Population metric |
| AC | Allele Count | Observed count |
| AN | Allele Number | Total chromosomes |
| SVM | Support Vector Machine | Filter method |
| NHLBI | National Heart, Lung, Blood Institute | Funding agency |
| dbGaP | Database of Genotypes and Phenotypes | Access portal |
| VCF | Variant Call Format | File format |
| CRAM | Compressed Reference-oriented Alignment Map | Alignment format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| gnomAD | Variant ID | Frequency comparison |
| dbGaP | phs number | Study access |
| ClinVar | VCV | Clinical significance |

---

## Freeze Versions

| Freeze | Samples | Variants | Date |
|--------|---------|----------|------|
| Freeze 5 | 53,800 | 410M | 2018 |
| Freeze 6 | 92,000 | 480M | 2019 |
| Freeze 8 | 132,000 | 570M | 2021 |
| Freeze 9 | 180,000+ | 600M+ | 2023 |

---

## Data Quality Notes

1. **Cardinality:** One frequency per variant across full cohort
2. **Sequencing:** Deep WGS (30x average coverage)
3. **Access Control:** Requires dbGaP authorization
4. **Ancestry:** More diverse than many reference panels
5. **Novel Variants:** ~97% of variants are rare (AF < 0.5%)
6. **Phenotypes:** Linked to rich clinical/phenotype data
