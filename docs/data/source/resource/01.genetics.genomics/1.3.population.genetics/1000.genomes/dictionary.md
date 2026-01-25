# 1000 Genomes - Data Dictionary

## Overview

This data dictionary documents the schema for 1000 Genomes Project population reference database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | 1000.genomes |
| **Name** | 1000 Genomes Project |
| **Parent** | 1.3.population.genetics |
| **Total Fields** | 18 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| CHROM | string | 1:1 | Yes | Chromosome | 1, X, MT |
| POS | integer | 1:1 | Yes | Genomic position (1-based) | 10177 |
| ID | string | 1:1 | No | Variant identifier | rs367896724 |
| REF | string | 1:1 | Yes | Reference allele | A |
| ALT | string | 1:N | Yes | Alternate allele(s) | C, C,G |
| QUAL | float | 1:1 | No | Quality score | 100 |
| FILTER | string | 1:1 | Yes | Filter status | PASS |

### Population Frequencies

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| AF | float | 1:1 | Yes | Global allele frequency | 0.425 |
| AC | integer | 1:1 | Yes | Allele count | 2130 |
| AN | integer | 1:1 | Yes | Allele number (chromosomes) | 5008 |
| EAS_AF | float | 1:1 | No | East Asian frequency | 0.3363 |
| EUR_AF | float | 1:1 | No | European frequency | 0.4970 |
| AFR_AF | float | 1:1 | No | African frequency | 0.4909 |
| AMR_AF | float | 1:1 | No | American frequency | 0.3602 |
| SAS_AF | float | 1:1 | No | South Asian frequency | 0.4949 |

### Genotype Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| GT | string | 1:1 | Yes | Genotype | 0/0, 0/1, 1/1 |
| DP | integer | 1:1 | No | Read depth | 30 |
| GQ | integer | 1:1 | No | Genotype quality | 99 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| rsID | rs + digits | rs367896724 | dbSNP identifier |
| Sample ID | Population + digits | NA12878 | Individual sample |
| Population | 3-letter code | CEU, YRI | Population code |
| Super-population | 3-letter code | EUR, AFR | Continental group |

---

## Enumerations

### Super-populations

| Code | Name | Sample Count |
|------|------|--------------|
| AFR | African | 661 |
| AMR | Admixed American | 347 |
| EAS | East Asian | 504 |
| EUR | European | 503 |
| SAS | South Asian | 489 |

### Populations (Selected)

| Code | Name | Super-pop |
|------|------|-----------|
| CEU | Utah residents (CEPH) | EUR |
| GBR | British | EUR |
| YRI | Yoruba in Ibadan | AFR |
| CHB | Han Chinese Beijing | EAS |
| JPT | Japanese Tokyo | EAS |
| PUR | Puerto Rican | AMR |
| GIH | Gujarati Indian | SAS |

### Filter Status

| Value | Meaning |
|-------|---------|
| PASS | Passed all filters |
| LowQual | Low quality score |
| VQSRTrancheSNP | Failed VQSR for SNPs |
| VQSRTrancheINDEL | Failed VQSR for indels |

---

## Entity Relationships

### Variant to Population
- **Cardinality:** 1:N
- **Description:** Each variant has frequencies across populations
- **Key Fields:** CHROM, POS, population code

### Variant to Sample
- **Cardinality:** 1:N
- **Description:** Each variant genotyped in multiple samples
- **Key Fields:** CHROM, POS, sample_id

### Sample to Population
- **Cardinality:** N:1
- **Description:** Multiple samples per population
- **Key Fields:** sample_id, population

### Population to Super-population
- **Cardinality:** N:1
- **Description:** Multiple populations per continental group
- **Key Fields:** population, super_population

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| 1KG | 1000 Genomes Project | Project abbreviation |
| AF | Allele Frequency | Population metric |
| AC | Allele Count | Observed count |
| AN | Allele Number | Total chromosomes |
| MAF | Minor Allele Frequency | Less common allele |
| VCF | Variant Call Format | File format |
| VQSR | Variant Quality Score Recalibration | QC method |
| GRCh38 | Genome Reference Consortium human 38 | Current assembly |
| LD | Linkage Disequilibrium | Haplotype structure |
| CEPH | Centre d'Etude du Polymorphisme Humain | Sample source |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| gnomAD | Variant ID | Extended populations |
| Ensembl | ENSG | Gene annotation |
| UCSC | Genome Browser | Visualization |

---

## Population Codes Reference

### Full Population List

| Code | Full Name | N | Super-pop |
|------|-----------|---|-----------|
| ACB | African Caribbeans Barbados | 96 | AFR |
| ASW | Americans African Southwest US | 61 | AFR |
| BEB | Bengali Bangladesh | 86 | SAS |
| CDX | Chinese Dai Xishuangbanna | 93 | EAS |
| CEU | Utah CEPH | 99 | EUR |
| CHB | Han Chinese Beijing | 103 | EAS |
| CHS | Southern Han Chinese | 105 | EAS |
| CLM | Colombians Medellin | 94 | AMR |
| ESN | Esan Nigeria | 99 | AFR |
| FIN | Finnish | 99 | EUR |
| GBR | British | 91 | EUR |
| GIH | Gujarati Houston | 103 | SAS |
| GWD | Gambian Western Divisions | 113 | AFR |
| IBS | Iberian Spain | 107 | EUR |
| ITU | Indian Telugu UK | 102 | SAS |
| JPT | Japanese Tokyo | 104 | EAS |
| KHV | Kinh Vietnam | 99 | EAS |
| LWK | Luhya Kenya | 99 | AFR |
| MSL | Mende Sierra Leone | 85 | AFR |
| MXL | Mexican Los Angeles | 64 | AMR |
| PEL | Peruvians Lima | 85 | AMR |
| PJL | Punjabi Pakistan | 96 | SAS |
| PUR | Puerto Rican | 104 | AMR |
| STU | Sri Lankan Tamil UK | 102 | SAS |
| TSI | Toscani Italy | 107 | EUR |
| YRI | Yoruba Ibadan | 108 | AFR |

---

## Data Quality Notes

1. **Cardinality:** One frequency set per variant per population
2. **Sample Size:** 2,504 individuals from 26 populations
3. **Coverage:** Average 7x whole genome sequencing
4. **Phases:** Phase 3 is current final release
5. **Ancestry Bias:** Some populations underrepresented
