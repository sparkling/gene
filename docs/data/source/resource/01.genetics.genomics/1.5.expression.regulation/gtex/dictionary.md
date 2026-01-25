# GTEx - Data Dictionary

## Overview

This data dictionary documents the schema for GTEx (Genotype-Tissue Expression) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gtex |
| **Name** | GTEx |
| **Parent** | 1.5.expression.regulation |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### eQTL Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| variant_id | string | 1:1 | Yes | GTEx variant ID | chr1_10177_A_C_b38 |
| gene_id | string | 1:1 | Yes | Ensembl gene ID | ENSG00000141510.17 |
| tss_distance | integer | 1:1 | Yes | Distance to TSS | 50000 |
| ma_samples | integer | 1:1 | Yes | Minor allele samples | 150 |
| ma_count | integer | 1:1 | Yes | Minor allele count | 200 |
| pval_nominal | float | 1:1 | Yes | Nominal p-value | 1.5e-10 |
| slope | float | 1:1 | Yes | Effect size (beta) | 0.35 |
| slope_se | float | 1:1 | Yes | Standard error | 0.05 |
| qval | float | 1:1 | No | q-value (FDR) | 0.001 |

### Gene Expression

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_id | string | 1:1 | Yes | Ensembl gene ID | ENSG00000141510 |
| gene_symbol | string | 1:1 | Yes | Gene symbol | TP53 |
| tissue | string | 1:1 | Yes | Tissue name | Whole_Blood |
| TPM | float | 1:1 | Yes | Transcripts per million | 25.5 |
| read_count | integer | 1:1 | Yes | Raw read count | 1500 |

### Sample/Subject

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| SUBJID | string | 1:1 | Yes | Subject identifier | GTEX-1117F |
| SAMPID | string | 1:1 | Yes | Sample identifier | GTEX-1117F-0003 |
| SMTSD | string | 1:1 | Yes | Tissue description | Heart - Left Ventricle |
| SMTS | string | 1:1 | Yes | Tissue abbreviation | Heart |
| SEX | integer | 1:1 | Yes | Sex (1=M, 2=F) | 1 |
| AGE | string | 1:1 | Yes | Age bracket | 50-59 |
| DTHHRDY | integer | 1:1 | No | Hardy scale death | 0-4 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Variant ID | chr_pos_ref_alt_build | chr1_10177_A_C_b38 | GTEx variant format |
| Gene ID | ENSG + digits + version | ENSG00000141510.17 | Ensembl with version |
| Subject ID | GTEX- + alphanumeric | GTEX-1117F | Donor identifier |
| Sample ID | Subject + suffix | GTEX-1117F-0003 | Tissue sample |
| rsID | rs + digits | rs12345 | dbSNP cross-reference |

---

## Enumerations

### Tissues (54 total, selected)

| Tissue Code | Full Name | N Samples |
|-------------|-----------|-----------|
| Whole_Blood | Whole Blood | 755 |
| Muscle_Skeletal | Muscle - Skeletal | 803 |
| Skin_Sun_Exposed | Skin - Sun Exposed (Lower leg) | 701 |
| Adipose_Subcutaneous | Adipose - Subcutaneous | 663 |
| Artery_Tibial | Artery - Tibial | 663 |
| Thyroid | Thyroid | 653 |
| Nerve_Tibial | Nerve - Tibial | 619 |
| Lung | Lung | 578 |
| Esophagus_Mucosa | Esophagus - Mucosa | 555 |
| Brain_Cortex | Brain - Cortex | 255 |

### Sex Encoding

| Code | Meaning |
|------|---------|
| 1 | Male |
| 2 | Female |

### Hardy Scale (Death Classification)

| Code | Description |
|------|-------------|
| 0 | Ventilator case |
| 1 | Violent and fast death |
| 2 | Fast death of natural causes |
| 3 | Intermediate death |
| 4 | Slow death |

### Age Brackets

| Bracket | Range |
|---------|-------|
| 20-29 | 20-29 years |
| 30-39 | 30-39 years |
| 40-49 | 40-49 years |
| 50-59 | 50-59 years |
| 60-69 | 60-69 years |
| 70-79 | 70-79 years |

---

## Entity Relationships

### Subject to Sample
- **Cardinality:** 1:N
- **Description:** One donor provides multiple tissue samples
- **Key Fields:** SUBJID, SAMPID

### Sample to Tissue
- **Cardinality:** N:1
- **Description:** Multiple samples per tissue type
- **Key Fields:** SAMPID, SMTSD

### Variant to eQTL
- **Cardinality:** 1:N
- **Description:** One variant can be eQTL in multiple tissues
- **Key Fields:** variant_id, tissue

### Gene to Expression
- **Cardinality:** 1:N
- **Description:** One gene expressed across tissues
- **Key Fields:** gene_id, tissue

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GTEx | Genotype-Tissue Expression | Project name |
| eQTL | Expression Quantitative Trait Locus | Regulatory variant |
| sQTL | Splicing QTL | Splice-affecting variant |
| TPM | Transcripts Per Million | Expression unit |
| TSS | Transcription Start Site | Gene start |
| FDR | False Discovery Rate | Multiple testing |
| MAF | Minor Allele Frequency | Population metric |
| SMTSD | Sample Tissue Detail | Tissue description |
| DTHHRDY | Death Hardy Scale | Death classification |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| Ensembl | ENSG | Gene reference |
| dbSNP | rsID | Variant identifier |
| eQTLGen | eQTL | Blood comparison |
| gnomAD | Variant ID | Population frequency |
| ENCODE | Regulatory | Functional annotations |

---

## Data Versions

| Version | Release | Samples | Donors |
|---------|---------|---------|--------|
| V6 | 2016 | 8,555 | 449 |
| V7 | 2017 | 11,688 | 714 |
| V8 | 2019 | 17,382 | 948 |
| V10 | 2024 | 19,000+ | 980 |

---

## Data Quality Notes

1. **Cardinality:** One eQTL effect per variant-gene-tissue combination
2. **Tissue Coverage:** 54 tissues from ~1,000 donors
3. **Cis-eQTL Window:** 1Mb from TSS
4. **Access Levels:** dbGaP controlled + public summary stats
5. **Normalization:** TMM + inverse normal transformation
6. **Covariates:** PEER factors, genotype PCs, sex, platform
