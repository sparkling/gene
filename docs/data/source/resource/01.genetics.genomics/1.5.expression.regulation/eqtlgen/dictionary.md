# eQTLGen - Data Dictionary

## Overview

This data dictionary documents the schema for eQTLGen expression quantitative trait loci database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | eqtlgen |
| **Name** | eQTLGen |
| **Parent** | 1.5.expression.regulation |
| **Total Fields** | 15 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### eQTL Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Pvalue | float | 1:1 | Yes | Association p-value | 1.23e-50 |
| SNP | string | 1:1 | Yes | SNP identifier | rs12345 |
| SNPChr | integer | 1:1 | Yes | SNP chromosome | 1 |
| SNPPos | integer | 1:1 | Yes | SNP position | 10177 |
| AssessedAllele | string | 1:1 | Yes | Effect allele | A |
| OtherAllele | string | 1:1 | Yes | Non-effect allele | G |
| Zscore | float | 1:1 | Yes | Z-score | 15.2 |
| Gene | string | 1:1 | Yes | Gene symbol | TP53 |
| GeneSymbol | string | 1:1 | Yes | HGNC symbol | TP53 |
| GeneChr | integer | 1:1 | Yes | Gene chromosome | 17 |
| GenePos | integer | 1:1 | Yes | Gene TSS position | 7668402 |
| NrCohorts | integer | 1:1 | Yes | Number of cohorts | 37 |
| NrSamples | integer | 1:1 | Yes | Total sample size | 31684 |
| FDR | float | 1:1 | Yes | False discovery rate | 0.001 |

### cis/trans Classification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| cis_trans | string | 1:1 | Yes | cis or trans eQTL | cis |
| distance | integer | 1:1 | No | SNP-gene distance | 50000 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| rsID | rs + digits | rs12345 | dbSNP identifier |
| Gene Symbol | HGNC symbol | TP53 | Gene identifier |
| Ensembl Gene | ENSG + digits | ENSG00000141510 | Ensembl gene ID |
| Probe ID | ILMN_ + digits | ILMN_1234567 | Illumina array probe |

---

## Enumerations

### eQTL Type

| Type | Definition | Distance Cutoff |
|------|------------|-----------------|
| cis | Local regulation | SNP within 1Mb of gene |
| trans | Distal regulation | SNP >5Mb from gene or different chromosome |

### Significance Thresholds

| Threshold | P-value | Use Case |
|-----------|---------|----------|
| Genome-wide | <5e-8 | Standard GWAS |
| eQTL suggestive | <1e-5 | Exploratory |
| FDR 5% | FDR<0.05 | Multiple testing corrected |
| FDR 1% | FDR<0.01 | Stringent correction |

---

## Entity Relationships

### SNP to Gene
- **Cardinality:** N:M
- **Description:** SNPs can affect multiple genes; genes have multiple eQTLs
- **Key Fields:** SNP, Gene

### eQTL to Cohort
- **Cardinality:** 1:N
- **Description:** Each eQTL derived from multiple cohorts
- **Key Fields:** SNP, Gene, cohort_id

### Gene to Tissue
- **Cardinality:** 1:1
- **Description:** eQTLGen is blood-only
- **Key Fields:** Gene, tissue (blood)

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| eQTL | Expression Quantitative Trait Locus | Variant affecting expression |
| cis-eQTL | Local eQTL | Within 1Mb of gene |
| trans-eQTL | Distal eQTL | >5Mb or different chromosome |
| FDR | False Discovery Rate | Multiple testing correction |
| TSS | Transcription Start Site | Gene position reference |
| MAF | Minor Allele Frequency | Population metric |
| GWAS | Genome-Wide Association Study | Analysis type |
| SMR | Summary-based Mendelian Randomization | Causal inference |
| HEIDI | Heterogeneity in Dependent Instruments | Pleiotropy test |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| Ensembl | ENSG | Gene reference |
| GTEx | eQTL | Tissue comparison |
| GWAS Catalog | rsID | GWAS overlap |
| gnomAD | Variant ID | Population frequency |

---

## Data Subsets

| Dataset | Description | Size |
|---------|-------------|------|
| cis-eQTLs | Local regulatory variants | 19,250 genes |
| trans-eQTLs | Distal regulatory variants | 6,298 SNP-gene pairs |
| Full summary stats | All tested associations | ~10M SNP-gene pairs |
| SMR-ready | Formatted for SMR | Available |

---

## Data Quality Notes

1. **Cardinality:** One effect size per SNP-gene pair
2. **Tissue:** Blood (whole blood/PBMCs) only
3. **Sample Size:** 31,684 individuals from 37 cohorts
4. **Meta-analysis:** Inverse-variance weighted fixed effects
5. **Cis Window:** 1Mb from gene center
6. **Trans Threshold:** P < 5e-8 (genome-wide significant)
