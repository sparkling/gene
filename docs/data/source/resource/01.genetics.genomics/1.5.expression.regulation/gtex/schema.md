---
id: schema-gtex
title: "GTEx Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: migrated
tags: [schema, database, expression, eqtl, tissue, rna-seq]
---

# GTEx Schema Documentation

**Document ID:** SCHEMA-GTEX
**Version:** 1.0
**Source Version:** V8 (2019)

---

## TL;DR

GTEx provides gene expression data and eQTL associations across 54 human tissues from 948 donors. Core data structures include expression matrices (GCT format), eQTL results (TSV), and comprehensive sample/subject annotations enabling tissue-specific regulatory analysis.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Donors | 948 | V8 release |
| Samples | 17,382 | RNA-seq |
| Tissues | 54 | Non-diseased |
| Genes Tested | 38,595 | Protein-coding + lncRNA |
| Cis-eQTLs | 3,300,000+ | All tissues |
| sQTLs | 2,800,000+ | Splicing QTLs |

---

## Entity Relationship Overview

```
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│    Subject    │────▶│    Sample     │────▶│  Expression   │
├───────────────┤     ├───────────────┤     ├───────────────┤
│ SUBJID        │     │ SAMPID        │     │ gene_id       │
│ SEX           │     │ SMTSD         │     │ TPM           │
│ AGE           │     │ SMRIN         │     │ read_counts   │
└───────────────┘     └───────────────┘     └───────────────┘
                             │
                             ▼
                      ┌───────────────┐
                      │     eQTL      │
                      ├───────────────┤
                      │ variant_id    │
                      │ gene_id       │
                      │ pval_nominal  │
                      └───────────────┘
```

---

## Core Tables/Entities

### Subject Phenotypes

**Description:** Donor-level demographic and phenotypic information

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| SUBJID | string | Yes | Subject identifier (GTEX-XXXXX) |
| SEX | integer | Yes | 1=Male, 2=Female |
| AGE | string | Yes | Age bracket (20-29, 30-39, etc.) |
| DTHHRDY | integer | No | Hardy scale death classification (0-4) |

### Sample Attributes

**Description:** Sample-level metadata including tissue and quality metrics

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| SAMPID | string | Yes | Sample identifier |
| SMTS | string | Yes | Tissue type (general) |
| SMTSD | string | Yes | Tissue type (detailed) |
| SMRIN | float | No | RIN (RNA Integrity Number) |
| SMNABTCH | string | No | Nucleic acid batch |
| SMGEBTCH | string | No | Genotype/expression batch |
| SMAFRZE | string | Yes | Analysis freeze (RNASEQ, WGS, etc.) |

### Expression Matrix (GCT Format)

**Description:** Gene-level expression values (TPM or read counts)

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Name | string | Yes | Ensembl gene ID (ENSG) |
| Description | string | Yes | Gene symbol |
| [SAMPID columns] | float | Yes | TPM or read count per sample |

### eQTL Summary Statistics

**Description:** Cis-eQTL associations per tissue

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| variant_id | string | Yes | chr_pos_ref_alt_build format |
| gene_id | string | Yes | Ensembl gene ID |
| tss_distance | integer | Yes | Distance to transcription start site |
| ma_samples | integer | Yes | Minor allele samples |
| ma_count | integer | Yes | Minor allele count |
| maf | float | Yes | Minor allele frequency |
| pval_nominal | float | Yes | Nominal p-value |
| slope | float | Yes | Effect size (beta) |
| slope_se | float | Yes | Standard error of slope |
| pval_nominal_threshold | float | No | Significance threshold |
| min_pval_nominal | float | No | Minimum p-value for gene |
| pval_beta | float | No | Beta-adjusted p-value |

### sQTL Summary Statistics

**Description:** Splicing QTL associations per tissue

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| variant_id | string | Yes | chr_pos_ref_alt_build format |
| phenotype_id | string | Yes | Intron cluster ID |
| tss_distance | integer | Yes | Distance to splice site |
| pval_nominal | float | Yes | Nominal p-value |
| slope | float | Yes | Effect size |
| slope_se | float | Yes | Standard error |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /gene/expression | GET | Gene expression by tissue |
| /eqtl/associations | GET | eQTL lookup by variant/gene |
| /tissue/list | GET | Available tissues |
| /variant/annotation | GET | Variant annotations |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | GCT (Gene Cluster Text) |
| Alternative | Parquet, BigQuery |
| Compression | gzip (.gz) |
| Encoding | UTF-8 |
| Reference | GRCh38/hg38 |

---

## Sample Record

```json
{
  "subject": {
    "SUBJID": "GTEX-1117F",
    "SEX": 2,
    "AGE": "60-69",
    "DTHHRDY": 0
  },
  "sample": {
    "SAMPID": "GTEX-1117F-0226-SM-5GZZ7",
    "SMTS": "Blood",
    "SMTSD": "Whole Blood",
    "SMRIN": 7.8,
    "SMAFRZE": "RNASEQ"
  },
  "expression": {
    "gene_id": "ENSG00000012048.23",
    "gene_symbol": "BRCA1",
    "TPM": 5.234
  },
  "eqtl": {
    "variant_id": "chr17_43044295_A_G_b38",
    "gene_id": "ENSG00000012048.23",
    "pval_nominal": 1.5e-12,
    "slope": 0.45,
    "maf": 0.15
  }
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| TPM | Transcripts Per Million - normalized expression unit |
| eQTL | Expression Quantitative Trait Locus - variant affecting gene expression |
| sQTL | Splicing QTL - variant affecting alternative splicing |
| RIN | RNA Integrity Number (1-10 scale) |
| GCT | Gene Cluster Text format for expression matrices |
| Cis-eQTL | Variant within 1Mb of target gene |
| LOEUF | Loss-of-function observed/expected upper bound fraction |

---

## References

1. https://gtexportal.org/
2. GTEx Consortium (2020) Science. DOI: 10.1126/science.aaz1776
3. https://storage.googleapis.com/adult-gtex/
