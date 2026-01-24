---
id: schema-eqtlgen
title: "eQTLGen Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, database, eqtl, expression, blood]
---

# eQTLGen Schema Documentation

**Document ID:** SCHEMA-EQTLGEN
**Version:** 1.0
**Source Version:** Phase I (2019)

---

## TL;DR

eQTLGen provides cis- and trans-eQTL summary statistics from 31,684 blood samples. Data includes SNP-gene associations with effect sizes, p-values, and allele frequencies for interpreting regulatory effects of genetic variants.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Individuals | 31,684 | Phase I |
| Cis-eQTLs | 16,987 | Significant |
| Trans-eQTLs | 16,989 | Significant |
| Genes Tested | 19,250 | Expressed in blood |
| Cohorts | 37 | Contributing |

---

## Entity Relationship Overview

```
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│     SNP       │────▶│    eQTL       │────▶│    Gene       │
├───────────────┤     ├───────────────┤     ├───────────────┤
│ rsid          │     │ p-value       │     │ ensembl_id    │
│ chr:pos       │     │ z-score       │     │ symbol        │
│ alleles       │     │ effect_size   │     │ distance      │
└───────────────┘     └───────────────┘     └───────────────┘
```

---

## Core Tables/Entities

### Cis-eQTL Summary Statistics

**Description:** Local regulatory associations (SNP within 1Mb of gene)

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| SNP | string | Yes | RS ID |
| SNPChr | integer | Yes | SNP chromosome |
| SNPPos | integer | Yes | SNP position (hg19) |
| AssessedAllele | string | Yes | Effect allele |
| OtherAllele | string | Yes | Non-effect allele |
| Zscore | float | Yes | Association Z-score |
| Gene | string | Yes | Ensembl gene ID |
| GeneSymbol | string | Yes | HGNC symbol |
| GeneChr | integer | Yes | Gene chromosome |
| GenePos | integer | Yes | Gene TSS position |
| NrCohorts | integer | Yes | Contributing cohorts |
| NrSamples | integer | Yes | Sample size |
| FDR | float | Yes | False discovery rate |

### Trans-eQTL Summary Statistics

**Description:** Distal regulatory associations (SNP >5Mb from gene or different chr)

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| SNP | string | Yes | RS ID |
| SNPChr | integer | Yes | SNP chromosome |
| SNPPos | integer | Yes | SNP position |
| AssessedAllele | string | Yes | Effect allele |
| OtherAllele | string | Yes | Non-effect allele |
| Zscore | float | Yes | Association Z-score |
| Gene | string | Yes | Ensembl gene ID |
| GeneSymbol | string | Yes | HGNC symbol |
| GeneChr | integer | Yes | Gene chromosome |
| Pvalue | float | Yes | Association p-value |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | gzipped TSV |
| Alternative | Parquet (large files) |
| Encoding | UTF-8 |
| Reference | hg19/GRCh37 |

---

## Sample Record

```tsv
SNP	SNPChr	SNPPos	AssessedAllele	OtherAllele	Zscore	Gene	GeneSymbol	GeneChr	GenePos	NrCohorts	NrSamples	FDR
rs12345	1	100000	A	G	15.5	ENSG00000012048	BRCA1	17	43044295	37	31684	1e-50
```

---

## Glossary

| Term | Definition |
|------|------------|
| Cis-eQTL | SNP affecting expression of nearby gene (<1Mb) |
| Trans-eQTL | SNP affecting expression of distant gene |
| Z-score | Standardized effect size |
| FDR | Benjamini-Hochberg false discovery rate |

---

## References

1. https://www.eqtlgen.org/
2. Vosa et al. (2021) Nat Genet. DOI: 10.1038/s41588-021-00913-z
