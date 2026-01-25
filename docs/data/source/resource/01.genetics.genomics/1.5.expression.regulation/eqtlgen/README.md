---
id: eqtlgen
title: "eQTLGen"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: expression.regulation
tier: 2
status: active
last_updated: 2026-01-25
tags:
  - eqtl
  - expression
  - blood
  - meta-analysis
  - trans-eqtl
---

# eQTLGen

eQTLGen is a large-scale meta-analysis project that has identified expression quantitative trait loci (eQTLs) in blood samples from 31,684 individuals. It provides the most comprehensive catalog of blood eQTLs, including both cis-eQTLs (local regulatory effects) and trans-eQTLs (distal regulatory effects).

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | eQTLGen Consortium |
| **Website** | https://www.eqtlgen.org/ |
| **Update Frequency** | Phase releases |
| **Records** | 16,000+ cis-eQTLs |
| **Latest Release** | Phase I |

Phase I of eQTLGen identified 16,987 significant cis-eQTLs affecting 12,176 genes and 16,989 trans-eQTL associations. The large sample size enables detection of small-effect and low-frequency variant associations not detectable in smaller studies.

eQTLGen data is particularly valuable for interpreting GWAS signals, Mendelian randomization studies, and understanding gene regulatory networks. The consortium continues to expand with Phase II adding single-cell resolution data.

## Key Statistics

| Metric | Value |
|--------|-------|
| Individuals | 31,684 |
| Cis-eQTLs | 16,987 |
| Trans-eQTLs | 16,989 |
| Genes Tested | 19,250 |
| Cohorts | 37 |

## Primary Use Cases

1. Blood eQTL lookup
2. GWAS signal interpretation
3. Mendelian randomization
4. Gene network analysis

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| SNP | `RS ID` | `rs12345` |
| Gene | `Ensembl ID` | `ENSG00000012048` |
| P-value | `Significance` | `1e-100` |

## Limitations

- Blood tissue only (no solid tissues)
- Primarily European ancestry samples
- Bulk expression (no cell-type resolution in Phase I)
- Trans-eQTLs have higher false positive rate
- Limited to variants tested on arrays

## Data Quality Notes

eQTLGen applies stringent quality control including sample and variant filtering, population stratification correction, and multiple testing adjustment. Cis-eQTLs are reported at FDR < 0.05, while trans-eQTLs use more conservative thresholds. Effect sizes and standard errors are provided for colocalization and Mendelian randomization analyses.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [GTEx](../gtex/README.md) - Multi-tissue eQTLs
- [ENCODE](../encode/README.md) - Regulatory elements
- [GWAS Catalog](../gwas.catalog/README.md) - GWAS associations
