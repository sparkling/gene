---
id: gtex
title: "GTEx"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: expression.regulation
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - expression
  - eqtl
  - tissue
  - rna-seq
  - splicing
---

# GTEx

GTEx (Genotype-Tissue Expression) is a comprehensive resource of gene expression and genetic regulation data from 54 non-diseased human tissues collected from nearly 1,000 individuals. It provides expression quantitative trait loci (eQTLs) and splicing QTLs (sQTLs) linking genetic variants to gene expression changes across tissues.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | GTEx Consortium / Broad Institute |
| **Website** | https://gtexportal.org/ |
| **Update Frequency** | Major releases |
| **Records** | 3,300,000+ eQTLs |
| **Latest Release** | V8 (August 2019) |

The V8 release includes RNA-seq data from 17,382 samples across 54 tissues with matched whole genome sequencing, enabling identification of tissue-specific regulatory effects. GTEx data is essential for interpreting GWAS signals and understanding tissue-specific disease mechanisms.

GTEx also provides gene expression levels (TPM values) for each tissue, enabling comparisons of gene activity across tissues and identification of tissue-specific genes. The portal allows interactive exploration of expression patterns and eQTL associations.

## Key Statistics

| Metric | Value |
|--------|-------|
| Donors | 948 |
| Samples | 17,382 |
| Tissues | 54 |
| eQTLs (V8) | 3,300,000+ |
| Genes Tested | 38,595 |

## Primary Use Cases

1. GWAS variant colocalization
2. Tissue-specific expression analysis
3. eQTL/sQTL lookup
4. Gene expression reference

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Sample | `GTEX-prefix` | `GTEX-1117F` |
| Gene | `Ensembl ID` | `ENSG00000012048` |
| Tissue | `Standardized name` | `Liver` |

## Limitations

- Non-diseased tissues only (no tumor samples)
- Deceased donors (postmortem ischemia effects)
- Sample sizes vary greatly by tissue
- European ancestry dominant
- Individual-level data requires dbGaP access

## Data Quality Notes

GTEx implements extensive quality control including RNA quality assessment (RIN > 6.0), sample identity verification, and expression outlier detection. eQTLs are reported at FDR < 0.05 with permutation-based significance testing. Aggregate data (eQTL summary statistics, median TPM) is freely available, while individual-level data requires dbGaP authorization.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [ENCODE](../encode/README.md) - Regulatory elements
- [eQTLGen](../eqtlgen/README.md) - Blood eQTLs
- [GWAS Catalog](../gwas.catalog/README.md) - GWAS colocalization
