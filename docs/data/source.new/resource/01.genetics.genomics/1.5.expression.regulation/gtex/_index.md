---
id: gtex
title: "GTEx"
type: data-source
category: genetics
subcategory: expression.regulation
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [expression, eqtl, tissue, rna-seq, splicing]
---

# GTEx

**Category:** [Genetics & Genomics](../../_index.md) > [Expression & Regulation](../_index.md)

## Overview

GTEx (Genotype-Tissue Expression) is a comprehensive resource of gene expression and genetic regulation data from 54 non-diseased human tissues collected from nearly 1,000 individuals. It provides expression quantitative trait loci (eQTLs) and splicing QTLs (sQTLs) linking genetic variants to gene expression changes across tissues.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Sample | GTEX-prefix | GTEX-1117F |
| Gene | Ensembl ID | ENSG00000012048 |
| Tissue | Standardized name | Liver |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Portal | https://gtexportal.org/ | Interactive browser |
| API | REST endpoints | Programmatic access |
| Downloads | Portal/dbGaP | Protected data via dbGaP |
| Google Cloud | BigQuery | Bulk analysis |

## License

| Aspect | Value |
|--------|-------|
| License | dbGaP agreement (individual) |
| Aggregate Data | Open access |
| Commercial Use | Per dbGaP terms |

## See Also

- [ENCODE](../encode/_index.md) - Regulatory elements
- [eQTLGen](../eqtlgen/_index.md) - Blood eQTLs
- [GWAS Catalog](../gwas.catalog/_index.md) - GWAS colocalization
