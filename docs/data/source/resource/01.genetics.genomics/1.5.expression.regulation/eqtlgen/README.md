---
id: eqtlgen
title: "eQTLGen"
type: data-source
category: genetics
subcategory: expression.regulation
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [eqtl, expression, blood, meta-analysis, trans-eqtl]
---

# eQTLGen

**Category:** [Genetics & Genomics](../../_index.md) > [Expression & Regulation](../_index.md)

## Overview

eQTLGen is a large-scale meta-analysis project that has identified expression quantitative trait loci (eQTLs) in blood samples from 31,684 individuals. It provides the most comprehensive catalog of blood eQTLs, including both cis-eQTLs (local regulatory effects) and trans-eQTLs (distal regulatory effects).

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| SNP | RS ID | rs12345 |
| Gene | Ensembl ID | ENSG00000012048 |
| P-value | Significance | 1e-100 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://www.eqtlgen.org/ | Portal |
| Downloads | Summary statistics | Full cis/trans results |
| API | REST endpoints | Programmatic access |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access |
| Commercial Use | Yes |
| Citation | Required |

## See Also

- [GTEx](../gtex/_index.md) - Multi-tissue eQTLs
- [ENCODE](../encode/_index.md) - Regulatory elements
- [GWAS Catalog](../gwas.catalog/_index.md) - GWAS associations
