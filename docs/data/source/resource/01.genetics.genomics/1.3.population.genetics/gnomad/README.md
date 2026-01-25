---
id: gnomad
title: "gnomAD"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: population.genetics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - population
  - frequency
  - exome
  - genome
  - constraint
---

# gnomAD

gnomAD (Genome Aggregation Database) is the largest publicly available collection of human exome and genome sequencing data, aggregating data from large-scale sequencing projects worldwide. It provides allele frequencies across diverse populations and gene-level constraint metrics essential for clinical variant interpretation.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Broad Institute |
| **Website** | https://gnomad.broadinstitute.org/ |
| **Update Frequency** | Major releases |
| **Records** | 786,000,000+ variants |
| **Latest Release** | v4.0 (November 2023) |

The database includes both exome sequencing (covering protein-coding regions) and whole genome sequencing data, with variants annotated using VEP and LOFTEE for loss-of-function assessment. gnomAD provides population-specific frequencies across multiple ancestry groups including African, Latino, Ashkenazi Jewish, East Asian, European, and South Asian populations.

Gene constraint metrics (pLI, LOEUF, missense Z-scores) quantify tolerance to different types of variants, enabling identification of genes under strong selective pressure. The structural variant dataset provides copy number and other SV frequencies.

## Key Statistics

| Metric | Value |
|--------|-------|
| Exomes | 730,947 |
| Genomes | 76,215 |
| Total Individuals | 807,162 |
| SNVs | 786,011,050 |
| Structural Variants | 1,108,839 |

## Primary Use Cases

1. Population allele frequency filtering
2. Rare variant identification
3. Gene constraint assessment
4. Structural variant frequency lookup

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Variant | `chr-pos-ref-alt` | `1-12345-A-G` |
| Gene | `HGNC/Ensembl` | `BRCA1, ENSG00000012048` |
| SV | `gnomAD-SV ID` | `gnomAD-SV_v3_DUP_1_1` |

## Limitations

- Excludes individuals with severe pediatric disease (selection bias)
- Ancestry groups unbalanced (historically European-heavy)
- No phenotype data available (aggregate frequencies only)
- Some regions have low coverage (telomeres, centromeres)
- Large file sizes for downloads (terabytes)

## Data Quality Notes

gnomAD applies stringent quality control filters including sample QC (contamination, sex concordance), variant QC (allele balance, call rate), and site-specific filters. Variants are annotated with quality flags (lcr, segdup) indicating potential issues. The v4 release improved ancestry inference and variant calling, particularly for non-European populations.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [dbSNP](../../1.1.variant.repositories/dbsnp/README.md) - RS identifiers
- [ClinVar](../../1.1.variant.repositories/clinvar/README.md) - Clinical significance
