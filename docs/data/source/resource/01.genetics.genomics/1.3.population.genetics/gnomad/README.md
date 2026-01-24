---
id: gnomad
title: "gnomAD"
type: data-source
category: genetics
subcategory: population.genetics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [population, frequency, exome, genome, constraint]
---

# gnomAD

**Category:** [Genetics & Genomics](../../_index.md) > [Population Genetics](../_index.md)

## Overview

gnomAD (Genome Aggregation Database) is the largest publicly available collection of human exome and genome sequencing data, aggregating data from large-scale sequencing projects worldwide. It provides allele frequencies across diverse populations and gene-level constraint metrics essential for clinical variant interpretation.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Variant | chr-pos-ref-alt | 1-12345-A-G |
| Gene | HGNC/Ensembl | BRCA1, ENSG00000012048 |
| SV | gnomAD-SV ID | gnomAD-SV_v3_DUP_1_1 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Browser | https://gnomad.broadinstitute.org | Interactive search |
| GraphQL | gnomAD API | Programmatic access |
| Downloads | Google Cloud, AWS | VCF, Hail tables |

## License

| Aspect | Value |
|--------|-------|
| License | ODC-ODbL 1.0 |
| Commercial Use | Yes (with attribution) |
| Citation | Required |

## See Also

- [Schema Documentation](./schema.md)
- [dbSNP](../../1.1.variant.repositories/dbsnp/_index.md) - RS identifiers
- [ClinVar](../../1.1.variant.repositories/clinvar/_index.md) - Clinical significance
