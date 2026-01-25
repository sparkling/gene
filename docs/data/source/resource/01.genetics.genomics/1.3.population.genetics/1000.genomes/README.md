---
id: 1000-genomes
title: "1000 Genomes Project"
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
  - reference
  - diversity
  - global
---

# 1000 Genomes Project

The 1000 Genomes Project was the first large-scale effort to catalogue human genetic variation, sequencing genomes from 2,504 individuals across 26 populations from five continental super-populations: Africa (AFR), Americas (AMR), East Asia (EAS), Europe (EUR), and South Asia (SAS).

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | IGSR (International Genome Sample Resource) |
| **Website** | https://www.internationalgenome.org/ |
| **Update Frequency** | Static (Phase 3 complete) |
| **Records** | 88,000,000+ variants |
| **Latest Release** | Phase 3 (2015), 30x resequencing (2020) |

The project established foundational methods for population genetics studies and created a reference panel widely used for genotype imputation. Phase 3 data includes approximately 88 million variants including SNPs, indels, and structural variants, providing allele frequencies across diverse global populations.

The 30x high-coverage dataset (released in 2020) provides improved accuracy for low-frequency variants and expanded structural variant detection. This resource remains essential for population stratification, ancestry inference, and variant filtering in both research and clinical applications.

## Key Statistics

| Metric | Value |
|--------|-------|
| Individuals | 2,504 |
| Populations | 26 |
| Super-populations | 5 |
| SNVs | 84,700,000+ |
| Indels | 3,600,000+ |

## Primary Use Cases

1. Genotype imputation reference panel
2. Population allele frequency filtering
3. Ancestry inference and PCA
4. Population genetics research

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Sample | `Population-based` | `HG00096` |
| Population | `3-letter code` | `GBR, YRI, CHB` |
| Super-population | `3-letter code` | `EUR, AFR, EAS` |

## Limitations

- Sample size small compared to newer resources (gnomAD, TOPMed)
- Limited rare variant detection (<0.5% MAF)
- Some populations underrepresented
- Low-coverage sequencing in Phase 3 (4-8x average)
- No phenotype or disease data available

## Data Quality Notes

The 30x high-coverage resequencing addresses many quality limitations of the original low-coverage Phase 3 data. Variant calls are validated through multiple calling pipelines and available in both GRCh37 and GRCh38 coordinates. The dataset remains the most widely used imputation reference panel due to its diverse population representation.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [gnomAD](../gnomad/README.md) - Expanded populations
- [TOPMed](../topmed/README.md) - Deep sequencing
- [UK Biobank](../uk.biobank/README.md) - Large cohort
