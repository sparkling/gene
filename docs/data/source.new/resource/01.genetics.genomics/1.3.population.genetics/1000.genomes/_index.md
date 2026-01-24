---
id: 1000-genomes
title: "1000 Genomes Project"
type: data-source
category: genetics
subcategory: population.genetics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [population, frequency, reference, diversity, global]
---

# 1000 Genomes Project

**Category:** [Genetics & Genomics](../../_index.md) > [Population Genetics](../_index.md)

## Overview

The 1000 Genomes Project was the first large-scale effort to catalogue human genetic variation, sequencing genomes from 2,504 individuals across 26 populations from five continental super-populations: Africa (AFR), Americas (AMR), East Asia (EAS), Europe (EUR), and South Asia (SAS).

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Sample | Population-based | HG00096 |
| Population | 3-letter code | GBR, YRI, CHB |
| Super-population | 3-letter code | EUR, AFR, EAS |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| IGSR | https://www.internationalgenome.org/ | Primary portal |
| FTP | ftp://ftp.1000genomes.ebi.ac.uk/ | Data files |
| AWS | s3://1000genomes | Cloud access |
| Ensembl | Data browser | Integrated view |

## License

| Aspect | Value |
|--------|-------|
| License | Fort Lauderdale Agreement |
| Commercial Use | Yes |
| Citation | Required |

## See Also

- [gnomAD](../gnomad/_index.md) - Expanded populations
- [TOPMed](../topmed/_index.md) - Deep sequencing
- [UK Biobank](../uk.biobank/_index.md) - Large cohort
