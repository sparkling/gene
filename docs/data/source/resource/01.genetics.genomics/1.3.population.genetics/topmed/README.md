---
id: topmed
title: "TOPMed"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: population.genetics
tier: 2
status: active
last_updated: 2026-01-25
tags:
  - population
  - frequency
  - deep-sequencing
  - rare-variants
  - nhlbi
---

# TOPMed

TOPMed (Trans-Omics for Precision Medicine) is an NHLBI program providing deep whole genome sequencing (~30x coverage) of over 180,000 individuals from diverse ancestral backgrounds. It was designed specifically to discover rare genetic variants associated with heart, lung, blood, and sleep disorders.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | NHLBI / University of Michigan |
| **Website** | https://topmed.nhlbi.nih.gov/ |
| **Update Frequency** | Periodic freezes |
| **Records** | 400,000,000+ variants |
| **Latest Release** | Freeze 10 |

The program integrates genomic data with phenotypic information from multiple NIH-funded cohort studies, enabling genotype-phenotype association studies at unprecedented scale. TOPMed's imputation server provides free imputation services using the TOPMed reference panel, which captures rare variants better than previous panels.

The Bravo browser provides public access to allele frequencies for variants discovered in TOPMed, while individual-level data requires dbGaP authorization due to the associated phenotypic information.

## Key Statistics

| Metric | Value |
|--------|-------|
| Individuals | 180,000+ |
| Sequencing Depth | ~30x WGS |
| Variants | 400,000,000+ |
| Rare Variants | 97% of SNVs |
| Studies | 80+ |

## Primary Use Cases

1. Rare variant discovery and analysis
2. Genotype imputation (TOPMed panel)
3. Disease association studies
4. Multi-ethnic genetic research

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Sample | `NWD prefix` | `NWD123456` |
| Study | `phs accession` | `phs000956` |
| Freeze | `Numbered version` | `Freeze 8` |

## Limitations

- Individual-level data requires dbGaP approval (months-long process)
- Phenotype access varies by contributing study
- Primarily US-based cohorts
- Complex data use agreements for some studies
- Large file sizes require substantial storage

## Data Quality Notes

TOPMed variant calls undergo extensive quality control including sample QC, variant QC, and batch effect assessment. The deep sequencing enables confident calling of rare variants (MAF < 0.1%). Bravo browser provides aggregate frequencies without individual-level data, suitable for most annotation needs.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [gnomAD](../gnomad/README.md) - Aggregated frequencies
- [1000 Genomes](../1000.genomes/README.md) - Reference panel
- [UK Biobank](../uk.biobank/README.md) - Large cohort
