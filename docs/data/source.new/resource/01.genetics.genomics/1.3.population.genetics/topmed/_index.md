---
id: topmed
title: "TOPMed"
type: data-source
category: genetics
subcategory: population.genetics
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [population, frequency, deep-sequencing, rare-variants, nhlbi]
---

# TOPMed

**Category:** [Genetics & Genomics](../../_index.md) > [Population Genetics](../_index.md)

## Overview

TOPMed (Trans-Omics for Precision Medicine) is an NHLBI program providing deep whole genome sequencing (~30x coverage) of over 180,000 individuals from diverse ancestral backgrounds. It was designed specifically to discover rare genetic variants associated with heart, lung, blood, and sleep disorders.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Sample | NWD prefix | NWD123456 |
| Study | phs accession | phs000956 |
| Freeze | Numbered version | Freeze 8 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Bravo | https://bravo.sph.umich.edu/ | Public frequencies |
| Imputation Server | https://imputation.biodatacatalyst.nhlbi.nih.gov/ | Free imputation |
| dbGaP | https://www.ncbi.nlm.nih.gov/gap/ | Individual data |
| BioData Catalyst | NHLBI platform | Analysis platform |

## License

| Aspect | Value |
|--------|-------|
| License | dbGaP Data Use Agreement |
| Commercial Use | Varies by study |
| Public Frequencies | Freely available |

## See Also

- [gnomAD](../gnomad/_index.md) - Aggregated frequencies
- [1000 Genomes](../1000.genomes/_index.md) - Reference panel
- [UK Biobank](../uk.biobank/_index.md) - Large cohort
