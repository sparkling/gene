---
id: hmp
title: "HMP - Human Microbiome Project"
type: data-source
category: microbiome
subcategory: body.site.microbiomes
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [microbiome, multi-site, reference, nih, metagenomics]
---

# HMP - Human Microbiome Project

**Category:** [Microbiome](../../../README.md) > [Body Site Microbiomes](../README.md)

## Overview

The Human Microbiome Project (HMP) is a NIH initiative that characterized the microbial communities from multiple body sites in healthy individuals. It established reference data and methods for studying the human microbiome.

HMP includes two phases: HMP1 (reference genomes, healthy adult microbiomes) and iHMP/HMP2 (longitudinal studies, disease states). The project characterized 5 major body areas: oral, skin, nasal, gastrointestinal, and urogenital.

HMP provides the foundational reference dataset for human microbiome research across all body sites.

## Key Statistics

| Metric | Value |
|--------|-------|
| Reference Genomes | 3,000+ |
| Body Sites | 18 |
| Subjects (HMP1) | 300 |
| Samples | 11,000+ |
| Sequence Data | 5.2 Tbp |

## Primary Use Cases

1. Reference microbiome characterization
2. Multi-body-site comparisons
3. Metagenomic assembly benchmarking
4. Species-level microbiome analysis
5. Method development and validation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Sample ID | HMP-specific | SRS011084 |
| Body Site | Controlled | Tongue dorsum |
| NCBI BioProject | `PRJNA[0-9]+` | PRJNA48479 |
| SRA Accession | `[SED]RR[0-9]+` | SRR346657 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| DACC Portal | https://hmpdacc.org | Main data access |
| HMP Cloud | https://portal.hmpdacc.org | Analysis portal |
| NCBI | https://www.ncbi.nlm.nih.gov/bioproject/PRJNA28331 | Raw data |
| Reference Genomes | https://hmpdacc.org/reference_genomes | Assembled genomes |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (NIH) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [Download Instructions](./download.md)
- [HOMD](../homd/README.md) - Oral microbiome detail
- [mBodyMap](../mbodymap/README.md) - Body site microbiome atlas
- [HMP Gut Data](../../9.1.gut.microbiome/hmp/README.md) - Stool-focused reference
- [MetaHIT](../../9.1.gut.microbiome/metahit/README.md) - European gut project
