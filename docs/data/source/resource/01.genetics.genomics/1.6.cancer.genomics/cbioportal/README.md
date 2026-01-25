---
id: cbioportal
title: "cBioPortal"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: cancer.genomics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - cancer
  - mutation
  - genomics
  - tcga
  - clinical
---

# cBioPortal

cBioPortal for Cancer Genomics is an open-source platform providing visualization, analysis, and download of large-scale cancer genomics datasets. It hosts data from major cancer sequencing projects including TCGA, AACR Project GENIE, and institutional studies, making multidimensional cancer genomics data accessible to researchers and clinicians.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Memorial Sloan Kettering Cancer Center / cBioPortal Team |
| **Website** | https://www.cbioportal.org/ |
| **Update Frequency** | Continuous study additions |
| **Records** | 220,890+ samples |
| **Latest Release** | Current (continuous) |

The portal integrates mutation, copy number alteration, expression, and clinical data across hundreds of cancer studies. Users can query genes across studies, visualize mutations in protein context, analyze co-occurrence patterns, and explore survival associations.

cBioPortal provides both web-based analysis tools and programmatic access through REST APIs, enabling integration into bioinformatics pipelines and clinical decision support systems.

## Key Statistics

| Metric | Value |
|--------|-------|
| Cancer Studies | 423 |
| Samples | 220,890 |
| Genes | 23,867 |
| Mutations | 20,000,000+ |
| Cancer Types | 100+ |

## Primary Use Cases

1. Cancer mutation profiling
2. Cross-study mutation analysis
3. Clinical trial matching
4. Biomarker discovery

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Study | `study_id` | `brca_tcga` |
| Sample | `sample_id` | `TCGA-A1-A0SB-01` |
| Gene | `HUGO symbol` | `TP53` |

## Limitations

- Study heterogeneity (different panels, platforms)
- Some studies have restricted clinical data
- Variant calling pipelines differ between studies
- Not all mutations are validated
- API rate limits apply for bulk queries

## Data Quality Notes

cBioPortal standardizes data from contributing studies through harmonization pipelines including variant normalization and clinical data mapping. Study-level metadata documents the sequencing platform, variant calling pipeline, and clinical annotation depth. Users should consider these factors when comparing across studies.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [COSMIC](../cosmic/README.md) - Somatic mutations
- [OncoKB](../oncokb/README.md) - Precision oncology
