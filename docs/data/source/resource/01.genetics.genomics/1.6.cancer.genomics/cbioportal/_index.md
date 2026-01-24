---
id: cbioportal
title: "cBioPortal"
type: data-source
category: genetics
subcategory: cancer.genomics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [cancer, mutation, genomics, tcga, clinical]
---

# cBioPortal

**Category:** [Genetics & Genomics](../../_index.md) > [Cancer Genomics](../_index.md)

## Overview

cBioPortal for Cancer Genomics is an open-source platform providing visualization, analysis, and download of large-scale cancer genomics datasets. It hosts data from major cancer sequencing projects including TCGA, AACR Project GENIE, and institutional studies, making multidimensional cancer genomics data accessible to researchers and clinicians.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Study | study_id | brca_tcga |
| Sample | sample_id | TCGA-A1-A0SB-01 |
| Gene | HUGO symbol | TP53 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://www.cbioportal.org/ | Interactive analysis |
| REST API | /api/ endpoints | JSON responses |
| R Package | cgdsr | R integration |
| Python | cbioportalR | Python client |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access (ODbL) |
| Commercial Use | Yes |
| Citation | Required |

## See Also

- [Schema Documentation](./schema.md)
- [COSMIC](../cosmic/_index.md) - Somatic mutations
- [OncoKB](../oncokb/_index.md) - Precision oncology
