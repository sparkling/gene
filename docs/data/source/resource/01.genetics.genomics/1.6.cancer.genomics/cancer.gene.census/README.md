---
id: cancer-gene-census
title: "Cancer Gene Census"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: cancer.genomics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - cancer
  - driver
  - oncogene
  - tumor-suppressor
  - census
---

# Cancer Gene Census

The Cancer Gene Census (CGC) is an ongoing effort to catalogue genes whose mutations are causally implicated in cancer. Maintained by the Wellcome Sanger Institute as part of COSMIC, it distinguishes between Tier 1 genes with documented cancer-causing mutations and Tier 2 genes with strong evidence of cancer involvement.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Wellcome Sanger Institute |
| **Website** | https://cancer.sanger.ac.uk/census |
| **Update Frequency** | Periodic updates |
| **Records** | 736 genes |
| **Latest Release** | Current (within COSMIC) |

Each CGC entry includes the gene's role (oncogene, tumor suppressor, or fusion partner), associated cancer types, mutation types commonly observed, and supporting evidence from the literature. The census provides a gold-standard reference for cancer gene identification and driver mutation analysis.

The CGC is continuously updated as new cancer genes are identified through research and clinical studies. It serves as a foundation for cancer gene panel design, variant interpretation, and research prioritization.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Genes | 736 |
| Tier 1 Genes | 588 |
| Tier 2 Genes | 148 |
| Oncogenes | 300+ |
| Tumor Suppressors | 300+ |

## Primary Use Cases

1. Cancer panel gene selection
2. Driver mutation identification
3. Oncogene/TSG classification
4. Research prioritization

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Gene | `HUGO symbol` | `TP53` |
| Tier | `1 or 2` | `Tier 1` |
| Role | `Oncogene/TSG/Fusion` | `TSG` |

## Limitations

- Registration required for full data access
- Commercial use requires COSMIC license
- Limited to genes with documented driver mutations
- Does not include passenger mutation genes
- Cancer type associations may be incomplete

## Data Quality Notes

CGC genes are curated based on published evidence of driver mutation status. Tier 1 genes have extensive documented evidence including functional studies, while Tier 2 genes have strong but less comprehensive evidence. The cancer type associations reflect predominant patterns and may not include all reported cancer contexts.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [COSMIC](../cosmic/README.md) - Somatic mutations
- [OncoKB](../oncokb/README.md) - Actionable genes
- [CIViC](../civic/README.md) - Clinical interpretations
