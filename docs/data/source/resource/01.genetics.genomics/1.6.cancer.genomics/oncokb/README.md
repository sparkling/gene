---
id: oncokb
title: "OncoKB"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: cancer.genomics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - cancer
  - precision-oncology
  - actionable
  - therapy
  - clinical
---

# OncoKB

OncoKB is a precision oncology knowledge base developed by Memorial Sloan Kettering Cancer Center. It provides expert-curated annotation of the biological and clinical effects of cancer variants, linking alterations to cancer types, therapeutic implications, and clinical evidence levels.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Memorial Sloan Kettering Cancer Center |
| **Website** | https://www.oncokb.org/ |
| **Update Frequency** | Continuous |
| **Records** | 5,000+ annotated variants |
| **Latest Release** | Current (continuous) |

Each variant annotation includes functional effect classification (oncogenic, likely oncogenic, etc.), associated cancer types, and therapeutic implications with FDA/guideline-based evidence levels (1-4). OncoKB covers mutations, copy number alterations, fusions, and tumor mutational burden across hundreds of genes.

OncoKB is integrated into clinical molecular tumor boards and commercial molecular profiling platforms. The FDA has recognized OncoKB as a public tumor mutation database supporting clinical decision-making.

## Key Statistics

| Metric | Value |
|--------|-------|
| Genes | 700+ |
| Variants Annotated | 5,000+ |
| Cancer Types | 100+ |
| Level 1 Evidence | 100+ |
| Therapeutic Implications | 3,000+ |

## Primary Use Cases

1. Clinical variant interpretation
2. Therapy selection support
3. Clinical trial matching
4. Molecular tumor board support

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Gene | `HUGO symbol` | `BRAF` |
| Variant | `Protein change` | `V600E` |
| Evidence Level | `1-4, R1-R2` | `Level 1` |

## Limitations

- Commercial use requires license agreement
- Annotations focused on actionable variants
- May not cover all rare variants
- Evidence levels change with FDA approvals
- API requires registration

## Data Quality Notes

OncoKB annotations are curated by domain experts at MSK and reviewed according to a structured evidence framework. Level 1 evidence indicates FDA-recognized biomarkers, while lower levels indicate emerging or investigational evidence. Annotations are updated as new therapies gain approval and clinical trial data mature.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [CIViC](../civic/README.md) - Community curation
- [cBioPortal](../cbioportal/README.md) - Integrated annotations
- [COSMIC](../cosmic/README.md) - Somatic mutations
