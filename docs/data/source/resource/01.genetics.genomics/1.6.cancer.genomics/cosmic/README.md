---
id: cosmic
title: "COSMIC"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: cancer.genomics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - cancer
  - somatic
  - mutation
  - signature
  - census
---

# COSMIC

COSMIC (Catalogue Of Somatic Mutations In Cancer) is the world's largest and most comprehensive resource for exploring somatic mutations in human cancer. Developed by the Wellcome Sanger Institute, it provides expert-curated data on somatic mutations, gene fusions, copy number variants, and mutational signatures across all cancer types.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Wellcome Sanger Institute |
| **Website** | https://cancer.sanger.ac.uk/cosmic |
| **Update Frequency** | Quarterly releases |
| **Records** | 17,000,000+ mutations |
| **Latest Release** | v99 |

The database includes the Cancer Gene Census, a catalog of genes with mutations causally implicated in cancer, and mutational signatures representing the patterns of mutations caused by specific mutagenic processes. COSMIC integrates data from literature curation, large-scale screens, and major cancer genome projects.

COSMIC v99+ provides genomic coordinates for all mutations on GRCh38, with legacy support for GRCh37. The data is essential for cancer research, clinical interpretation of tumor sequencing, and understanding cancer biology.

## Key Statistics

| Metric | Value |
|--------|-------|
| Coding Mutations | 17,000,000+ |
| Samples | 1,500,000+ |
| Tumor Types | 500+ |
| Cancer Genes (Census) | 736 |
| Mutational Signatures | 100+ |

## Primary Use Cases

1. Somatic variant annotation
2. Driver gene identification
3. Mutational signature analysis
4. Cancer type comparison

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| COSMIC ID | `COSM/COSV + integer` | `COSM476` |
| Gene | `HUGO symbol` | `BRAF` |
| Signature | `SBS/DBS/ID + number` | `SBS1` |

## Limitations

- Registration required for data access
- Commercial use requires paid license
- Literature curation may lag behind primary research
- Sample overlap exists between studies
- Not all mutations have functional characterization

## Data Quality Notes

COSMIC data is curated from peer-reviewed publications and major cancer sequencing projects (TCGA, ICGC). Variants are classified as confirmed somatic (from matched tumor-normal sequencing) or reported in cancer (from tumor-only studies). The COSMIC ID system provides stable identifiers for tracking variants across releases.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [Cancer Gene Census](../cancer.gene.census/README.md) - Driver genes
- [cBioPortal](../cbioportal/README.md) - Study data
- [OncoKB](../oncokb/README.md) - Clinical annotations
