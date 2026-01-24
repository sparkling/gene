---
id: cosmic
title: "COSMIC"
type: data-source
category: genetics
subcategory: cancer.genomics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [cancer, somatic, mutation, signature, census]
---

# COSMIC

**Category:** [Genetics & Genomics](../../_index.md) > [Cancer Genomics](../_index.md)

## Overview

COSMIC (Catalogue Of Somatic Mutations In Cancer) is the world's largest and most comprehensive resource for exploring somatic mutations in human cancer. Developed by the Wellcome Sanger Institute, it provides expert-curated data on somatic mutations, gene fusions, copy number variants, and mutational signatures across all cancer types.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| COSMIC ID | COSM/COSV + integer | COSM476 |
| Gene | HUGO symbol | BRAF |
| Signature | SBS/DBS/ID + number | SBS1 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://cancer.sanger.ac.uk/cosmic | Interactive search |
| Downloads | Registered access | VCF, TSV files |
| API | REST endpoints | Programmatic access |
| Genome Browser | GRCh38 tracks | Visualization |

## License

| Aspect | Value |
|--------|-------|
| License | Free for academic; commercial license |
| Commercial Use | Requires license |
| Registration | Required for downloads |

## See Also

- [Cancer Gene Census](../cancer.gene.census/_index.md) - Driver genes
- [cBioPortal](../cbioportal/_index.md) - Study data
- [OncoKB](../oncokb/_index.md) - Clinical annotations
