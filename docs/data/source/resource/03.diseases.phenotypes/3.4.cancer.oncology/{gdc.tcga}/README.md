---
id: gdc.tcga
title: "GDC/TCGA - Genomic Data Commons and The Cancer Genome Atlas"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: cancer.oncology
tags:
  - cancer
  - genomics
  - tcga
  - nci
  - multi-omics
  - somatic-mutations
---

# GDC/TCGA - Genomic Data Commons and The Cancer Genome Atlas

## Overview

The Genomic Data Commons (GDC) is the National Cancer Institute's centralized repository for cancer genomic data, including the landmark datasets from The Cancer Genome Atlas (TCGA) program. Together, GDC and TCGA represent the most comprehensive resource for cancer molecular characterization, enabling discovery of cancer driver genes, molecular subtypes, and therapeutic targets.

TCGA, which ran from 2006 to 2017, comprehensively characterized over 20,000 primary cancer samples across 33 cancer types using genomic, epigenomic, transcriptomic, and proteomic profiling. The project identified thousands of cancer-associated mutations, characterized molecular subtypes with clinical significance, and established benchmarks for cancer genomics research.

The GDC continues to host TCGA data alongside newer cancer genomics projects, providing harmonized data processing pipelines and a unified data model. The platform supports controlled-access for protected genomic data while making summary mutations, copy number alterations, and gene expression data publicly available.

## Key Statistics

| Metric | Value |
|--------|-------|
| TCGA Samples | 20,000+ |
| Cancer Types (TCGA) | 33 |
| Total GDC Cases | 80,000+ |
| Somatic Mutations | 3M+ |
| Projects in GDC | 70+ |

## Primary Use Cases

1. Pan-cancer analysis of driver mutations and pathways
2. Molecular subtype discovery and validation
3. Biomarker identification for precision oncology
4. Drug target discovery in cancer
5. Multi-omics integration for cancer research

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| TCGA Barcode | `TCGA-[A-Z]{2}-[A-Z0-9]{4}` | TCGA-A1-A0SK |
| GDC UUID | UUID format | 8a3a3f26-... |
| Case ID | `[A-Z]+-[A-Z0-9]+` | TCGA-A1-A0SK |
| Project ID | `[A-Z]+-[A-Z]+` | TCGA-BRCA |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| GDC Portal | https://portal.gdc.cancer.gov/ | Web interface |
| GDC API | https://api.gdc.cancer.gov/ | REST API |
| GDC Data Transfer | Command line tool | Bulk downloads |
| cBioPortal | https://cbioportal.org/ | Analysis interface |

## Data Types

| Data Type | Description | Access |
|-----------|-------------|--------|
| Somatic Mutations | MAF files | Open |
| Copy Number | Gene-level CNV | Open |
| Gene Expression | RNA-seq FPKM | Open |
| Clinical Data | Demographics, outcomes | Open |
| Aligned Reads | BAM files | Controlled |
| Raw Sequences | FASTQ files | Controlled |

## Cancer Types (TCGA)

| Abbreviation | Cancer Type | Samples |
|--------------|-------------|---------|
| BRCA | Breast invasive carcinoma | 1,098 |
| LUAD | Lung adenocarcinoma | 585 |
| PRAD | Prostate adenocarcinoma | 500 |
| COAD | Colon adenocarcinoma | 478 |
| KIRC | Kidney clear cell carcinoma | 537 |

## Limitations

- Individual-level data requires dbGaP controlled access approval
- TCGA samples are treatment-naive primary tumors; recurrence data limited
- Sample demographics reflect US populations; global diversity limited
- Data processing pipelines may differ from other cancer databases

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
