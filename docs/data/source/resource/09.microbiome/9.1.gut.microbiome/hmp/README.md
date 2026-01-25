---
id: hmp
title: "HMP - Human Microbiome Project"
type: data-source
category: microbiome
subcategory: gut.microbiome
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [microbiome, metagenomics, reference, multi-omics, NIH]
---

# HMP - Human Microbiome Project

**Category:** [Microbiome](../../../README.md) > [Gut Microbiome](../README.md)

## Overview

The Human Microbiome Project (HMP) and integrative HMP (iHMP) characterize the human microbiome and its role in health and disease. This NIH-funded initiative provides comprehensive multi-omic data including 16S rRNA sequencing, whole-genome shotgun (WGS) sequencing, transcriptomics, proteomics, metabolomics, and lipidome data.

The project encompasses multiple major studies: HMP1 (healthy cohort baseline), IBDMDB (inflammatory bowel disease), T2D (prediabetes/type 2 diabetes), and MOMS-PI (pregnancy and preterm birth).

HMP data is essential for microbiome reference, disease association studies, and understanding the role of microbial communities in human health.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Data Volume | 48+ TB |
| Studies | 4 major studies |
| Data Types | 16S, WGS, RNA-seq, Proteomics, Metabolomics |
| File Formats | FASTA, FASTQ, SFF, BIOM, TSV |
| Cloud Availability | AWS Open Data Registry |

## Primary Use Cases

1. Reference microbiome characterization
2. Disease-microbiome associations
3. Multi-omic integration studies
4. Longitudinal microbiome analysis
5. Method benchmarking

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Sample ID | `HMP_[A-Z0-9]+` | HMP_2012SN-1234 |
| Subject ID | `[0-9]{3}-[0-9]{6}` | 159-494937 |
| Visit Number | Numeric | 1, 2, 3 |
| Body Site | Controlled vocabulary | Stool, Anterior_nares |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Data Portal | https://portal.hmpdacc.org | Main data access |
| AWS Open Data | s3://human-microbiome-project | Cloud access |
| GitHub Schemas | https://github.com/ihmpdcc/osdf-schemas | OSDF schemas |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access (NIH policy) |
| Commercial Use | Permitted |
| Attribution | Required per NIH |

## See Also

- [Schema Documentation](./schema.md)
- [GMrepo](../gmrepo/README.md) - Curated gut microbiome data
- [MetaHIT](../metahit/README.md) - European reference metagenomes
