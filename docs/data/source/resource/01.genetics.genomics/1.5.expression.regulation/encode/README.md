---
id: encode
title: "ENCODE"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: expression.regulation
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - functional
  - regulatory
  - enhancer
  - chromatin
  - epigenomics
---

# ENCODE

ENCODE (Encyclopedia of DNA Elements) is a comprehensive project to identify all functional elements in the human and mouse genomes. It provides extensive data on transcription factor binding sites, chromatin accessibility, histone modifications, DNA methylation, and RNA expression across hundreds of cell types and tissues.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | ENCODE Consortium |
| **Website** | https://www.encodeproject.org/ |
| **Update Frequency** | Continuous |
| **Records** | 926,535+ human cCREs |
| **Latest Release** | ENCODE 4 (ongoing) |

The ENCODE Registry of candidate cis-Regulatory Elements (cCREs) provides a standardized set of regulatory elements including promoters, enhancers, and CTCF-binding sites, classified by chromatin signatures. The SCREEN database allows interactive exploration of these elements and their associated signals.

ENCODE data is essential for interpreting non-coding variants, understanding gene regulation, and identifying potential disease mechanisms. The project continues to generate new data types and expand tissue coverage.

## Key Statistics

| Metric | Value |
|--------|-------|
| Human cCREs | 926,535 |
| Mouse cCREs | 339,815 |
| Experiments | 20,000+ |
| Biosample Types | 500+ |
| Assay Types | 30+ |

## Primary Use Cases

1. Regulatory variant interpretation
2. Enhancer and promoter identification
3. Cell type-specific regulation analysis
4. Epigenomic data integration

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| cCRE | `EH38E + 7 digits` | `EH38E1234567` |
| Experiment | `ENCSR + 6 chars` | `ENCSR000AAA` |
| File | `ENCFF + 6 chars` | `ENCFF123ABC` |

## Limitations

- Cell line data may not reflect in vivo biology
- Tissue coverage incomplete (some tissues underrepresented)
- Functional validation limited for most elements
- Computational predictions may include false positives
- Large data volumes require significant storage

## Data Quality Notes

ENCODE implements standardized data processing pipelines with quality metrics for each assay type. Experiments are assigned audit flags indicating potential quality issues. The cCRE registry integrates multiple assays to reduce false positives, with classification based on chromatin accessibility and histone modification patterns.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [GTEx](../gtex/README.md) - Expression data
- [eQTLGen](../eqtlgen/README.md) - eQTL data
