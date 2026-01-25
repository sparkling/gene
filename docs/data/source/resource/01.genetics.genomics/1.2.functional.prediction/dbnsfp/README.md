---
id: dbnsfp
title: "dbNSFP"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: functional.prediction
tier: 2
status: active
last_updated: 2026-01-25
tags:
  - pathogenicity
  - prediction
  - annotation
  - ensemble
  - conservation
---

# dbNSFP

dbNSFP (database for Non-synonymous SNPs' Functional Predictions) is a comprehensive database integrating functional predictions from 36+ algorithms for all possible non-synonymous SNVs in the human genome. It serves as a one-stop resource combining deleteriousness scores, conservation metrics, and population frequencies with gene-level annotations.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Liu Lab (USF) |
| **Website** | https://www.dbnsfp.org/ |
| **Update Frequency** | Periodic releases |
| **Records** | 85,000,000+ variants |
| **Latest Release** | v5.3.1 |

The database covers approximately 85 million variants including non-synonymous SNVs and splice-site variants. It integrates predictions from SIFT, PolyPhen-2, CADD, REVEL, AlphaMissense, and many other algorithms, along with conservation scores (PhyloP, GERP++, phastCons) and population frequencies from gnomAD, 1000 Genomes, and TOPMed.

dbNSFP provides normalized rankscores for cross-algorithm comparison and includes gene-level annotations such as pLI, LOEUF, and disease associations from OMIM and Orphanet.

## Key Statistics

| Metric | Value |
|--------|-------|
| Non-synonymous SNVs | 83,049,507 |
| Splice-site SNVs | 2,446,464 |
| Prediction Algorithms | 36+ |
| Conservation Scores | 9 |
| Current Version | 5.3.1 |

## Primary Use Cases

1. Multi-algorithm pathogenicity assessment
2. Clinical variant annotation
3. Variant prioritization pipelines
4. Consensus-based interpretation

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Position | `chr:pos (1-based)` | `17:7673803` |
| Gene | `HGNC symbol` | `TP53` |
| Ensembl | `ENSG/ENST IDs` | `ENSG00000141510` |

## Limitations

- File size is very large (~40GB uncompressed)
- Some component scores have commercial licensing restrictions
- Scores may disagree; no single consensus metric
- Limited to coding and splice-site variants
- Updates lag behind component databases

## Data Quality Notes

dbNSFP provides rankscores (0-1) normalized across all variants for each prediction algorithm, enabling fair comparison between methods. The MetaSVM and MetaLR ensemble scores combine multiple predictors with validated performance. Users should be aware that some component scores (e.g., CADD, REVEL) may have separate licensing requirements for commercial use.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [AlphaMissense](../alphamissense/README.md) - Included predictor
- [CADD](../cadd/README.md) - Included predictor
