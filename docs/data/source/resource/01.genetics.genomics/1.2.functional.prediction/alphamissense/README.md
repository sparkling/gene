---
id: alphamissense
title: "AlphaMissense"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: functional.prediction
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - pathogenicity
  - missense
  - deep-learning
  - alphafold
  - prediction
---

# AlphaMissense

AlphaMissense is Google DeepMind's AI-powered pathogenicity predictor for all possible human missense variants. Built on AlphaFold2's protein structure prediction architecture, it incorporates evolutionary conservation, structural context, and population frequency data to generate calibrated pathogenicity scores.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Google DeepMind |
| **Website** | https://alphamissense.hegelab.org/ |
| **Update Frequency** | Static (model v1) |
| **Records** | 71,000,000+ missense variants |
| **Latest Release** | September 2023 |

Pre-computed predictions cover 71 million single nucleotide missense variants across 19,233 canonical transcripts, achieving state-of-the-art performance with 89% of ClinVar pathogenic variants correctly classified as likely_pathogenic and 90% of benign variants classified as likely_benign.

The model provides continuous pathogenicity scores from 0.0 (benign) to 1.0 (pathogenic) with three-class classifications based on validated thresholds. AlphaMissense scores are integrated into dbNSFP and available through Ensembl VEP plugins.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Variants | 71,000,000+ |
| Canonical Transcripts | 19,233 |
| AA Substitutions | 216,000,000+ |
| AUROC (ClinVar) | 0.940 |
| File Size | 643 MB (hg38 compressed) |

## Primary Use Cases

1. Missense variant pathogenicity prediction
2. Clinical variant prioritization
3. Research variant filtering
4. Gene-level pathogenicity assessment

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Genomic | `chr:pos:ref:alt` | `chr7:140753336:A:T` |
| Protein | `UniProt + AA change` | `P15056 V600E` |
| Transcript | `ENST + version` | `ENST00000288602.11` |

## Limitations

- Limited to missense variants only (no indels, nonsense, splice)
- Model trained primarily on known pathogenic variants (potential bias)
- Predictions are per-transcript; isoform effects may differ
- Model weights not publicly released (pre-computed scores only)
- Novel protein folds may have reduced accuracy

## Data Quality Notes

AlphaMissense achieves AUROC of 0.940 on ClinVar benchmark variants and outperforms existing methods including REVEL and CADD for missense classification. Scores are calibrated such that thresholds of 0.34 (benign) and 0.564 (pathogenic) achieve 90% precision for respective classes. Users should consider combining with other evidence types for clinical interpretation.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [dbNSFP](../dbnsfp/README.md) - Integrated scores
- [CADD](../cadd/README.md) - Alternative predictor
