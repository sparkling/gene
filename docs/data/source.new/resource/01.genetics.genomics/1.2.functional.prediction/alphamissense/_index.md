---
id: alphamissense
title: "AlphaMissense"
type: data-source
category: genetics
subcategory: functional.prediction
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [pathogenicity, missense, deep-learning, alphafold, prediction]
---

# AlphaMissense

**Category:** [Genetics & Genomics](../../_index.md) > [Functional Prediction](../_index.md)

## Overview

AlphaMissense is Google DeepMind's AI-powered pathogenicity predictor for all possible human missense variants. Built on AlphaFold2's protein structure prediction architecture, it incorporates evolutionary conservation, structural context, and population frequency data to generate calibrated pathogenicity scores.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Genomic | chr:pos:ref:alt | chr7:140753336:A:T |
| Protein | UniProt + AA change | P15056 V600E |
| Transcript | ENST + version | ENST00000288602.11 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Google Cloud | gs://dm_alphamissense/ | Official source |
| Zenodo | https://zenodo.org/records/8360242 | Alternative |
| VEP Plugin | Ensembl VEP | Annotation integration |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 (Google Cloud) |
| Commercial Use | Yes (with attribution) |
| Model Weights | Not released |

## See Also

- [Schema Documentation](./schema.md)
- [dbNSFP](../dbnsfp/_index.md) - Integrated scores
- [CADD](../cadd/_index.md) - Alternative predictor
