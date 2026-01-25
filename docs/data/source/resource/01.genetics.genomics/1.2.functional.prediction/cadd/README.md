---
id: cadd
title: "CADD"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: functional.prediction
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - pathogenicity
  - deleteriousness
  - annotation
  - scoring
  - prediction
---

# CADD

CADD (Combined Annotation Dependent Depletion) is an integrative annotation method that scores the deleteriousness of single nucleotide variants and insertion/deletions in the human genome. It combines diverse annotations into a single score by contrasting variants that survived natural selection with simulated mutations.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | University of Washington / Berlin Institute of Health |
| **Website** | https://cadd.gs.washington.edu/ |
| **Update Frequency** | Periodic releases |
| **Records** | ~9 billion possible SNVs |
| **Latest Release** | v1.7 |

CADD integrates over 60 genomic features including conservation scores, regulatory information, protein-level annotations, and transcript information. The resulting C-score is PHRED-scaled, where higher scores indicate more likely deleteriousness. A score of 20 indicates the variant is in the top 1% of deleterious variants genome-wide.

The tool provides pre-computed scores for all possible SNVs in the reference genome, as well as scoring services for novel variants and indels. CADD scores are widely used in clinical genetics and research for variant prioritization.

## Key Statistics

| Metric | Value |
|--------|-------|
| SNV Coverage | All possible (~9 billion) |
| Indel Support | Yes (via scoring) |
| Features Integrated | 60+ |
| Current Version | v1.7 |
| Reference Build | GRCh38, GRCh37 |

## Primary Use Cases

1. Variant deleteriousness scoring
2. Clinical variant prioritization
3. Non-coding variant assessment
4. Indel pathogenicity prediction

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Genomic | `chr-pos-ref-alt` | `1-12345-A-G` |
| CADD Score | `PHRED-scaled` | `25.3` |
| Raw Score | `Continuous` | `4.52` |

## Limitations

- Commercial use requires license agreement
- Not specific to disease or variant type
- May score common benign variants as deleterious
- Indel scoring less validated than SNV scoring
- Does not distinguish between different disease mechanisms

## Data Quality Notes

CADD PHRED scores are calibrated such that a score of 10 indicates top 10% deleteriousness, 20 indicates top 1%, and 30 indicates top 0.1%. Scores above 15-20 are commonly used as pathogenicity thresholds, though optimal cutoffs vary by application. CADD performs best for coding variants and is less validated for regulatory regions.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [dbNSFP](../dbnsfp/README.md) - Integrated CADD scores
- [AlphaMissense](../alphamissense/README.md) - Missense predictor
- [SpliceAI](../spliceai/README.md) - Splice prediction
