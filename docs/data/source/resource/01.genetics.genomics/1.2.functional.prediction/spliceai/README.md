---
id: spliceai
title: "SpliceAI"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: functional.prediction
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - splicing
  - prediction
  - deep-learning
  - variant
  - rna
---

# SpliceAI

SpliceAI is a deep learning-based tool developed by Illumina that predicts splicing alterations caused by single nucleotide variants. It uses a 32-layer deep neural network trained on over 10,000 human genes to predict the effect of variants on RNA splicing with high accuracy.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Illumina |
| **Website** | https://spliceailookup.broadinstitute.org/ |
| **Update Frequency** | Static (pre-computed) |
| **Records** | All SNVs in genes |
| **Latest Release** | v1.3 |

The model generates four delta scores representing the probability of splice site gain or loss at acceptor and donor positions: delta_score (DS) for acceptor gain, acceptor loss, donor gain, and donor loss. Scores range from 0 to 1, with higher values indicating stronger predicted splicing effects.

SpliceAI pre-computed scores are available for all possible SNVs within genes, and the model can score novel variants including indels up to 10kb from splice sites. It has become a standard tool in clinical genetics for evaluating potential splicing effects of variants.

## Key Statistics

| Metric | Value |
|--------|-------|
| Coverage | All SNVs in genes |
| Training Genes | 10,000+ |
| Network Depth | 32 layers |
| Context Window | 10kb each direction |
| Delta Score Range | 0.0 - 1.0 |

## Primary Use Cases

1. Splice variant effect prediction
2. Cryptic splice site identification
3. Clinical variant interpretation
4. Intronic variant assessment

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Position | `chr-pos-ref-alt` | `2-179446218-A-G` |
| Gene | `HGNC symbol` | `TTN` |
| Delta Position | `Integer offset` | `-2 to +2` |

## Limitations

- Predictions are sequence-based only (no tissue-specific splicing)
- May miss complex splicing effects requiring cellular context
- Indel predictions less validated than SNV predictions
- Does not predict extent of exon skipping or inclusion
- Commercial use requires Illumina licensing

## Data Quality Notes

SpliceAI achieves >95% accuracy for detecting known splice-altering variants in validation datasets. Delta scores above 0.2 are considered potentially splice-altering, with scores above 0.5 indicating high confidence. Predictions should be validated with RNA-seq or RT-PCR when available, particularly for clinical applications.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [dbNSFP](../dbnsfp/README.md) - Integrated scores
- [CADD](../cadd/README.md) - Complementary scoring
