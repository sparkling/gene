---
id: spliceai
title: "SpliceAI"
type: data-source
category: genetics
subcategory: functional.prediction
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [splicing, prediction, deep-learning, variant, rna]
---

# SpliceAI

**Category:** [Genetics & Genomics](../../_index.md) > [Functional Prediction](../_index.md)

## Overview

SpliceAI is a deep learning-based tool developed by Illumina that predicts splicing alterations caused by single nucleotide variants. It uses a 32-layer deep neural network trained on over 10,000 human genes to predict the effect of variants on RNA splicing with high accuracy.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Position | chr-pos-ref-alt | 2-179446218-A-G |
| Gene | HGNC symbol | TTN |
| Delta Position | Integer offset | -2 to +2 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Lookup API | https://spliceailookup-api.broadinstitute.org | Pre-computed scores |
| GitHub | https://github.com/Illumina/SpliceAI | Source code |
| VEP Plugin | Ensembl VEP | Annotation |
| Pre-computed | Illumina downloads | VCF files |

## License

| Aspect | Value |
|--------|-------|
| License | GPLv3 (code); data varies |
| Commercial Use | Contact Illumina |
| Pre-computed | Freely available |

## See Also

- [Schema Documentation](./schema.md)
- [dbNSFP](../dbnsfp/_index.md) - Integrated scores
- [CADD](../cadd/_index.md) - Complementary scoring
