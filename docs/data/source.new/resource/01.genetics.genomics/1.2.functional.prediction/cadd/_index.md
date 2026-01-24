---
id: cadd
title: "CADD"
type: data-source
category: genetics
subcategory: functional.prediction
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [pathogenicity, deleteriousness, annotation, scoring, prediction]
---

# CADD

**Category:** [Genetics & Genomics](../../_index.md) > [Functional Prediction](../_index.md)

## Overview

CADD (Combined Annotation Dependent Depletion) is an integrative annotation method that scores the deleteriousness of single nucleotide variants and insertion/deletions in the human genome. It combines diverse annotations into a single score by contrasting variants that survived natural selection with simulated mutations.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Genomic | chr-pos-ref-alt | 1-12345-A-G |
| CADD Score | PHRED-scaled | 25.3 |
| Raw Score | Continuous | 4.52 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://cadd.gs.washington.edu/ | Online scoring |
| Download | https://cadd.gs.washington.edu/download | Pre-computed scores |
| API | REST endpoint | Variant scoring |

## License

| Aspect | Value |
|--------|-------|
| License | Academic free; commercial contact |
| Commercial Use | Requires license |
| Citation | Required |

## See Also

- [dbNSFP](../dbnsfp/_index.md) - Integrated CADD scores
- [AlphaMissense](../alphamissense/_index.md) - Missense predictor
- [SpliceAI](../spliceai/_index.md) - Splice prediction
