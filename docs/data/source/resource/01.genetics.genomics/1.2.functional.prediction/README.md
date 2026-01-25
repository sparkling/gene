---
id: functional.prediction
title: "Functional Prediction"
type: subcategory
parent: ../README.md
last_updated: 2026-01-23
status: active
tags: [prediction, pathogenicity, splicing, missense, machine-learning]
---

# Functional Prediction

**Parent:** [Genetics & Genomics](../README.md)

## Overview

Functional prediction resources provide computational assessments of variant effects on protein function, splicing, and disease risk. These tools use evolutionary conservation, protein structure, and machine learning to predict variant pathogenicity.

These scores are essential for prioritizing variants of uncertain significance and for research-grade variant filtering. Different predictors are optimized for different variant types (missense, splice site, indels).

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [CADD](./cadd/README.md) | 1 | Combined annotation-dependent depletion scores |
| [dbNSFP](./dbnsfp/README.md) | 1 | Aggregated functional prediction database |
| [SpliceAI](./spliceai/README.md) | 1 | Deep learning splice variant predictions |
| [AlphaMissense](./alphamissense/README.md) | 1 | AlphaFold-based missense predictions |

## Integration Notes

dbNSFP aggregates 30+ prediction scores including CADD, REVEL, and others. Use CADD for genome-wide variant prioritization. SpliceAI specifically for splice effects. AlphaMissense for structure-aware missense predictions. Combine multiple predictors for consensus scoring.
