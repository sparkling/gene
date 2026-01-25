---
id: herb
title: "HERB - High-throughput Experiment and Reference Database"
type: source
parent: ../README.md
tier: 1
status: active
category: traditional.medicine
subcategory: traditional.chinese.medicine
tags:
  - tcm
  - traditional-chinese-medicine
  - gene-expression
  - transcriptomics
  - herbs
---

# HERB - High-throughput Experiment and Reference Database

## Overview

HERB (High-throughput Experiment and Reference Database for TCM) is a unique database that integrates high-throughput gene expression data with Traditional Chinese Medicine information. Unlike other TCM databases that focus primarily on compound-target predictions, HERB provides experimentally-derived transcriptomic signatures for herbs, formulas, and ingredients.

The database connects TCM entities to their molecular effects through curated gene expression profiles from microarray and RNA-seq experiments. This enables systematic analysis of how TCM treatments affect cellular pathways and provides experimental validation for computational predictions from other databases.

HERB 2.0 significantly expanded coverage to include more herbs, experiments, and integration with disease signatures for connectivity analysis.

## Key Statistics

| Metric | Value |
|--------|-------|
| Herbs | 1,037 |
| Ingredients | 12,933 |
| Targets | 2,064 |
| Gene Expression Experiments | 2,000+ |
| Diseases | 866 |
| Herb-Ingredient Links | 49,258 |
| Ingredient-Target Links | 28,212 |

## Primary Use Cases

1. Finding gene expression signatures for TCM herbs and formulas
2. Connectivity analysis between TCM treatments and disease signatures
3. Experimental validation of predicted herb-target relationships
4. Pathway-level analysis of TCM mechanisms
5. Multi-omics integration for systems pharmacology

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Herb ID | HERB internal | HERB001234 |
| Ingredient ID | HERB internal | ING001234 |
| Target ID | Gene Symbol / Entrez | TP53 / 7157 |
| Experiment ID | GEO accession | GSE12345 |
| Disease ID | Internal / MESH | DIS001 / D001234 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://herb.ac.cn/ | Interactive browsing |
| Download | Bulk data available | TSV/Excel files |
| Search | Multi-entity search | Herbs, ingredients, genes |

## Data Model

```
Gene Expression Experiments (GEO)
            |
            v
+-----> Herbs (1,037) <-----> Diseases (866)
|           |
|           v
|    Ingredients (12,933)
|           |
|           v
+------ Targets (2,064) <---- Experimental Validation
```

## Expression Data Features

| Feature | Description |
|---------|-------------|
| Differential Expression | Up/down-regulated genes per treatment |
| Connectivity Scores | Similarity to disease signatures |
| Pathway Enrichment | KEGG/GO enrichment results |
| Dose-Response | Multiple concentration data when available |

## Integration Value

HERB provides unique experimental validation data that complements prediction-based databases:

| Pair With | Purpose |
|-----------|---------|
| BATMAN-TCM | Validate predicted targets with expression |
| SymMap | Add expression data to symptom mapping |
| DisGeNET | Connect herb effects to disease genes |

## Limitations

- Gene expression data biased toward commonly studied herbs
- Not all herbs have associated experiments
- Western experimental conditions may not reflect traditional usage

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
