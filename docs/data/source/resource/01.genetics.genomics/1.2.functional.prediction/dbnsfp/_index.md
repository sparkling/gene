---
id: dbnsfp
title: "dbNSFP"
type: data-source
category: genetics
subcategory: functional.prediction
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [pathogenicity, prediction, annotation, ensemble, conservation]
---

# dbNSFP

**Category:** [Genetics & Genomics](../../_index.md) > [Functional Prediction](../_index.md)

## Overview

dbNSFP (database for Non-synonymous SNPs' Functional Predictions) is a comprehensive database integrating functional predictions from 36+ algorithms for all possible non-synonymous SNVs in the human genome. It serves as a one-stop resource combining deleteriousness scores, conservation metrics, and population frequencies with gene-level annotations.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Position | chr:pos(1-based) | 17:7673803 |
| Gene | HGNC symbol | TP53 |
| Ensembl | ENSG/ENST IDs | ENSG00000141510 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Amazon S3 | dbnsfp.s3.amazonaws.com | Primary download |
| Website | https://www.dbnsfp.org/ | Documentation |
| VEP Plugin | Ensembl VEP | Annotation |
| SnpSift | SnpEff suite | Annotation |

## License

| Aspect | Value |
|--------|-------|
| License | Academic free; commercial contact |
| Commercial Use | Requires contact (some scores excluded) |
| Citation | Required (version-specific) |

## See Also

- [Schema Documentation](./schema.md)
- [AlphaMissense](../alphamissense/_index.md) - Included predictor
- [CADD](../cadd/_index.md) - Included predictor
