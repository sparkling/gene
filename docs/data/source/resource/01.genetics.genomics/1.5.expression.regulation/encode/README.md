---
id: encode
title: "ENCODE"
type: data-source
category: genetics
subcategory: expression.regulation
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [functional, regulatory, enhancer, chromatin, epigenomics]
---

# ENCODE

**Category:** [Genetics & Genomics](../../_index.md) > [Expression & Regulation](../_index.md)

## Overview

ENCODE (Encyclopedia of DNA Elements) is a comprehensive project to identify all functional elements in the human and mouse genomes. It provides extensive data on transcription factor binding sites, chromatin accessibility, histone modifications, DNA methylation, and RNA expression across hundreds of cell types and tissues.

The ENCODE Registry of candidate cis-Regulatory Elements (cCREs) provides a standardized set of regulatory elements including promoters, enhancers, and CTCF-binding sites, classified by chromatin signatures. The SCREEN database allows interactive exploration of these elements and their associated signals.

ENCODE data is essential for interpreting non-coding variants, understanding gene regulation, and identifying potential disease mechanisms. The project continues to generate new data types and expand tissue coverage.

## Key Statistics

| Metric | Value |
|--------|-------|
| Human cCREs | 926,535 |
| Mouse cCREs | 339,815 |
| Experiments | 20,000+ |
| Biosample Types | 500+ |
| Assay Types | 30+ |

## Primary Use Cases

1. Regulatory variant interpretation
2. Enhancer and promoter identification
3. Cell type-specific regulation analysis
4. Epigenomic data integration

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| cCRE | EH38E + 7 digits | EH38E1234567 |
| Experiment | ENCSR + 6 chars | ENCSR000AAA |
| File | ENCFF + 6 chars | ENCFF123ABC |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Portal | https://www.encodeproject.org/ | Interactive search |
| SCREEN | https://screen.encodeproject.org/ | cCRE browser |
| REST API | /api/ endpoints | Programmatic access |
| Cloud | AWS, Google Cloud | Data downloads |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access |
| Commercial Use | Yes |
| Citation | Required per data use |

## See Also

- [Schema Documentation](./schema.md)
- [GTEx](../gtex/_index.md) - Expression data
- [eQTLGen](../eqtlgen/_index.md) - eQTL data
