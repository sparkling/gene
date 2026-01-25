---
id: dbvar
title: "dbVar"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: variant.repositories
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - variant
  - structural
  - cnv
  - ncbi
  - deletion
  - duplication
---

# dbVar

dbVar is NCBI's public database of human genomic structural variation, archiving variants greater than 50 base pairs including deletions, duplications, insertions, inversions, translocations, and complex rearrangements. It serves as the primary repository for structural variant (SV) data from research and clinical studies.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | NCBI/NIH |
| **Website** | https://www.ncbi.nlm.nih.gov/dbvar |
| **Update Frequency** | Continuous |
| **Records** | 6,000,000+ structural variants |
| **Latest Release** | Current (continuous updates) |

The database uses a hierarchical identifier system: studies (nstd/estd/dstd), variant regions (nsv/esv/dsv), and variant calls (nssv/essv/dssv). This structure captures both the genomic locations of structural variants and the individual experimental observations supporting them.

dbVar synchronizes with the European Variation Archive (EVA/DGVa) through monthly data exchanges, ensuring global coverage of structural variation data. Non-redundant files provide deduplicated variant sets for common analysis workflows.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total SVs | 6,000,000+ |
| Deletions | 1.9M (non-redundant) |
| Duplications | 659K (non-redundant) |
| Insertions | 1.7M (non-redundant) |
| Studies | 185+ |

## Primary Use Cases

1. Structural variant annotation and lookup
2. Copy number variation analysis
3. Clinical CNV interpretation
4. Population SV frequency assessment

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Variant Region | `nsv/esv/dsv + integer` | `nsv1234567` |
| Variant Call | `nssv/essv/dssv + integer` | `nssv14580340` |
| Study | `nstd/estd/dstd + integer` | `nstd166` |

## Limitations

- Breakpoint precision varies by detection method
- Limited clinical interpretation for most SVs
- Population frequencies less comprehensive than SNV databases
- Complex rearrangements may be incompletely characterized
- Different studies use varying SV calling thresholds

## Data Quality Notes

dbVar applies quality filters based on submission validation status and supporting evidence. Variant regions represent consensus boundaries across multiple observations, while variant calls preserve original study-specific coordinates. Users should consider validation status and supporting read depth when interpreting SV calls.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [ClinVar](../clinvar/README.md) - Clinical SVs forwarded here
- [gnomAD](../../1.3.population.genetics/gnomad/README.md) - SV frequencies
