---
id: dbvar
title: "dbVar"
type: data-source
category: genetics
subcategory: variant.repositories
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [variant, structural, cnv, ncbi, deletion, duplication]
---

# dbVar

**Category:** [Genetics & Genomics](../../_index.md) > [Variant Repositories](../_index.md)

## Overview

dbVar is NCBI's public database of human genomic structural variation, archiving variants greater than 50 base pairs including deletions, duplications, insertions, inversions, translocations, and complex rearrangements. It serves as the primary repository for structural variant (SV) data from research and clinical studies.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Variant Region | nsv/esv/dsv + integer | nsv1234567 |
| Variant Call | nssv/essv/dssv + integer | nssv14580340 |
| Study | nstd/estd/dstd + integer | nstd166 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| FTP | ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/ | TSV, BED, VCF, GVF |
| Web | https://www.ncbi.nlm.nih.gov/dbvar | Interactive search |
| GitHub | https://github.com/ncbi/dbvar | Documentation, tools |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [ClinVar](../clinvar/_index.md) - Clinical SVs forwarded here
- [gnomAD](../../1.3.population.genetics/gnomad/_index.md) - SV frequencies
