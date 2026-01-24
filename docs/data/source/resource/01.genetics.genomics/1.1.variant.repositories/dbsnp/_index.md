---
id: dbsnp
title: "dbSNP"
type: data-source
category: genetics
subcategory: variant.repositories
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [variant, snp, polymorphism, ncbi, population]
---

# dbSNP

**Category:** [Genetics & Genomics](../../_index.md) > [Variant Repositories](../_index.md)

## Overview

dbSNP (Database of Single Nucleotide Polymorphisms) is NCBI's public archive for genetic variation data. It serves as the central repository for short genetic variations including single nucleotide variants (SNVs), small insertions/deletions, and microsatellites discovered across all organisms.

The database provides the canonical RS identifier system used universally in genetic research and clinical practice. dbSNP integrates with the NCBI Variation Services API, supporting modern variant notation standards including SPDI (Sequence-Position-Deletion-Insertion), HGVS, and GA4GH VR formats.

The ALFA (Allele Frequency Aggregator) component provides population-level allele frequencies from over 400,000 subjects across 12 populations, making it essential for variant filtering and population genetics research.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Variants | 1,000,000,000+ |
| ALFA Variants | 900M+ |
| ALFA Subjects | 400K+ |
| Populations | 12 |
| Reference Build | GRCh38.p14 |

## Primary Use Cases

1. Variant identification and RS number lookup
2. Population allele frequency queries
3. Variant notation conversion (SPDI, HGVS, VCF)
4. Clinical variant annotation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| RS ID | rs + integer | rs334 |
| RefSNP ID | Integer (no prefix) | 334 |
| SPDI | seq:pos:del:ins | NC_000011.10:5227001:T:A |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| REST API | https://api.ncbi.nlm.nih.gov/variation/v0 | JSON responses |
| FTP | ftp://ftp.ncbi.nlm.nih.gov/snp/ | VCF, frequency files |
| Web | https://www.ncbi.nlm.nih.gov/snp/ | Interactive search |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (NCBI) |
| Commercial Use | Yes |
| Rate Limit | 1 request/second recommended |

## See Also

- [Schema Documentation](./schema.md)
- [ClinVar](../clinvar/_index.md) - Clinical interpretations
- [gnomAD](../../1.3.population.genetics/gnomad/_index.md) - Population frequencies
