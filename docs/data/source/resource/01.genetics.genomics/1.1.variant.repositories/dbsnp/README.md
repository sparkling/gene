---
id: dbsnp
title: "dbSNP"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: variant.repositories
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - variant
  - snp
  - polymorphism
  - ncbi
  - population
---

# dbSNP

dbSNP (Database of Single Nucleotide Polymorphisms) is NCBI's public archive for genetic variation data. It serves as the central repository for short genetic variations including single nucleotide variants (SNVs), small insertions/deletions, and microsatellites discovered across all organisms.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | NCBI/NIH |
| **Website** | https://www.ncbi.nlm.nih.gov/snp/ |
| **Update Frequency** | Continuous |
| **Records** | 1,000,000,000+ variants |
| **Latest Release** | Build 156 (GRCh38.p14) |

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

| Identifier | Format | Example |
|------------|--------|---------|
| RS ID | `rs + integer` | `rs334` |
| RefSNP ID | `Integer (no prefix)` | `334` |
| SPDI | `seq:pos:del:ins` | `NC_000011.10:5227001:T:A` |

## Limitations

- Limited clinical interpretation (use ClinVar for pathogenicity)
- Some legacy RS IDs have been merged or withdrawn
- Population frequency data limited to ALFA-included studies
- Rate limits apply to API access (1 req/sec recommended)
- Complex variants may have multiple RS IDs

## Data Quality Notes

dbSNP undergoes continuous quality control including validation status tracking, allele frequency verification, and mapping accuracy assessment. Variants are flagged for clinical significance when available from ClinVar. Users should verify variant mapping on current reference assemblies as coordinates may shift between builds.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [ClinVar](../clinvar/README.md) - Clinical interpretations
- [gnomAD](../../1.3.population.genetics/gnomad/README.md) - Population frequencies
