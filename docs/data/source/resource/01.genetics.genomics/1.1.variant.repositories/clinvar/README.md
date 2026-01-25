---
id: clinvar
title: "ClinVar"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: variant.repositories
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - variant
  - clinical
  - pathogenicity
  - ncbi
  - germline
  - somatic
---

# ClinVar

ClinVar is NCBI's public archive of interpretations of clinically relevant variants. It aggregates information about genomic variation and its relationship to human health, collecting submissions from clinical laboratories, research groups, and expert panels worldwide.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | NCBI/NIH |
| **Website** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Update Frequency** | Weekly (Mondays) |
| **Records** | 10,000,000+ submissions |
| **Latest Release** | Current (weekly updates) |

The database uses a three-tier accession system: VCV (Variation Archive) for aggregate interpretations across all submitters, RCV (Record) for individual variant-condition pairs, and SCV (Submitted Record) for individual submissions. This structure enables tracking of both consensus interpretations and the underlying evidence from multiple sources.

ClinVar integrates with other NCBI resources including dbSNP, dbVar, and MedGen, providing comprehensive variant annotation including genomic coordinates, HGVS nomenclature, gene associations, and disease relationships with standardized phenotype terms.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Submissions | 10,000,000+ |
| Unique Variants | 1,500,000+ |
| Top Submitter | LabCorp Genetics (1.89M) |
| Reference Genome | GRCh38, GRCh37 |
| Update Frequency | Weekly (Mondays) |

## Primary Use Cases

1. Clinical variant interpretation and pathogenicity assessment
2. Diagnostic laboratory variant classification lookup
3. Research variant prioritization and filtering
4. Gene-disease association discovery

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| VCV | `VCV + 9 digits + version` | `VCV000000123.4` |
| RCV | `RCV + 9 digits + version` | `RCV000012345.6` |
| SCV | `SCV + 9 digits + version` | `SCV000123456.7` |
| AlleleID | `Integer` | `12345` |
| VariationID | `Integer` | `67890` |

## Limitations

- Submitter interpretations may conflict; review status indicates consensus level
- Classification criteria vary between submitters and over time
- Some variants lack functional evidence (computational predictions only)
- Updates may lag behind current literature
- Limited representation of rare populations in frequency data

## Data Quality Notes

ClinVar uses a star-rating system (0-4 stars) to indicate review status and evidence quality. Four-star ratings indicate expert panel review or practice guidelines, while lower ratings reflect varying levels of submitter consensus. Users should prioritize variants with multiple concordant submissions and higher star ratings for clinical decision-making.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [dbSNP](../dbsnp/README.md) - SNP identifiers
- [dbVar](../dbvar/README.md) - Structural variants
