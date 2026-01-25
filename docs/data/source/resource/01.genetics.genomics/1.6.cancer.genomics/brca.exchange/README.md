---
id: brca-exchange
title: "BRCA Exchange"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: cancer.genomics
tier: 2
status: active
last_updated: 2026-01-25
tags:
  - cancer
  - brca1
  - brca2
  - hereditary
  - variant
---

# BRCA Exchange

BRCA Exchange is a comprehensive resource for BRCA1 and BRCA2 variant data aggregation and classification. As part of the Global Alliance for Genomics and Health (GA4GH), it brings together data from clinical laboratories, research databases, and population resources to provide expert-reviewed classifications for variants in these critical cancer predisposition genes.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | BRCA Exchange / GA4GH |
| **Website** | https://brcaexchange.org/ |
| **Update Frequency** | Monthly |
| **Records** | 52,000+ BRCA variants |
| **Latest Release** | Current (monthly updates) |

The database aggregates variant classifications from ClinVar, LOVD, ENIGMA consortium, and other sources, providing a unified view of clinical significance with expert review. The ENIGMA consortium provides the most authoritative classifications based on extensive functional and clinical evidence.

BRCA Exchange is essential for clinical laboratories interpreting BRCA variants, genetic counselors, and researchers studying hereditary breast and ovarian cancer.

## Key Statistics

| Metric | Value |
|--------|-------|
| BRCA1 Variants | 22,000+ |
| BRCA2 Variants | 30,000+ |
| Expert Reviewed | 15,000+ |
| Data Sources | 10+ |
| ENIGMA Classified | 5,000+ |

## Primary Use Cases

1. BRCA variant classification lookup
2. Clinical variant interpretation
3. Hereditary cancer risk assessment
4. Research variant annotation

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Gene | `BRCA1/BRCA2` | `BRCA1` |
| HGVS | `Protein/DNA` | `c.68_69delAG` |
| Source | `Database name` | `ClinVar` |

## Limitations

- Limited to BRCA1 and BRCA2 only
- Classification discordances exist between sources
- Some variants lack expert review
- Population frequency data incomplete for some variants
- Historical nomenclature variations may exist

## Data Quality Notes

BRCA Exchange prioritizes ENIGMA expert-reviewed classifications when available, as these represent the highest confidence assessments. Discordances between source databases are flagged for user review. The aggregation approach may show conflicting classifications, which should be interpreted in context of the underlying evidence and source credibility.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [ClinVar](../../1.1.variant.repositories/clinvar/README.md) - Source data
- [OncoKB](../oncokb/README.md) - Actionability
- [CIViC](../civic/README.md) - Clinical evidence
