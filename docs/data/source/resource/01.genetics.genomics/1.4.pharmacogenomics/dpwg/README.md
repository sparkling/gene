---
id: dpwg
title: "DPWG"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: pharmacogenomics
tier: 2
status: active
last_updated: 2026-01-25
tags:
  - pharmacogenomics
  - guidelines
  - dutch
  - prescribing
  - european
---

# DPWG

DPWG (Dutch Pharmacogenetics Working Group) develops pharmacogenetics-based therapeutic recommendations for drug-gene interactions. Operating since 2005, DPWG guidelines are integrated into the Dutch national drug database (G-Standaard) and are used throughout the Netherlands healthcare system.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | KNMP (Royal Dutch Pharmacists Association) |
| **Website** | https://www.pharmgkb.org/page/dpwg |
| **Update Frequency** | Periodic updates |
| **Records** | 100+ gene-drug pairs |
| **Latest Release** | Current (continuous) |

DPWG guidelines follow a systematic methodology assessing clinical relevance of gene-drug interactions based on published evidence. Recommendations are categorized by urgency (requiring immediate action vs. optional adjustments) and provide specific dosing guidance or alternative drug recommendations.

The guidelines are available through PharmGKB and complement CPIC guidelines, with efforts ongoing to harmonize recommendations between the two organizations. DPWG guidelines are particularly valuable for European clinical implementation.

## Key Statistics

| Metric | Value |
|--------|-------|
| Gene-Drug Pairs | 100+ |
| Guidelines Active | 90+ |
| Genes Covered | 15+ |
| Evidence Levels | Strong, Moderate, Limited |
| Implementation | G-Standaard |

## Primary Use Cases

1. European clinical prescribing
2. Dutch EHR integration
3. Pharmacogenomics implementation
4. Cross-reference with CPIC

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Gene | `HGNC symbol` | `CYP2C19` |
| Drug | `Generic name` | `clopidogrel` |
| ATC Code | `WHO ATC` | `B01AC04` |

## Limitations

- Primarily European population data
- Dutch language for some primary documentation
- May differ from CPIC recommendations
- Limited to actionable gene-drug pairs
- Updates may lag behind new evidence

## Data Quality Notes

DPWG guidelines are developed through systematic literature review with structured evidence assessment. Recommendations are graded by clinical relevance and evidence strength. Cross-referencing with CPIC guidelines is recommended when discrepancies exist, as methodological differences may lead to different recommendations for the same gene-drug pair.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [PharmGKB](../pharmgkb/README.md) - English access
- [CPIC](../cpic/README.md) - US guidelines
- [PharmVar](../pharmvar/README.md) - Allele definitions
