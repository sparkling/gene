---
id: cpic
title: "CPIC"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: pharmacogenomics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - pharmacogenomics
  - guidelines
  - clinical
  - prescribing
  - dosing
---

# CPIC

CPIC (Clinical Pharmacogenetics Implementation Consortium) is an international consortium that creates freely available, peer-reviewed, evidence-based clinical practice guidelines for gene-drug pairs. CPIC guidelines help clinicians translate pharmacogenomic test results into actionable prescribing decisions.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | CPIC (hosted at cpicpgx.org) |
| **Website** | https://cpicpgx.org/ |
| **Update Frequency** | Continuous guideline updates |
| **Records** | 75+ gene-drug pairs |
| **Latest Release** | Current (continuous) |

Each CPIC guideline provides gene-specific information, drug-specific dosing recommendations, and clinical decision support resources including allele functionality tables, phenotype-to-diplotype translation, and therapeutic recommendations. Guidelines are regularly updated as new evidence emerges.

CPIC guidelines are developed through collaboration between PharmGKB (curation), CPIC members (guideline development), and clinical experts. The guidelines are designed to be implemented directly into electronic health record systems and clinical decision support tools.

## Key Statistics

| Metric | Value |
|--------|-------|
| Published Guidelines | 33 |
| Gene-Drug Pairs | 75+ |
| Priority Genes | 25 |
| Evidence Levels | A, B, C, D |
| Update Cycle | Continuous |

## Primary Use Cases

1. Clinical prescribing decisions
2. EHR clinical decision support
3. Laboratory test interpretation
4. Pharmacy dosing protocols

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Guideline | `DOI` | `10.1002/cpt.1234` |
| Gene | `HGNC symbol` | `CYP2D6` |
| Drug | `Generic name` | `codeine` |

## Limitations

- Limited to genes with actionable evidence
- May not cover all drugs in a class
- Recommendations primarily based on European ancestry data
- Some guidelines lack pediatric dosing
- Implementation requires EHR customization

## Data Quality Notes

CPIC guidelines undergo rigorous peer review and are published in Clinical Pharmacology & Therapeutics. Evidence is graded from A (strong) to D (weak), with corresponding recommendation strength. Guidelines with Level A evidence and strong recommendations represent the highest confidence for clinical implementation.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [PharmGKB](../pharmgkb/README.md) - Data curation
- [PharmVar](../pharmvar/README.md) - Allele definitions
- [DPWG](../dpwg/README.md) - Dutch guidelines
