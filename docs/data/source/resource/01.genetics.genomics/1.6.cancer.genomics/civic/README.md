---
id: civic
title: "CIViC"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: cancer.genomics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - cancer
  - clinical
  - interpretation
  - community
  - evidence
---

# CIViC

CIViC (Clinical Interpretation of Variants in Cancer) is an open-source, community-driven knowledgebase of clinical interpretations for cancer variants. Developed by Washington University, it provides crowdsourced curation of variant-level evidence linking genomic alterations to clinical outcomes, therapies, and diagnostic/prognostic implications.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Washington University School of Medicine |
| **Website** | https://civicdb.org/ |
| **Update Frequency** | Continuous community curation |
| **Records** | 8,000+ evidence items |
| **Latest Release** | Current (continuous) |

Each evidence item in CIViC links a variant to a specific disease, evidence type (predictive, diagnostic, prognostic, predisposing, functional, oncogenic), evidence level (validated, clinical, case study, preclinical), and clinical significance. The wiki-style platform allows community contributions with expert review.

CIViC follows the AMP/ASCO/CAP guidelines for somatic variant interpretation and integrates with other knowledgebases through the VICC (Variant Interpretation for Cancer Consortium) harmonization efforts.

## Key Statistics

| Metric | Value |
|--------|-------|
| Genes | 500+ |
| Variants | 3,000+ |
| Evidence Items | 8,000+ |
| Assertions | 200+ |
| Contributors | 400+ |

## Primary Use Cases

1. Cancer variant interpretation
2. Clinical evidence lookup
3. Knowledge curation
4. Guideline-based classification

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Gene | `CIViC ID` | `5 (BRAF)` |
| Variant | `CIViC ID` | `12 (V600E)` |
| Evidence | `CIViC ID` | `EID1234` |

## Limitations

- Community curation may include errors
- Evidence quality varies by item
- Coverage incomplete for rare variants
- Some evidence items awaiting expert review
- Therapy recommendations may not reflect current practice

## Data Quality Notes

CIViC implements a trust-based curation model with evidence ratings and expert review. Assertions represent higher-confidence curated statements following AMP/ASCO/CAP guidelines. Evidence items should be evaluated based on their evidence level (A-E) and rating for clinical use. The CC0 license enables unrestricted use and redistribution.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [OncoKB](../oncokb/README.md) - Expert curation
- [COSMIC](../cosmic/README.md) - Somatic mutations
- [cBioPortal](../cbioportal/README.md) - Genomics data
