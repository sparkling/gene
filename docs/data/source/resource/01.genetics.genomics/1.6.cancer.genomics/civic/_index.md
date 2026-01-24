---
id: civic
title: "CIViC"
type: data-source
category: genetics
subcategory: cancer.genomics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [cancer, clinical, interpretation, community, evidence]
---

# CIViC

**Category:** [Genetics & Genomics](../../_index.md) > [Cancer Genomics](../_index.md)

## Overview

CIViC (Clinical Interpretation of Variants in Cancer) is an open-source, community-driven knowledgebase of clinical interpretations for cancer variants. Developed by Washington University, it provides crowdsourced curation of variant-level evidence linking genomic alterations to clinical outcomes, therapies, and diagnostic/prognostic implications.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Gene | CIViC ID | 5 (BRAF) |
| Variant | CIViC ID | 12 (V600E) |
| Evidence | CIViC ID | EID1234 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://civicdb.org/ | Interactive curation |
| API | GraphQL | Programmatic access |
| Downloads | Nightly releases | TSV, JSON |
| GitHub | civic-server | Open source |

## License

| Aspect | Value |
|--------|-------|
| License | CC0 1.0 (Public Domain) |
| Commercial Use | Yes |
| Contribution | Open to all |

## See Also

- [OncoKB](../oncokb/_index.md) - Expert curation
- [COSMIC](../cosmic/_index.md) - Somatic mutations
- [cBioPortal](../cbioportal/_index.md) - Genomics data
