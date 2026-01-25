---
id: pharmgkb
title: "PharmGKB"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: pharmacogenomics
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - pharmacogenomics
  - drug-response
  - clinical
  - guidelines
  - pgx
---

# PharmGKB

PharmGKB (Pharmacogenomics Knowledge Base) is the leading resource for pharmacogenomics information, curating knowledge about the impact of genetic variation on drug response. It aggregates and annotates pharmacogenomic data from the published literature, clinical guidelines, and drug labels.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | Stanford University |
| **Website** | https://www.pharmgkb.org/ |
| **Update Frequency** | Continuous |
| **Records** | 20,000+ clinical annotations |
| **Latest Release** | Current (continuous updates) |

The database provides clinical annotations with evidence levels (1A-4) indicating the strength of gene-drug associations, variant annotations describing functional effects, drug pathway diagrams, and links to clinical guidelines from CPIC and DPWG. PharmGKB also maintains curated data on drug labels containing pharmacogenomic information from FDA, EMA, PMDA, and HCSC.

PharmGKB serves as the curation hub for CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines, which provide prescribing recommendations based on genetic test results.

## Key Statistics

| Metric | Value |
|--------|-------|
| Clinical Annotations | 20,000+ |
| Variant Annotations | 170,000+ |
| Drugs | 800+ |
| Genes | 150+ (VIP genes) |
| CPIC Guidelines | 33 |

## Primary Use Cases

1. Drug-gene interaction lookup
2. Clinical pharmacogenomics implementation
3. Dosing guideline reference
4. Pharmacogenomic research

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Drug | `PharmGKB Accession` | `PA449053` |
| Gene | `PharmGKB Accession` | `PA124` |
| Pathway | `PharmGKB Accession` | `PA165959313` |

## Limitations

- Evidence levels vary; not all annotations are high-confidence
- Some variant annotations are computationally predicted
- Drug labels may conflict with clinical guidelines
- Coverage biased toward well-studied genes (CYP450 family)
- API rate limits apply

## Data Quality Notes

PharmGKB uses a structured evidence hierarchy with Level 1A (CPIC/DPWG guideline) representing the highest confidence and Level 4 representing preliminary evidence. Users should prioritize Level 1-2 annotations for clinical applications. All curated content is reviewed by PharmGKB staff with regular updates as new evidence emerges.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [CPIC](../cpic/README.md) - Clinical guidelines
- [PharmVar](../pharmvar/README.md) - Haplotype definitions
