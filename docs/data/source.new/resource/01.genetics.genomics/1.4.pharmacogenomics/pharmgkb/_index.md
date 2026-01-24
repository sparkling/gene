---
id: pharmgkb
title: "PharmGKB"
type: data-source
category: genetics
subcategory: pharmacogenomics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [pharmacogenomics, drug-response, clinical, guidelines, pgx]
---

# PharmGKB

**Category:** [Genetics & Genomics](../../_index.md) > [Pharmacogenomics](../_index.md)

## Overview

PharmGKB (Pharmacogenomics Knowledge Base) is the leading resource for pharmacogenomics information, curating knowledge about the impact of genetic variation on drug response. It aggregates and annotates pharmacogenomic data from the published literature, clinical guidelines, and drug labels.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Drug | PharmGKB Accession | PA449053 |
| Gene | PharmGKB Accession | PA124 |
| Pathway | PharmGKB Accession | PA165959313 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://www.pharmgkb.org/ | Interactive search |
| API | REST endpoints | Programmatic access |
| Downloads | Bulk files | TSV format |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY-SA 4.0 |
| Commercial Use | Yes (with attribution) |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [CPIC](../cpic/_index.md) - Clinical guidelines
- [PharmVar](../pharmvar/_index.md) - Haplotype definitions
