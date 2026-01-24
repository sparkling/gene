---
id: pharmvar
title: "PharmVar"
type: data-source
category: genetics
subcategory: pharmacogenomics
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [pharmacogenomics, haplotype, star-allele, nomenclature, cyp]
---

# PharmVar

**Category:** [Genetics & Genomics](../../_index.md) > [Pharmacogenomics](../_index.md)

## Overview

PharmVar (Pharmacogene Variation Consortium) is the central repository for pharmacogene haplotype (star allele) definitions and nomenclature. It standardizes the naming and definition of pharmacogene alleles, ensuring consistent interpretation of pharmacogenomic test results worldwide.

PharmVar maintains the official nomenclature for key pharmacogenes including CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP3A5, and other genes involved in drug metabolism and response. Each allele definition includes the defining variants, their positions on reference sequences, and functional status when known.

The database evolved from the CYP Allele Nomenclature Database and now provides a comprehensive resource for all pharmacogenes. PharmVar allele definitions are used by CPIC, DPWG, and clinical laboratories for consistent genotype reporting.

## Key Statistics

| Metric | Value |
|--------|-------|
| Genes | 30+ |
| Total Alleles | 3,000+ |
| CYP2D6 Alleles | 150+ |
| CYP2C19 Alleles | 40+ |
| Reference Sequences | Gene-specific |

## Primary Use Cases

1. Star allele definition lookup
2. Genotyping assay design
3. Haplotype-to-phenotype translation
4. Variant-to-allele assignment

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Allele | Gene*number | CYP2D6*4 |
| Core Allele | Base definition | CYP2D6*4.001 |
| Sub-allele | With variants | CYP2D6*4.013 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://www.pharmvar.org/ | Interactive search |
| API | REST endpoints | Programmatic access |
| Downloads | Gene-specific files | TSV, FASTA |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access |
| Commercial Use | Yes |
| Citation | Required |

## See Also

- [PharmGKB](../pharmgkb/_index.md) - Clinical annotations
- [CPIC](../cpic/_index.md) - Prescribing guidelines
- [DPWG](../dpwg/_index.md) - Dutch guidelines
