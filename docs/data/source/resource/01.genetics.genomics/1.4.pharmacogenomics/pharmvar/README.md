---
id: pharmvar
title: "PharmVar"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: pharmacogenomics
tier: 2
status: active
last_updated: 2026-01-25
tags:
  - pharmacogenomics
  - haplotype
  - star-allele
  - nomenclature
  - cyp
---

# PharmVar

PharmVar (Pharmacogene Variation Consortium) is the central repository for pharmacogene haplotype (star allele) definitions and nomenclature. It standardizes the naming and definition of pharmacogene alleles, ensuring consistent interpretation of pharmacogenomic test results worldwide.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | PharmVar Consortium |
| **Website** | https://www.pharmvar.org/ |
| **Update Frequency** | Continuous |
| **Records** | 3,000+ alleles |
| **Latest Release** | Current (continuous) |

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

| Identifier | Format | Example |
|------------|--------|---------|
| Allele | `Gene*number` | `CYP2D6*4` |
| Core Allele | `Base definition` | `CYP2D6*4.001` |
| Sub-allele | `With variants` | `CYP2D6*4.013` |

## Limitations

- Allele function not always established
- Some alleles population-specific (rare in other groups)
- Complex rearrangements difficult to represent
- Gene conversion events may complicate calling
- Historical nomenclature inconsistencies exist

## Data Quality Notes

PharmVar allele definitions are curated by expert working groups with variant validation and reference sequence alignment. Function assignments are based on published literature and consensus expert opinion. The sub-allele system captures additional variants on core haplotype backgrounds, enabling precise genotype reporting.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [PharmGKB](../pharmgkb/README.md) - Clinical annotations
- [CPIC](../cpic/README.md) - Prescribing guidelines
- [DPWG](../dpwg/README.md) - Dutch guidelines
