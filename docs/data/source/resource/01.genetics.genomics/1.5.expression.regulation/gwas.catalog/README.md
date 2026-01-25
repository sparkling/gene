---
id: gwas-catalog
title: "GWAS Catalog"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: expression.regulation
tier: 1
status: active
last_updated: 2026-01-25
tags:
  - gwas
  - association
  - trait
  - snp
  - phenotype
---

# GWAS Catalog

The NHGRI-EBI GWAS Catalog is the comprehensive curated collection of published genome-wide association studies and their associated SNP-trait associations. It provides standardized, quality-controlled GWAS results with consistent trait mapping using the Experimental Factor Ontology (EFO).

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | NHGRI-EBI |
| **Website** | https://www.ebi.ac.uk/gwas/ |
| **Update Frequency** | Weekly |
| **Records** | 1,058,471+ associations |
| **Latest Release** | Current (weekly updates) |

The catalog includes study-level information (sample sizes, ancestry, publication details), association-level data (p-values, effect sizes, risk alleles), and mapped genes and traits. Summary statistics from many studies are also available for download and imputation.

The GWAS Catalog is essential for variant interpretation, identifying known disease associations, and meta-analysis. EFO trait mapping enables cross-study comparisons and disease hierarchy navigation.

## Key Statistics

| Metric | Value |
|--------|-------|
| Studies | 186,237 |
| Associations | 1,058,471 |
| SNPs | 512,069 |
| Summary Statistics | 155,485 |
| EFO Traits | 21,004 |

## Primary Use Cases

1. Variant-trait association lookup
2. Disease gene discovery
3. Polygenic risk score development
4. GWAS meta-analysis

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Study | `GCST + 6 digits` | `GCST000854` |
| SNP | `RS ID` | `rs7329174` |
| EFO Trait | `EFO_XXXXXXX` | `EFO_0001060` |

## Limitations

- Publication bias toward significant results
- Ancestry bias (historically European-dominant)
- Variable effect size reporting across studies
- Some associations not replicated
- Summary statistics not available for all studies

## Data Quality Notes

The GWAS Catalog applies standardized curation including p-value threshold enforcement (p < 1e-5), variant mapping to current genome builds, and EFO trait harmonization. Quality flags indicate potential issues such as small sample size or limited ancestry diversity. Summary statistics undergo additional QC before release.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [GTEx](../gtex/README.md) - eQTL colocalization
- [UK Biobank](../../1.3.population.genetics/uk.biobank/README.md) - Source studies
