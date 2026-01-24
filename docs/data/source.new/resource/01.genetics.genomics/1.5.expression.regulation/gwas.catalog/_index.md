---
id: gwas-catalog
title: "GWAS Catalog"
type: data-source
category: genetics
subcategory: expression.regulation
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [gwas, association, trait, snp, phenotype]
---

# GWAS Catalog

**Category:** [Genetics & Genomics](../../_index.md) > [Expression & Regulation](../_index.md)

## Overview

The NHGRI-EBI GWAS Catalog is the comprehensive curated collection of published genome-wide association studies and their associated SNP-trait associations. It provides standardized, quality-controlled GWAS results with consistent trait mapping using the Experimental Factor Ontology (EFO).

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Study | GCST + 6 digits | GCST000854 |
| SNP | RS ID | rs7329174 |
| EFO Trait | EFO_XXXXXXX | EFO_0001060 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://www.ebi.ac.uk/gwas/ | Interactive search |
| REST API | https://www.ebi.ac.uk/gwas/rest/api | JSON (HAL format) |
| FTP | http://ftp.ebi.ac.uk/pub/databases/gwas/ | Bulk downloads |
| Summary Stats | API and FTP | Full GWAS results |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Citation | Required |

## See Also

- [Schema Documentation](./schema.md)
- [GTEx](../gtex/_index.md) - eQTL colocalization
- [UK Biobank](../../1.3.population.genetics/uk.biobank/_index.md) - Source studies
