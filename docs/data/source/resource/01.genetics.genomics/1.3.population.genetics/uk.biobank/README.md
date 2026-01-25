---
id: uk-biobank
title: "UK Biobank"
type: source
parent: ../README.md
category: genetics.genomics
subcategory: population.genetics
tier: 2
status: active
last_updated: 2026-01-25
tags:
  - population
  - biobank
  - phenotype
  - gwas
  - longitudinal
---

# UK Biobank

UK Biobank is a large-scale biomedical database containing genetic, physical, and health information from 500,000 volunteer participants aged 40-69 from across the United Kingdom. It is one of the most comprehensive prospective cohort studies ever conducted, linking genetic data with extensive phenotypic measurements and health outcomes.

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | UK Biobank Ltd |
| **Website** | https://www.ukbiobank.ac.uk/ |
| **Update Frequency** | Continuous data additions |
| **Records** | 500,000 participants |
| **Latest Release** | Ongoing releases |

The genetic data includes genotyping array data (all participants), exome sequencing (all participants), and whole genome sequencing (expanding coverage). Phenotypic data encompasses physical measurements, lifestyle factors, biomarker assays, imaging data, and longitudinal health records through linkage with NHS and death registries.

UK Biobank has enabled thousands of genetic association studies and serves as a critical resource for understanding the genetic basis of common diseases, pharmacogenomics, and population health research.

## Key Statistics

| Metric | Value |
|--------|-------|
| Participants | 500,000 |
| Genotyped | 500,000 |
| Exome Sequenced | 470,000+ |
| WGS | 200,000+ |
| Phenotypes | 10,000+ |

## Primary Use Cases

1. Genome-wide association studies
2. Phenome-wide association studies
3. Mendelian randomization
4. Polygenic risk score development

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| Participant | `eid (encrypted)` | `1234567` |
| Field | `field-instance-array` | `21001-0-0` |
| Data Category | `Category ID` | `100001` |

## Limitations

- Access requires approved research application (weeks to months)
- Participants predominantly white British (limited diversity)
- Age range 40-69 at recruitment (no pediatric data)
- Healthy volunteer bias (healthier than general population)
- Commercial access requires separate licensing

## Data Quality Notes

UK Biobank implements extensive quality control including genotype QC, sample QC, and phenotype validation. Genetic data is imputed to the Haplotype Reference Consortium panel. The Research Analysis Platform (RAP) provides cloud-based analysis without requiring data download, addressing storage and security concerns.

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
- [gnomAD](../gnomad/README.md) - Population frequencies
- [GWAS Catalog](../../1.5.expression.regulation/gwas.catalog/README.md) - GWAS results
- [TOPMed](../topmed/README.md) - US equivalent
