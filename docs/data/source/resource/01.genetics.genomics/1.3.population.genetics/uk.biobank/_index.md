---
id: uk-biobank
title: "UK Biobank"
type: data-source
category: genetics
subcategory: population.genetics
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [population, biobank, phenotype, gwas, longitudinal]
---

# UK Biobank

**Category:** [Genetics & Genomics](../../_index.md) > [Population Genetics](../_index.md)

## Overview

UK Biobank is a large-scale biomedical database containing genetic, physical, and health information from 500,000 volunteer participants aged 40-69 from across the United Kingdom. It is one of the most comprehensive prospective cohort studies ever conducted, linking genetic data with extensive phenotypic measurements and health outcomes.

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

| Identifier | Pattern | Example |
|------------|---------|---------|
| Participant | eid (encrypted) | 1234567 |
| Field | field-instance-array | 21001-0-0 |
| Data Category | Category ID | 100001 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Portal | https://www.ukbiobank.ac.uk/ | Application required |
| RAP | Research Analysis Platform | Cloud analysis |
| DNA Nexus | UKB-RAP | Approved researchers |
| Bulk Data | Downloads | Large file access |

## License

| Aspect | Value |
|--------|-------|
| License | UK Biobank Access Agreement |
| Commercial Use | Requires approval |
| Application | Required for all access |

## See Also

- [gnomAD](../gnomad/_index.md) - Population frequencies
- [GWAS Catalog](../../1.5.expression.regulation/gwas.catalog/_index.md) - GWAS results
- [TOPMed](../topmed/_index.md) - US equivalent
