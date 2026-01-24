---
id: decipher
title: "DECIPHER - Database of Chromosomal Imbalance and Phenotype in Humans using Ensembl Resources"
type: data-source
category: diseases
subcategory: rare.orphan.diseases
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [rare-diseases, cnv, developmental-disorders, clinical-genetics, genomic-variants]
---

# DECIPHER

**Category:** [Diseases & Phenotypes](../../_index.md) > [Rare & Orphan Diseases](../_index.md)

## Overview

DECIPHER (DatabasE of genomiC varIation and Phenotype in Humans using Ensembl Resources) is an interactive web-based database that incorporates clinical interpretation of rare copy number variants (CNVs) and sequence variants. The database facilitates the identification of genes critical to human development and health by enabling clinical data sharing between rare disease specialists worldwide.

Hosted by the Wellcome Sanger Institute, DECIPHER collects and shares anonymized patient phenotypes and genomic variants from clinical genetics centers across 35 countries. The platform uses HPO terminology for phenotype annotation and provides tools for variant interpretation including gene dosage sensitivity scores, constraint metrics, and overlap with known pathogenic variants.

DECIPHER is particularly valuable for interpreting novel CNVs in developmental disorders, allowing clinicians to compare patients with similar genetic lesions and identify candidate disease genes within large chromosomal regions. The database also integrates with DDD (Deciphering Developmental Disorders) study data, providing research-grade variant interpretations.

## Key Statistics

| Metric | Value |
|--------|-------|
| Patient Records | 50,000+ |
| Contributing Centers | 300+ |
| Countries | 35+ |
| CNV Records | 40,000+ |
| SNV Records | 60,000+ |

## Primary Use Cases

1. CNV interpretation in clinical genetics
2. Novel gene-disease relationship discovery
3. Patient phenotype matching for diagnosis
4. Dosage sensitivity assessment
5. Developmental disorder gene identification

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| DECIPHER ID | `[0-9]+` | 250000 |
| Patient ID | `DECIPHER:[0-9]+` | DECIPHER:250001 |
| HPO Term | `HP:[0-9]{7}` | HP:0001249 |
| Ensembl Gene | `ENSG[0-9]{11}` | ENSG00000146648 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| DECIPHER Portal | https://www.deciphergenomics.org/ | Web interface |
| Gene Browser | https://www.deciphergenomics.org/browser | Genome view |
| API | https://www.deciphergenomics.org/api | Data access |
| Variant Search | https://www.deciphergenomics.org/search | Query interface |

## Data Types

| Data Type | Description | Access |
|-----------|-------------|--------|
| CNV | Copy number variants with coordinates | Open summary |
| SNV | Sequence variants | Open summary |
| Phenotypes | HPO-coded clinical features | Open summary |
| Inheritance | Mode of inheritance | Open |
| Patient Details | Individual-level data | Controlled |

## Scoring Metrics

| Metric | Description |
|--------|-------------|
| HI Score | Haploinsufficiency prediction |
| pLI | Probability of LoF intolerance |
| LOEUF | LoF observed/expected upper bound |
| Pathogenicity | Variant classification |

## License

| Aspect | Value |
|--------|-------|
| License | DECIPHER Data Sharing Agreement |
| Public Data | Aggregate statistics, gene-level |
| Controlled Data | Patient-level (requires agreement) |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [Orphanet](../orphanet/_index.md) - Rare disease portal
- [HPO](../../3.2.phenotype.databases/hpo/_index.md) - Phenotype terminology
