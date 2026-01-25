---
id: orphanet.phenotypes
title: "Orphanet Phenotype Annotations"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: phenotype.databases
tags:
  - phenotypes
  - rare-diseases
  - hpo
  - clinical-features
  - orphanet
---

# Orphanet Phenotype Annotations

## Overview

Orphanet Phenotype Annotations (product4) provides comprehensive mappings between rare diseases and their associated clinical features using Human Phenotype Ontology (HPO) terms. This dataset enables systematic phenotyping of rare diseases and supports computational approaches to differential diagnosis and patient matching.

Each disease-phenotype association in the Orphanet dataset includes frequency information (obligate, very frequent, frequent, occasional, very rare) and validation status. The annotations are curated by clinical experts in rare diseases across the Orphanet network of 41 countries, ensuring high quality and clinical relevance.

The phenotype data integrates with Orphanet's broader rare disease knowledge base including epidemiology, inheritance patterns, age of onset, and expert center information. This integration makes Orphanet Phenotypes particularly valuable for rare disease diagnosis algorithms that combine phenotypic and molecular evidence.

## Key Statistics

| Metric | Value |
|--------|-------|
| Annotated Diseases | 4,337+ |
| HPO Annotations | 100,000+ |
| Frequency Categories | 6 |
| Languages | 8 |
| Update Frequency | Monthly |

## Primary Use Cases

1. Rare disease differential diagnosis based on phenotype matching
2. Clinical decision support for genetic testing
3. Research cohort identification by phenotype
4. Phenotype frequency analysis across rare diseases
5. Training phenotype prediction models

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Orphanet ID | `Orphanet:[0-9]+` | Orphanet:558 |
| ORPHA Code | `[0-9]+` | 558 (Marfan) |
| HPO Term | `HP:[0-9]{7}` | HP:0001166 |
| Frequency | HPO frequency terms | HP:0040281 (Very frequent) |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Orphadata | https://www.orphadata.com/ | Science data portal |
| Product 4 | https://www.orphadata.com/en_product4.xml | Phenotype file |
| REST API | https://api.orphanet.org/ | Programmatic access |
| HOOM | https://www.orphadata.com/hoom/ | HPO-Orphanet mappings |

## Data Formats

| Format | File | Size |
|--------|------|------|
| XML | en_product4.xml | ~15 MB |
| JSON | API responses | Variable |
| RDF | Available via HOOM | Linked data |

## Frequency Categories

| Category | HPO Term | Percentage |
|----------|----------|------------|
| Obligate | HP:0040280 | 100% |
| Very frequent | HP:0040281 | 80-99% |
| Frequent | HP:0040282 | 30-79% |
| Occasional | HP:0040283 | 5-29% |
| Very rare | HP:0040284 | 1-4% |
| Excluded | HP:0040285 | 0% |

## Limitations

- Rare diseases only; common conditions not covered
- Frequency estimates based on literature reports
- Not all diseases have complete phenotype annotation
- Expert center coverage varies by country

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [HPO](../hpo/README.md) - Human Phenotype Ontology
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/README.md) - Main Orphanet database
