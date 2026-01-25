---
id: orphanet
title: "Orphanet - The Portal for Rare Diseases and Orphan Drugs"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: rare.orphan.diseases
tags:
  - rare-diseases
  - orphan-drugs
  - clinical-resources
  - gene-disease
  - ordo
---

# Orphanet - The Portal for Rare Diseases and Orphan Drugs

## Overview

Orphanet is the reference portal for information on rare diseases and orphan drugs, serving patients, healthcare professionals, and researchers worldwide. Established in 1997 by the French National Institute for Health and Medical Research (INSERM), Orphanet now operates as a consortium of 41 countries coordinated by the INSERM.

The database provides expert-validated information on over 6,500 rare diseases, including clinical descriptions, epidemiology, gene-disease associations, diagnostic tests, and expert centers. Orphanet uses a unique nomenclature (ORPHA codes) that is cross-referenced to other terminologies including ICD-10, OMIM, UMLS, and MeSH.

Orphanet data is made available through Orphadata for research use, including the ORDO ontology (Orphanet Rare Disease Ontology) which provides an OWL representation of the disease classification. The platform also maintains a directory of orphan drugs at various stages of development, supporting drug development and access for rare disease patients.

## Key Statistics

| Metric | Value |
|--------|-------|
| Rare Diseases | 6,528 |
| Linked Genes | 4,512 |
| Diagnostic Tests | 36,595 |
| Expert Centres | 8,722 |
| Patient Organizations | 3,500+ |

## Primary Use Cases

1. Rare disease diagnosis and clinical guidance
2. Gene-disease association research
3. Expert center and patient organization lookup
4. Orphan drug development tracking
5. Disease nomenclature standardization

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| ORPHA Code | `Orphanet:[0-9]+` | Orphanet:558 |
| ORDO ID | `ORDO:[0-9]+` | ORDO:558 |
| OMIM XREF | `OMIM:[0-9]{6}` | OMIM:154700 |
| ICD-10 XREF | `ICD-10:[A-Z][0-9]+` | ICD-10:Q87.4 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Orphanet Portal | https://www.orpha.net/ | Web interface |
| Orphadata | https://www.orphadata.com/ | Science downloads |
| ORDO Ontology | https://www.ebi.ac.uk/ols4/ontologies/ordo | OWL format |
| REST API | https://api.orphanet.org/ | Programmatic access |

## Data Products

| Product | File | Description |
|---------|------|-------------|
| Disease List | en_product1.xml | All diseases |
| Gene-Disease | en_product6.xml | Gene associations |
| Phenotypes | en_product4.xml | HPO annotations |
| Epidemiology | en_product9.xml | Prevalence data |
| Classifications | en_product3.xml | Disease hierarchy |

## Limitations

- Rare diseases only (prevalence < 1 in 2000)
- Expert center information varies by country
- English translations may lag native language updates
- Orphan drug data focused on European market

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [OMIM](../../3.2.phenotype.databases/omim/README.md) - Mendelian inheritance database
- [HPO](../../3.2.phenotype.databases/hpo/README.md) - Phenotype ontology
