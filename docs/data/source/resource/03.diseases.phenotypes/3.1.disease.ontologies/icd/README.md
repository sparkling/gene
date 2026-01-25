---
id: icd
title: "International Classification of Diseases (ICD)"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: disease.ontologies
tags:
  - ontology
  - classification
  - diseases
  - who
  - clinical-coding
  - icd-10
  - icd-11
---

# International Classification of Diseases (ICD)

## Overview

The International Classification of Diseases (ICD) is the global standard for diagnostic health information maintained by the World Health Organization (WHO). ICD provides a common language for recording, reporting, and monitoring diseases, enabling the comparison of health data across populations, time periods, and countries.

ICD-10, the tenth revision, has been in widespread clinical use since 1994 and contains approximately 70,000 codes for diseases, signs, symptoms, and external causes of injury or disease. ICD-11, released in 2019 and effective from January 2022, represents a major modernization with improved digital functionality, more detailed coding, and better integration with electronic health records.

The classification system is hierarchical, with chapters covering broad disease categories (e.g., infectious diseases, neoplasms, circulatory system disorders) that subdivide into increasingly specific codes. ICD codes are essential for healthcare administration, epidemiological research, clinical documentation, and insurance reimbursement worldwide.

## Key Statistics

| Metric | Value |
|--------|-------|
| ICD-10 Codes | ~70,000 |
| ICD-11 Entities | ~80,000 |
| Chapters (ICD-10) | 22 |
| Chapters (ICD-11) | 26 |
| Countries Using | 194 |

## Primary Use Cases

1. Clinical diagnosis coding in healthcare settings
2. Health statistics reporting and epidemiological surveillance
3. Insurance claims processing and reimbursement
4. Cross-referencing clinical phenotypes to research databases
5. Mortality and morbidity statistics compilation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| ICD-10-CM | `[A-Z][0-9]{2}(\.[0-9]{1,4})?` | E11.9 (Type 2 diabetes) |
| ICD-10 | `[A-Z][0-9]{2}(\.[0-9])?` | I21.0 (MI) |
| ICD-11 | `[A-Z0-9]{4,6}` | 5A11 (Type 2 diabetes) |
| Block Code | `[A-Z][0-9]{2}-[A-Z][0-9]{2}` | E10-E14 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| ICD-11 Browser | https://icd.who.int/browse11/ | Official WHO browser |
| ICD-10 Browser | https://icd.who.int/browse10/ | ICD-10 version |
| ICD API | https://icd.who.int/icdapi | RESTful API access |
| ICD-10-CM | https://www.cdc.gov/nchs/icd/icd10cm.htm | US clinical modification |

## Data Formats

| Format | Source | Notes |
|--------|--------|-------|
| XML | WHO | Full classification export |
| JSON | ICD API | API responses |
| Tabular | CDC | ICD-10-CM text files |
| Foundation | WHO | ICD-11 foundation layer |

## Limitations

- Commercial use requires WHO license agreement
- API access requires registration and approval
- Clinical modifications (ICD-10-CM) vary by country
- Rare diseases less granularly coded than common conditions

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [MeSH](../mesh/README.md) - Medical vocabulary cross-reference
- [MONDO](../mondo/README.md) - Disease ontology with ICD mappings
