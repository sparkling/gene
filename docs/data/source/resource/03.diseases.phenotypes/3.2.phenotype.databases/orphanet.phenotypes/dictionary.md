# Orphanet Phenotypes - Data Dictionary

## Overview

This data dictionary documents the schema for Orphanet Phenotype Annotations.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | orphanet.phenotypes |
| **Name** | Orphanet Phenotypes |
| **Parent** | 3.2.phenotype.databases |
| **Total Fields** | 15+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### HPO Disease Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Orphanet disease ID | 558 |
| DiseaseName | string | 1:1 | Yes | Disease name | Marfan syndrome |
| HPOId | string | 1:1 | Yes | HPO term ID | HP:0001166 |
| HPOTerm | string | 1:1 | Yes | Phenotype name | Arachnodactyly |
| HPOFrequencyId | string | 1:1 | No | Frequency HPO term | HP:0040281 |
| HPOFrequencyName | string | 1:1 | No | Frequency label | Very frequent (99-80%) |
| Diagnostic | boolean | 1:1 | No | Diagnostic criterion | true |
| ValidationStatus | string | 1:1 | No | Curation status | Validated |

### Disorder Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Orphanet ID | 558 |
| Name | string | 1:1 | Yes | Disease name | Marfan syndrome |
| DisorderType | string | 1:1 | Yes | Entity type | Disease |
| DisorderGroup | string | 1:1 | Yes | Group type | Disorder |
| AverageAgeOfOnset | array | 1:N | No | Onset ages | [Infancy, Childhood] |
| TypeOfInheritance | array | 1:N | No | Inheritance modes | [AD] |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| ORPHAcode | [0-9]+ | 558 | Orphanet disease ID |
| Orphanet ID | Orphanet:[0-9]+ | Orphanet:558 | CURIE format |
| HPO ID | HP:[0-9]{7} | HP:0001166 | Phenotype term |
| HPO Frequency | HP:004028[0-5] | HP:0040281 | Frequency term |

---

## Enumerations

### Frequency Categories

| Category | HPO Term | Percentage | Description |
|----------|----------|------------|-------------|
| Obligate | HP:0040280 | 100% | Always present |
| Very frequent | HP:0040281 | 80-99% | Almost always |
| Frequent | HP:0040282 | 30-79% | Common |
| Occasional | HP:0040283 | 5-29% | Sometimes |
| Very rare | HP:0040284 | 1-4% | Rarely |
| Excluded | HP:0040285 | 0% | Never present |

### Diagnostic Criteria

| Value | Description |
|-------|-------------|
| true | Phenotype is diagnostic criterion |
| false | Not a diagnostic criterion |
| null | Not specified |

### Validation Status

| Status | Description |
|--------|-------------|
| Validated | Expert validated |
| Not yet validated | Pending validation |

### Disorder Types

| Type | Description |
|------|-------------|
| Disease | Clinical disease entity |
| Malformation syndrome | Congenital malformation |
| Morphological anomaly | Structural anomaly |
| Clinical syndrome | Clinical pattern |

### Age of Onset Categories

| Category | Description |
|----------|-------------|
| Antenatal | Before birth |
| Neonatal | First 28 days |
| Infancy | 1 month - 1 year |
| Childhood | 1-11 years |
| Adolescent | 12-18 years |
| Adult | 18+ years |
| Elderly | 65+ years |
| All ages | Any age |

### Inheritance Types

| Code | Name |
|------|------|
| AD | Autosomal dominant |
| AR | Autosomal recessive |
| XLD | X-linked dominant |
| XLR | X-linked recessive |
| MT | Mitochondrial |
| MF | Multifactorial |

---

## Entity Relationships

### Disease to Phenotypes
- **Cardinality:** 1:N
- **Description:** Each disease has multiple phenotype annotations
- **Key Fields:** ORPHAcode, HPOId

### Phenotype to Frequency
- **Cardinality:** 1:1
- **Description:** Each annotation has one frequency
- **Key Fields:** HPOId, HPOFrequencyId

### Disease to Onset Ages
- **Cardinality:** 1:N
- **Description:** Diseases can have multiple onset age ranges
- **Key Fields:** ORPHAcode, AverageAgeOfOnset

### Disease to Inheritance
- **Cardinality:** 1:N
- **Description:** Diseases can have multiple inheritance patterns
- **Key Fields:** ORPHAcode, TypeOfInheritance

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HPO | Human Phenotype Ontology | Phenotype vocabulary |
| ORPHA | Orphanet | Disease identifier prefix |
| AD | Autosomal Dominant | Inheritance pattern |
| AR | Autosomal Recessive | Inheritance pattern |
| XLR | X-Linked Recessive | Inheritance pattern |
| XLD | X-Linked Dominant | Inheritance pattern |
| MT | Mitochondrial | Inheritance pattern |
| MF | Multifactorial | Inheritance pattern |
| HOOM | HPO-Orphanet Ontological Module | Integration resource |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| HPO | HP ID | Phenotype terms |
| OMIM | MIM number | Disease cross-reference |
| MONDO | MONDO ID | Disease ontology |
| UMLS | CUI | Concept mapping |

---

## Data Quality Notes

1. **Annotated Diseases:** 4,337+ rare diseases
2. **HPO Annotations:** 100,000+ phenotype links
3. **Frequency Data:** 6 standardized frequency categories
4. **Expert Curation:** 41-country expert network
5. **Diagnostic Criteria:** Diagnostic phenotypes flagged
6. **Monthly Updates:** Regular data releases
7. **XML Download:** Product 4 from Orphadata
8. **CC BY 4.0:** Free for commercial use

