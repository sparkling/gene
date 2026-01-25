# Orphanet Rare Diseases - Data Dictionary

## Overview

This data dictionary documents the schema for Orphanet rare disease database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | orphanet |
| **Name** | Orphanet |
| **Parent** | 3.5.rare.orphan.diseases |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Rare Disease

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Orphanet identifier | 558 |
| Name | string | 1:1 | Yes | Disease name | Marfan syndrome |
| ORPHAgroup | enum | 1:1 | Yes | Entity type | Disorder |
| DisorderType | object | 1:1 | Yes | Disease classification | Disease |
| Definition | string | 1:1 | No | Disease definition | A hereditary disorder... |
| SynonymList | array | 1:N | No | Alternative names | [MFS] |
| AverageAgeOfOnset | array | 1:N | No | Onset categories | [Childhood] |
| AverageAgeOfDeath | array | 1:N | No | Death categories | [Adult] |
| TypeOfInheritance | array | 1:N | No | Inheritance modes | [Autosomal dominant] |

### Gene Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Disease ID | 558 |
| GeneSymbol | string | 1:1 | Yes | HGNC symbol | FBN1 |
| GeneName | string | 1:1 | Yes | Full gene name | fibrillin 1 |
| HGNC | integer | 1:1 | No | HGNC ID | 3603 |
| GeneLocus | string | 1:1 | No | Chromosomal location | 15q21.1 |
| AssociationType | string | 1:1 | Yes | Relationship type | Disease-causing germline mutation |
| AssociationStatus | string | 1:1 | Yes | Curation status | Assessed |

### Epidemiological Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Disease ID | 558 |
| PrevalenceType | string | 1:1 | No | Prevalence type | Point prevalence |
| PrevalenceClass | string | 1:1 | No | Prevalence category | 1-5 / 10 000 |
| ValMoy | float | 1:1 | No | Mean value | 1.5 |
| PrevalenceGeographic | string | 1:1 | No | Geographic region | Europe |
| PrevalenceValidationStatus | string | 1:1 | No | Validation status | Validated |

### External References

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Disease ID | 558 |
| Source | string | 1:1 | Yes | External database | OMIM |
| Reference | string | 1:1 | Yes | External ID | 154700 |
| DisorderMappingRelation | string | 1:1 | No | Mapping type | E (exact) |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| ORPHAcode | [0-9]+ | 558 | Numeric identifier |
| Orphanet ID | Orphanet:[0-9]+ | Orphanet:558 | CURIE format |
| ORDO URI | http://www.orpha.net/ORDO/Orphanet_[0-9]+ | Orphanet_558 | Ontology URI |
| OMIM | [0-9]{6} | 154700 | OMIM cross-reference |
| ICD-10 | [A-Z][0-9]{2}\.[0-9] | Q87.4 | ICD cross-reference |
| HGNC | [0-9]+ | 3603 | Gene ID |

---

## Enumerations

### ORPHAgroup Types

| Type | Description |
|------|-------------|
| Disorder | Single clinical entity |
| Group of disorders | Related disease collection |
| Subtype of a disorder | Variant of parent disease |

### Disorder Types

| Type | Description |
|------|-------------|
| Disease | Clinical disease |
| Malformation syndrome | Congenital malformation |
| Morphological anomaly | Structural anomaly |
| Clinical syndrome | Clinical pattern |
| Biological anomaly | Biological abnormality |

### Gene-Disease Association Types

| Type | Description |
|------|-------------|
| Disease-causing germline mutation(s) in | Direct causative |
| Disease-causing germline mutation(s) (loss of function) | LOF mechanism |
| Disease-causing germline mutation(s) (gain of function) | GOF mechanism |
| Major susceptibility factor in | Risk factor |
| Modifying germline mutation in | Disease modifier |
| Role in phenotype of | Phenotype contribution |
| Biomarker tested in | Diagnostic marker |
| Candidate gene tested in | Research candidate |

### Inheritance Types

| Type | Description |
|------|-------------|
| Autosomal dominant | AD inheritance |
| Autosomal recessive | AR inheritance |
| X-linked dominant | XLD inheritance |
| X-linked recessive | XLR inheritance |
| Mitochondrial inheritance | MT inheritance |
| Multigenic/multifactorial | Complex inheritance |
| Y-linked | YL inheritance |
| Unknown | Not determined |

### Prevalence Classes

| Class | Range |
|-------|-------|
| >1 / 1000 | Common |
| 1-5 / 10 000 | Relatively common |
| 1-9 / 100 000 | Rare |
| 1-9 / 1 000 000 | Very rare |
| <1 / 1 000 000 | Ultra-rare |
| Unknown | Not available |

### Age of Onset Categories

| Category | Description |
|----------|-------------|
| Antenatal | Before birth |
| Neonatal | First 28 days |
| Infancy | 1 month - 2 years |
| Childhood | 2-11 years |
| Adolescent | 12-18 years |
| Adult | 18-65 years |
| Elderly | >65 years |
| All ages | Variable onset |

### Mapping Relations

| Code | Meaning |
|------|---------|
| E | Exact mapping |
| NTBT | Narrower term, broader term |
| BTNT | Broader term, narrower term |
| ND | Not decided |

---

## Entity Relationships

### Disease to Genes
- **Cardinality:** N:M
- **Description:** Diseases linked to multiple genes
- **Key Fields:** ORPHAcode, GeneSymbol

### Disease to External References
- **Cardinality:** 1:N
- **Description:** Cross-references to OMIM, ICD
- **Key Fields:** ORPHAcode, Source, Reference

### Disease Hierarchy
- **Cardinality:** N:1
- **Description:** Disease classification hierarchy
- **Key Fields:** ORPHAcode, parent_ORPHAcode

### Disease to Epidemiology
- **Cardinality:** 1:N
- **Description:** Prevalence by region
- **Key Fields:** ORPHAcode, PrevalenceGeographic

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ORPHA | Orphanet | Disease identifier |
| ORDO | Orphanet Rare Disease Ontology | OWL format |
| AD | Autosomal Dominant | Inheritance |
| AR | Autosomal Recessive | Inheritance |
| XLR | X-Linked Recessive | Inheritance |
| XLD | X-Linked Dominant | Inheritance |
| LOF | Loss of Function | Mutation type |
| GOF | Gain of Function | Mutation type |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| ICD | International Classification of Diseases | Cross-reference |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| OMIM | MIM number | Gene/phenotype |
| ICD-10 | ICD code | Clinical coding |
| ICD-11 | ICD-11 code | New WHO version |
| UMLS | CUI | Concept mapping |
| MeSH | MeSH ID | NLM vocabulary |
| SNOMED CT | SCTID | Clinical terms |
| MONDO | MONDO ID | Disease ontology |
| HPO | HP ID | Phenotype terms |

---

## Data Quality Notes

1. **Rare Disease Coverage:** 6,528 rare diseases
2. **Gene Associations:** 4,512 linked genes
3. **Expert Curation:** 41-country network
4. **Epidemiology Data:** Prevalence and incidence
5. **Multi-Language:** 8 language translations
6. **Monthly Updates:** Regular data releases
7. **CC BY 4.0:** Free for commercial use
8. **REST API:** Programmatic access available

