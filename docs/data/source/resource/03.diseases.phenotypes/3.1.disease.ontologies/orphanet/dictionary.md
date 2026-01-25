# Orphanet - Data Dictionary

## Overview

This data dictionary documents the schema for Orphanet rare disease ontology.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | orphanet |
| **Name** | Orphanet |
| **Parent** | 3.1.disease.ontologies |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Disorder (Disease)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Primary identifier | 558 |
| Name | string | 1:1 | Yes | Preferred disease name | Marfan syndrome |
| ORPHAgroup | enum | 1:1 | Yes | Entity type | Disorder |
| DisorderType | object | 1:1 | Yes | Disease/malformation type | Disease |
| SynonymList | array | 1:N | No | Alternative names | [MFS, Marfan's] |
| ExternalReferenceList | array | 1:N | No | Cross-references | OMIM:154700 |
| TextualInformation | object | 1:1 | No | Definition, description | - |
| AverageAgeOfOnset | array | 1:N | No | Onset age categories | [Childhood, Infancy] |
| AverageAgeOfDeath | array | 1:N | No | Death age categories | - |
| TypeOfInheritanceList | array | 1:N | No | Inheritance patterns | [AD] |

### Gene Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Disease identifier | 558 |
| Gene.Symbol | string | 1:1 | Yes | HGNC symbol | FBN1 |
| Gene.Name | string | 1:1 | Yes | Gene name | fibrillin 1 |
| Gene.HGNC | integer | 1:1 | No | HGNC ID | 3603 |
| Gene.GeneLocus | string | 1:1 | No | Chromosomal location | 15q21.1 |
| DisorderGeneAssociationType | object | 1:1 | Yes | Association type | Disease-causing |
| DisorderGeneAssociationStatus | object | 1:1 | Yes | Validation status | Assessed |
| SourceOfValidation | string | 1:1 | No | Literature source | PMID:1852208 |

### HPO Phenotype Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Disease identifier | 558 |
| HPOId | string | 1:1 | Yes | HPO term ID | HP:0001166 |
| HPOTerm | string | 1:1 | Yes | Phenotype name | Arachnodactyly |
| HPOFrequency | object | 1:1 | No | Frequency classification | Very frequent |
| Diagnostic | boolean | 1:1 | No | Diagnostic criterion | true |

### Epidemiological Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ORPHAcode | integer | 1:1 | Yes | Disease identifier | 558 |
| PrevalenceType | string | 1:1 | No | Point/Birth/Lifetime | Point prevalence |
| PrevalenceClass | string | 1:1 | No | Categorical prevalence | 1-9 / 100 000 |
| ValMoy | float | 1:1 | No | Mean prevalence value | 1.5 |
| PrevalenceGeographic | string | 1:1 | No | Geographic scope | Europe |
| PrevalenceValidationStatus | string | 1:1 | No | Validation status | Validated |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| ORPHAcode | [0-9]+ | 558 | Numeric Orphanet ID |
| Orphanet ID | Orphanet:[0-9]+ | Orphanet:558 | CURIE format |
| ORDO URI | http://www.orpha.net/ORDO/Orphanet_[0-9]+ | Orphanet_558 | Ontology URI |
| OMIM | OMIM:[0-9]{6} | OMIM:154700 | Mendelian disease |
| HGNC | HGNC:[0-9]+ | HGNC:3603 | Gene identifier |
| ICD-10 | [A-Z][0-9]{2}\.[0-9] | Q87.4 | Clinical coding |

---

## Enumerations

### ORPHAgroup Types

| Type | Description |
|------|-------------|
| Disorder | Single clinical entity |
| Group of disorders | Collection of related diseases |
| Subtype of a disorder | Variant of parent disorder |

### Disorder Types

| Type | Description |
|------|-------------|
| Disease | Clinical disease |
| Malformation syndrome | Congenital malformation |
| Morphological anomaly | Structural anomaly |
| Particular clinical situation | Special circumstances |
| Biological anomaly | Biological abnormality |
| Clinical syndrome | Clinical pattern |
| Etiological subtype | Cause-based subtype |
| Histopathological subtype | Pathology-based subtype |
| Clinical subtype | Presentation-based subtype |

### Gene-Disease Association Types

| ID | Name | Description |
|----|------|-------------|
| 1 | Disease-causing germline mutation(s) in | Direct causative |
| 2 | Disease-causing germline mutation(s) (loss of function) | LOF mechanism |
| 3 | Disease-causing germline mutation(s) (gain of function) | GOF mechanism |
| 4 | Major susceptibility factor in | Strong risk factor |
| 5 | Role in phenotype of | Modifier gene |
| 6 | Biomarker tested in | Diagnostic marker |

### Inheritance Types

| ID | Name |
|----|------|
| 1 | Autosomal dominant |
| 2 | Autosomal recessive |
| 3 | X-linked dominant |
| 4 | X-linked recessive |
| 5 | Mitochondrial inheritance |
| 6 | Multigenic/multifactorial |
| 7 | Not applicable |
| 8 | Unknown |

### HPO Frequency Categories

| Category | HPO Term | Percentage |
|----------|----------|------------|
| Obligate | HP:0040280 | 100% |
| Very frequent | HP:0040281 | 80-99% |
| Frequent | HP:0040282 | 30-79% |
| Occasional | HP:0040283 | 5-29% |
| Very rare | HP:0040284 | 1-4% |
| Excluded | HP:0040285 | 0% |

### Prevalence Classes

| Class | Description |
|-------|-------------|
| >1 / 1000 | Common |
| 1-5 / 10 000 | Relatively common |
| 1-9 / 100 000 | Rare |
| <1 / 1 000 000 | Very rare |
| Unknown | Data not available |

---

## Entity Relationships

### Disease to Genes
- **Cardinality:** N:M
- **Description:** Diseases associated with multiple genes
- **Key Fields:** ORPHAcode, Gene.Symbol

### Disease to Phenotypes
- **Cardinality:** N:M
- **Description:** Diseases linked to HPO phenotypes
- **Key Fields:** ORPHAcode, HPOId

### Disease to External References
- **Cardinality:** 1:N
- **Description:** Cross-references to OMIM, ICD, etc.
- **Key Fields:** ORPHAcode, ExternalReferenceList

### Disease Hierarchy
- **Cardinality:** N:1
- **Description:** Disease classification hierarchy
- **Key Fields:** ORPHAcode, parent

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ORPHA | Orphanet | Disease identifier prefix |
| ORDO | Orphanet Rare Disease Ontology | OWL format |
| HPO | Human Phenotype Ontology | Phenotype terms |
| OMIM | Online Mendelian Inheritance in Man | Cross-reference |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| LOF | Loss of Function | Mutation type |
| GOF | Gain of Function | Mutation type |
| AD | Autosomal Dominant | Inheritance |
| AR | Autosomal Recessive | Inheritance |
| XLR | X-linked Recessive | Inheritance |
| XLD | X-linked Dominant | Inheritance |
| HOOM | HPO-Orphanet Ontological Module | Integration |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| OMIM | MIM number | Gene and phenotype |
| ICD-10 | ICD code | WHO classification |
| ICD-11 | ICD-11 code | New WHO version |
| UMLS | CUI | CUI mappings |
| MeSH | MeSH ID | NLM vocabulary |
| SNOMED CT | SCTID | Clinical terms |
| MedDRA | PT code | Drug safety |
| GARD | GARD ID | NIH rare diseases |
| MONDO | MONDO ID | Disease ontology |

---

## Data Quality Notes

1. **Rare Disease Coverage:** 6,528 rare diseases
2. **Gene Associations:** 4,512 linked genes
3. **HPO Annotations:** 100,000+ phenotype links
4. **Expert Curation:** 41-country expert network
5. **Epidemiology:** Prevalence and incidence data
6. **Multi-Language:** 8 language translations
7. **Data Products:** Orphadata XML downloads
8. **CC BY 4.0:** Free for commercial use

