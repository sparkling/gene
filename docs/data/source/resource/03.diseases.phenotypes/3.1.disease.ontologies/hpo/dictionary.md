# HPO - Data Dictionary

## Overview

This data dictionary documents the schema for HPO (Human Phenotype Ontology).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | hpo |
| **Name** | HPO |
| **Parent** | 3.1.disease.ontologies |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### HPO Term

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | HPO identifier | HP:0001250 |
| name | string | 1:1 | Yes | Phenotype name | Seizure |
| def | string | 1:1 | No | Definition with source | A seizure is an event... |
| synonym | array | 1:N | No | Alternative names with types | Epileptic seizure |
| xref | array | 1:N | No | Cross-references | UMLS:C0036572 |
| is_a | array | 1:N | Yes | Parent term relationships | HP:0012638 |
| comment | string | 1:1 | No | Additional notes | - |
| created_by | string | 1:1 | No | Curator ID | peter |
| creation_date | datetime | 1:1 | No | Creation timestamp | 2008-02-27T02:20:00Z |

### Disease Annotation (phenotype.hpoa)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| database_id | string | 1:1 | Yes | Disease ID | OMIM:154700 |
| disease_name | string | 1:1 | Yes | Disease name | Marfan syndrome |
| qualifier | string | 1:1 | No | NOT if phenotype excluded | NOT |
| hpo_id | string | 1:1 | Yes | HPO term ID | HP:0001166 |
| reference | string | 1:1 | Yes | Evidence source | OMIM:154700 |
| evidence | enum | 1:1 | Yes | Evidence code | TAS |
| onset | string | 1:1 | No | Age of onset HPO term | HP:0003577 |
| frequency | string | 1:1 | No | Frequency term or fraction | HP:0040283 |
| sex | enum | 1:1 | No | Sex specificity | MALE, FEMALE |
| modifier | string | 1:1 | No | Clinical modifiers | HP:0012828 |
| aspect | enum | 1:1 | Yes | HPO branch | P |
| biocuration | string | 1:1 | Yes | Curator and date | HPO:skoehler[2024-01-15] |

### Gene-Phenotype Link

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_id | integer | 1:1 | Yes | NCBI Gene ID | 7157 |
| gene_symbol | string | 1:1 | Yes | HGNC symbol | TP53 |
| hpo_id | string | 1:1 | Yes | HPO term ID | HP:0002664 |
| hpo_name | string | 1:1 | Yes | Term name | Neoplasm |
| frequency | string | 1:1 | No | Frequency info | HP:0040283 |
| disease_id | string | 1:1 | Yes | Source disease | OMIM:191170 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| HP ID | HP:[0-9]{7} | HP:0001250 | HPO term identifier |
| OMIM | OMIM:[0-9]{6} | OMIM:154700 | Disease database |
| ORPHA | ORPHA:[0-9]+ | ORPHA:558 | Rare disease |
| DECIPHER | DECIPHER:[0-9]+ | DECIPHER:16 | CNV database |
| UMLS CUI | UMLS:C[0-9]+ | UMLS:C0036572 | Unified vocabulary |
| SNOMED CT | SNOMEDCT_US:[0-9]+ | SNOMEDCT_US:91175000 | Clinical terms |
| MeSH | MSH:D[0-9]+ | MSH:D012640 | Medical subjects |

---

## Enumerations

### Evidence Codes

| Code | Name | Description |
|------|------|-------------|
| IEA | Inferred from Electronic Annotation | Automated text mining from OMIM |
| TAS | Traceable Author Statement | Expert curation with citation |
| PCS | Published Clinical Study | From clinical research |

### Aspect Codes

| Code | Branch | Description |
|------|--------|-------------|
| P | Phenotypic abnormality | Clinical features |
| I | Inheritance | Mode of inheritance |
| C | Clinical course | Onset, progression |
| M | Clinical modifier | Severity, laterality |
| H | Past medical history | Historical features |

### Frequency Terms

| HP ID | Name | Percentage |
|-------|------|------------|
| HP:0040280 | Obligate | 100% |
| HP:0040281 | Very frequent | 80-99% |
| HP:0040282 | Frequent | 30-79% |
| HP:0040283 | Occasional | 5-29% |
| HP:0040284 | Very rare | 1-4% |
| HP:0040285 | Excluded | 0% |

### Synonym Types

| Type | Description |
|------|-------------|
| EXACT | Exact synonym |
| RELATED | Related term |
| NARROW | More specific |
| BROAD | More general |
| layperson | Layperson-friendly term |
| abbreviation | Abbreviated form |
| plural_form | Plural variant |
| uk_spelling | UK spelling variant |

### Top-Level Categories

| HP ID | Name | Description |
|-------|------|-------------|
| HP:0000118 | Phenotypic abnormality | Root for clinical features |
| HP:0000005 | Mode of inheritance | Inheritance patterns |
| HP:0012823 | Clinical modifier | Severity, modifiers |
| HP:0040279 | Frequency | Occurrence frequency |
| HP:0032443 | Past medical history | Historical features |

---

## Entity Relationships

### Term to Parents
- **Cardinality:** N:M
- **Description:** Hierarchical is_a relationships
- **Key Fields:** id, is_a

### Term to Disease Annotations
- **Cardinality:** 1:N
- **Description:** Links phenotypes to diseases
- **Key Fields:** hpo_id, database_id

### Term to Gene Associations
- **Cardinality:** N:M
- **Description:** Phenotypes associated with genes
- **Key Fields:** hpo_id, gene_id

### Term to Cross-References
- **Cardinality:** 1:N
- **Description:** Mappings to external terminologies
- **Key Fields:** id, xref

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HPO | Human Phenotype Ontology | Database name |
| HP | Human Phenotype | Term prefix |
| HPOA | HPO Annotation | Annotation file format |
| IEA | Inferred from Electronic Annotation | Evidence code |
| TAS | Traceable Author Statement | Evidence code |
| PCS | Published Clinical Study | Evidence code |
| OMIM | Online Mendelian Inheritance in Man | Disease database |
| ORPHA | Orphanet | Rare disease database |
| SNOMED CT | Systematized Nomenclature of Medicine Clinical Terms | Clinical terminology |
| UMLS | Unified Medical Language System | Vocabulary |
| MeSH | Medical Subject Headings | NLM vocabulary |
| GA4GH | Global Alliance for Genomics and Health | Standards body |
| OBO | Open Biological Ontologies | Ontology format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UMLS | CUI | Concept mapping |
| SNOMED CT | SCTID | Clinical terms |
| MeSH | MeSH ID | Medical vocabulary |
| MedDRA | PT code | Drug safety |
| ICD-10 | ICD code | Clinical coding |
| MONDO | MONDO ID | Disease mapping |
| Uberon | Uberon ID | Anatomy |
| ChEBI | ChEBI ID | Chemical/metabolite |

---

## Data Quality Notes

1. **Term Coverage:** 13,000+ phenotype terms
2. **Disease Annotations:** 156,000+ disease-phenotype links
3. **Annotated Diseases:** ~9,500 rare diseases
4. **Gene Associations:** ~4,300 genes
5. **Multi-Language:** 10+ language translations
6. **Evidence Tracking:** IEA, TAS, PCS evidence codes
7. **Phenopacket Standard:** Foundation for GA4GH data exchange
8. **Monthly Updates:** Regular releases with new terms

