# HPO Phenotype Database - Data Dictionary

## Overview

This data dictionary documents the schema for HPO phenotype annotation database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | hpo |
| **Name** | HPO Phenotype Database |
| **Parent** | 3.2.phenotype.databases |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### HPO Term

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | HPO identifier | HP:0001250 |
| name | string | 1:1 | Yes | Term name | Seizure |
| def | string | 1:1 | No | Definition with source | A seizure is an event... |
| synonym | array | 1:N | No | Alternative names | Epileptic seizure |
| xref | array | 1:N | No | Cross-references | UMLS:C0036572 |
| is_a | array | 1:N | Yes | Parent term relationships | HP:0012638 |
| comment | string | 1:1 | No | Additional notes | - |
| created_by | string | 1:1 | No | Curator ID | peter |
| creation_date | datetime | 1:1 | No | When created | 2008-02-27 |

### Disease Annotation (phenotype.hpoa)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| database_id | string | 1:1 | Yes | Disease database ID | OMIM:154700 |
| disease_name | string | 1:1 | Yes | Disease name | Marfan syndrome |
| qualifier | string | 1:1 | No | NOT if negated | NOT |
| hpo_id | string | 1:1 | Yes | HPO term ID | HP:0001166 |
| reference | string | 1:1 | Yes | Evidence source | OMIM:154700 |
| evidence | enum | 1:1 | Yes | Evidence code | IEA, TAS, PCS |
| onset | string | 1:1 | No | Age of onset | HP:0003577 |
| frequency | string | 1:1 | No | Frequency info | HP:0040283 |
| sex | enum | 1:1 | No | Sex specificity | MALE, FEMALE |
| modifier | string | 1:1 | No | Clinical modifiers | HP:0012828 |
| aspect | enum | 1:1 | Yes | HPO branch | P, I, C, M, H |
| biocuration | string | 1:1 | Yes | Curator and date | HPO:skoehler[2024-01-15] |

### Gene-Phenotype Link

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_id | integer | 1:1 | Yes | NCBI Gene ID | 7157 |
| gene_symbol | string | 1:1 | Yes | HGNC symbol | TP53 |
| hpo_id | string | 1:1 | Yes | HPO term ID | HP:0002664 |
| hpo_name | string | 1:1 | Yes | Term name | Neoplasm |
| frequency | string | 1:1 | No | Frequency info | - |
| disease_id | string | 1:1 | Yes | Source disease | OMIM:191170 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| HP ID | HP:[0-9]{7} | HP:0001250 | HPO term identifier |
| OMIM ID | OMIM:[0-9]{6} | OMIM:154700 | Disease reference |
| ORPHA ID | ORPHA:[0-9]+ | ORPHA:558 | Rare disease reference |
| DECIPHER ID | DECIPHER:[0-9]+ | DECIPHER:16 | CNV database |
| MONDO ID | MONDO:[0-9]{7} | MONDO:0007947 | Disease ontology |
| Gene ID | Integer | 7157 | NCBI Gene ID |

---

## Enumerations

### Evidence Codes

| Code | Name | Description |
|------|------|-------------|
| IEA | Inferred from Electronic Annotation | Automated text mining |
| TAS | Traceable Author Statement | Expert curation |
| PCS | Published Clinical Study | Research evidence |

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

### Onset Categories

| HP ID | Name | Age Range |
|-------|------|-----------|
| HP:0030674 | Antenatal onset | Before birth |
| HP:0003577 | Congenital onset | At birth |
| HP:0003593 | Infantile onset | 0-1 year |
| HP:0011463 | Childhood onset | 1-5 years |
| HP:0003621 | Juvenile onset | 5-15 years |
| HP:0003581 | Adult onset | 16+ years |
| HP:0003584 | Late onset | 40+ years |

### Database Prefixes

| Prefix | Source | Description |
|--------|--------|-------------|
| OMIM: | OMIM | Mendelian diseases |
| ORPHA: | Orphanet | Rare diseases |
| DECIPHER: | DECIPHER | CNV syndromes |
| MONDO: | MONDO | Integrated diseases |

---

## Entity Relationships

### Term to Parents
- **Cardinality:** N:M
- **Description:** Hierarchical is_a relationships
- **Key Fields:** id, is_a

### Disease to Phenotypes
- **Cardinality:** N:M
- **Description:** Disease-phenotype associations
- **Key Fields:** database_id, hpo_id

### Gene to Phenotypes
- **Cardinality:** N:M
- **Description:** Gene-phenotype associations via diseases
- **Key Fields:** gene_id, hpo_id

### Term to Cross-References
- **Cardinality:** 1:N
- **Description:** Mappings to external terminologies
- **Key Fields:** id, xref

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HPO | Human Phenotype Ontology | Database name |
| HPOA | HPO Annotation | Annotation format |
| IEA | Inferred from Electronic Annotation | Evidence code |
| TAS | Traceable Author Statement | Evidence code |
| PCS | Published Clinical Study | Evidence code |
| OMIM | Online Mendelian Inheritance in Man | Disease source |
| ORPHA | Orphanet | Rare disease source |
| DECIPHER | - | CNV database |
| NCBI | National Center for Biotechnology Information | Gene IDs |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| OMIM | MIM number | Disease annotations |
| Orphanet | ORPHA ID | Rare disease annotations |
| MONDO | MONDO ID | Disease ontology |
| UMLS | CUI | Concept mapping |
| SNOMED CT | SCTID | Clinical terms |
| Gene (NCBI) | Gene ID | Gene associations |

---

## Data Quality Notes

1. **Term Coverage:** 13,000+ phenotype terms
2. **Disease Annotations:** 156,000+ annotations
3. **Annotated Diseases:** ~9,500 diseases
4. **Gene Associations:** ~4,300 genes
5. **Evidence Tracking:** IEA/TAS/PCS evidence codes
6. **Frequency Data:** HPO frequency terms
7. **Phenopacket Support:** GA4GH standard
8. **Monthly Updates:** Regular releases

