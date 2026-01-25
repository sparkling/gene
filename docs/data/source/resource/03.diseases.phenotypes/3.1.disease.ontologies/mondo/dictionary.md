# MONDO - Data Dictionary

## Overview

This data dictionary documents the schema for MONDO (Monarch Disease Ontology).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | mondo |
| **Name** | MONDO |
| **Parent** | 3.1.disease.ontologies |
| **Total Fields** | 15+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Disease Term

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | MONDO identifier | MONDO:0005015 |
| name | string | 1:1 | Yes | Disease name | diabetes mellitus |
| definition | string | 1:1 | No | Textual definition | A metabolic disorder... |
| synonyms | array | 1:N | No | Alternative names | [DM, sugar diabetes] |
| xrefs | array | 1:N | No | Cross-references | [OMIM:154700, Orphanet:558] |
| is_a | array | 1:N | Yes | Parent term relationships | MONDO:0000001 |
| subset | array | 1:N | No | Subset membership | mondo_rare |
| created_by | string | 1:1 | No | Curator ID | orcid:0000-0002-... |
| creation_date | datetime | 1:1 | No | When created | 2020-01-15 |

### Equivalence Axioms

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| exactMatch | array | 1:N | No | 1:1 equivalences | OMIM:154700 |
| closeMatch | array | 1:N | No | Near equivalences | DOID:14323 |
| narrowMatch | array | 1:N | No | More specific mappings | NCIt:C34807 |
| broadMatch | array | 1:N | No | More general mappings | MESH:D008382 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| MONDO ID | MONDO:[0-9]{7} | MONDO:0005015 | Primary identifier |
| OMIM | OMIM:[0-9]{6} | OMIM:154700 | Mendelian disease |
| Orphanet | Orphanet:[0-9]+ | Orphanet:558 | Rare disease |
| DOID | DOID:[0-9]+ | DOID:14323 | Disease ontology |
| NCIt | NCIT:C[0-9]+ | NCIT:C34807 | NCI cancer terms |
| MeSH | MESH:D[0-9]{6} | MESH:D008382 | Medical subjects |
| EFO | EFO_[0-9]{7} | EFO_0000685 | Experimental factors |

---

## Enumerations

### Mapping Predicates

| Predicate | Description |
|-----------|-------------|
| exactMatch | 1:1 equivalence, safe for data propagation |
| closeMatch | Very similar concept |
| narrowMatch | Maps to more specific concept |
| broadMatch | Maps to more general concept |
| relatedMatch | Related but different concept |

### Subsets

| Subset | Description |
|--------|-------------|
| mondo_rare | Rare diseases (< 1:2000 prevalence) |
| mondo_core | Core disease terms |
| gard_rare | GARD rare diseases |
| ordo_disease | Orphanet rare diseases |
| ordo_group | Orphanet disease groups |

### Disease Categories

| Category | Description | Examples |
|----------|-------------|----------|
| Human disease | Human disease terms | MONDO:0005015 |
| Cancer | Neoplastic diseases | 4,728 classes |
| Mendelian | Single-gene disorders | 11,639 conditions |
| Rare disease | Low prevalence | 15,901 conditions |
| Infectious | Infectious diseases | Various |
| Autoimmune | Immune-mediated | Various |

### Record Types

| Type | Description |
|------|-------------|
| disease | Clinical disease entity |
| disorder | Disorder (synonym) |
| syndrome | Syndrome entity |
| condition | Medical condition |

---

## Entity Relationships

### Term to Parents
- **Cardinality:** N:1
- **Description:** Hierarchical is_a relationships
- **Key Fields:** id, is_a

### Term to Cross-References
- **Cardinality:** N:M
- **Description:** Mappings to external disease databases
- **Key Fields:** id, xrefs

### Term to Equivalences
- **Cardinality:** 1:1
- **Description:** Exact match equivalences validated by OWL reasoning
- **Key Fields:** id, exactMatch

### Term to Subsets
- **Cardinality:** N:M
- **Description:** Membership in filtering categories
- **Key Fields:** id, subset

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| MONDO | Monarch Disease Ontology | Database name |
| OMIM | Online Mendelian Inheritance in Man | Cross-reference |
| DOID | Disease Ontology Identifier | Alternative ontology |
| NCIt | NCI Thesaurus | Cancer terminology |
| EFO | Experimental Factor Ontology | Trait ontology |
| MeSH | Medical Subject Headings | NLM vocabulary |
| OWL | Web Ontology Language | Ontology format |
| OBO | Open Biological Ontologies | Ontology format |
| CURIE | Compact URI Expression | Identifier format |
| GARD | Genetic and Rare Diseases | NIH program |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| OMIM | MIM number | Mendelian diseases |
| Orphanet | ORPHA code | Rare diseases |
| DOID | DOID ID | Disease ontology |
| NCIt | NCIt code | Cancer terms |
| MeSH | MeSH ID | Medical vocabulary |
| EFO | EFO ID | Traits/diseases |
| UMLS | CUI | Unified vocabulary |
| ICD-10 | ICD code | Clinical coding |

---

## Data Quality Notes

1. **Disease Coverage:** 25,938 total diseases
2. **Human Diseases:** 22,977 human-specific terms
3. **Equivalence Axioms:** OWL-reasoned 1:1 mappings
4. **Cross-References:** 129,914 database xrefs
5. **Pivot Ontology:** Designed for data integration
6. **Logic-Based:** OWL reasoning validates mappings
7. **Monthly Updates:** Regular releases
8. **CC BY 4.0:** Free for commercial use with attribution

