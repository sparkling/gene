# MeSH - Data Dictionary

## Overview

This data dictionary documents the schema for MeSH (Medical Subject Headings).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | mesh |
| **Name** | MeSH |
| **Parent** | 3.1.disease.ontologies |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Descriptor Record

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| DescriptorUI | string | 1:1 | Yes | MeSH unique identifier | D003920 |
| DescriptorName | string | 1:1 | Yes | Preferred term name | Diabetes Mellitus |
| TreeNumber | array | 1:N | Yes | Hierarchical positions | C18.452.394.750 |
| ScopeNote | string | 1:1 | No | Definition/scope | A heterogeneous group... |
| HistoryNote | string | 1:1 | No | Term history | Year introduced, changes |
| PublicMeSHNote | string | 1:1 | No | Usage guidance | See also specific types |

### Term Records

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Term | string | 1:1 | Yes | Entry term | Type 2 Diabetes |
| TermUI | string | 1:1 | Yes | Term identifier | T000001 |
| IsPermuted | boolean | 1:1 | No | Permuted term | false |
| LexicalTag | string | 1:1 | No | Lexical type | NON |
| RecordPreferred | boolean | 1:1 | No | Preferred term | true |

### Qualifiers

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| QualifierUI | string | 1:1 | Yes | Qualifier ID | Q000188 |
| QualifierName | string | 1:1 | Yes | Qualifier name | drug therapy |
| Abbreviation | string | 1:1 | Yes | Short form | DT |
| TreeNumber | string | 1:1 | Yes | Qualifier tree | Y02.407 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Descriptor UI | D[0-9]{6} | D003920 | Main heading ID |
| Supplemental UI | C[0-9]{6} | C537775 | Supplemental concept |
| Tree Number | [A-Z][0-9]{2}(\.[0-9]{3})+ | C18.452.394.750 | Hierarchy position |
| Qualifier UI | Q[0-9]{6} | Q000188 | Subheading ID |
| Term UI | T[0-9]{6} | T000001 | Entry term ID |

---

## Enumerations

### Major Tree Categories

| Category | Code | Description |
|----------|------|-------------|
| Anatomy | A | Body structures |
| Organisms | B | Living organisms |
| Diseases | C | Disease conditions |
| Chemicals | D | Drugs and chemicals |
| Analytical | E | Techniques, equipment |
| Psychiatry | F | Psychology, behavior |
| Phenomena | G | Biological phenomena |
| Disciplines | H | Fields and occupations |
| Anthropology | I | Culture, society |
| Technology | J | Food, technology |
| Humanities | K | Humanities |
| Information | L | Information science |
| Named Groups | M | Person categories |
| Healthcare | N | Health services |
| Geographicals | Z | Geographic locations |

### Disease Tree (C Category)

| SubCategory | Code | Description |
|-------------|------|-------------|
| Infections | C01 | Bacterial, viral, etc. |
| Neoplasms | C04 | Cancer |
| Musculoskeletal | C05 | Bone, joint diseases |
| Digestive | C06 | GI disorders |
| Stomatognathic | C07 | Oral diseases |
| Respiratory | C08 | Lung diseases |
| Otorhinolaryngologic | C09 | ENT diseases |
| Nervous | C10 | Neurological |
| Eye | C11 | Ophthalmologic |
| Male Urogenital | C12 | Male GU diseases |
| Female Urogenital | C13 | Female GU diseases |
| Cardiovascular | C14 | Heart, vessels |
| Hemic | C15 | Blood diseases |
| Congenital | C16 | Birth defects |
| Skin | C17 | Dermatologic |
| Nutritional | C18 | Metabolic |
| Endocrine | C19 | Hormonal |
| Immune | C20 | Immunological |
| Environmental | C21 | Environmental |
| Animal | C22 | Veterinary |
| Pathological | C23 | General pathology |
| Signs | C26 | Symptoms and signs |

### Record Types

| Type | Description |
|------|-------------|
| Descriptor | Main headings |
| Qualifier | Subheadings |
| Supplementary Concept | Chemicals, drugs, diseases |
| Scope Note | Definitions |
| Entry Term | Synonyms |

### Supplementary Concept Types

| Type | Description |
|------|-------------|
| Chemical | Drug/chemical names |
| Disease | Rare disease names |
| Protocol | Treatment protocols |

---

## Entity Relationships

### Descriptor to Tree
- **Cardinality:** 1:N
- **Description:** One descriptor can appear in multiple trees
- **Key Fields:** DescriptorUI, TreeNumber

### Descriptor to Terms
- **Cardinality:** 1:N
- **Description:** Descriptors have multiple entry terms
- **Key Fields:** DescriptorUI, TermUI

### Descriptor to Qualifiers
- **Cardinality:** N:M
- **Description:** Allowable descriptor-qualifier pairs
- **Key Fields:** DescriptorUI, QualifierUI

### Tree Hierarchy
- **Cardinality:** N:1
- **Description:** Tree number implies parent
- **Key Fields:** TreeNumber (parent derived)

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| MeSH | Medical Subject Headings | Database name |
| NLM | National Library of Medicine | Publisher |
| UI | Unique Identifier | ID format |
| SCR | Supplementary Concept Record | Additional concepts |
| PA | Pharmacological Action | Drug class |
| EC | Entry Combination | Term combination |
| RN | Registry Number | CAS number |
| NM | Non-MeSH | Supplementary concept |
| HM | Heading Mapped To | Mapping target |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UMLS | CUI | Concept mapping |
| SNOMED CT | SCTID | Clinical terms |
| ICD-10 | ICD code | Disease codes |
| RxNorm | RXCUI | Drug concepts |
| ChEBI | ChEBI ID | Chemical entities |
| DrugBank | DrugBank ID | Drug information |

---

## Data Quality Notes

1. **Descriptor Coverage:** ~30,000 descriptors
2. **Disease Descriptors:** ~5,000 in C tree
3. **Supplementary Concepts:** ~280,000 chemicals/diseases
4. **PubMed Indexing:** Primary vocabulary for literature
5. **Annual Updates:** New terms added yearly
6. **Multi-Format:** XML, ASCII, RDF/N-Triples, JSON-LD
7. **SPARQL Endpoint:** Linked data access
8. **Public Domain:** Free to use

