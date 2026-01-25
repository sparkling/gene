# RxNorm - Data Dictionary

## Overview

This data dictionary documents the schema for RxNorm normalized drug nomenclature.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | rxnorm |
| **Name** | RxNorm |
| **Parent** | 2.2.pharmaceuticals |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Concept Information (RXNCONSO)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| RXCUI | string | 1:1 | Yes | RxNorm Concept Unique Identifier | 2670 |
| STR | string | 1:1 | Yes | Drug name string | Acetaminophen 500 MG Oral Tablet |
| TTY | string | 1:1 | Yes | Term type | IN, SCD, SBD |
| SAB | string | 1:1 | Yes | Source abbreviation | RXNORM |
| CODE | string | 1:1 | No | Source code | 2670 |
| RXAUI | string | 1:1 | Yes | RxNorm Atom Unique Identifier | A0043984 |
| LUI | string | 1:1 | Yes | Lexical Unique Identifier | L0014479 |
| SUI | string | 1:1 | Yes | String Unique Identifier | S0042379 |
| ISPREF | string | 1:1 | No | Preferred atom indicator | Y/N |
| LAT | string | 1:1 | Yes | Language | ENG |

### Relationship Information (RXNREL)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| RXCUI1 | string | 1:1 | Yes | First concept RXCUI | 2670 |
| RXCUI2 | string | 1:1 | Yes | Second concept RXCUI | 161 |
| REL | string | 1:1 | Yes | Relationship type | RN, RO, RB |
| RELA | string | 1:1 | No | Relationship attribute | has_ingredient |
| RUI | string | 1:1 | Yes | Relationship identifier | R12345 |
| SAB | string | 1:1 | Yes | Source abbreviation | RXNORM |

### Attribute Information (RXNSAT)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| RXCUI | string | 1:1 | Yes | Concept identifier | 2670 |
| ATN | string | 1:1 | Yes | Attribute name | NDC |
| ATV | string | 1:1 | Yes | Attribute value | 00363-0109-01 |
| ATUI | string | 1:1 | Yes | Attribute identifier | AT12345 |
| SAB | string | 1:1 | Yes | Source abbreviation | RXNORM |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| RXCUI | Integer | 2670 | Concept identifier |
| RXAUI | A + digits | A0043984 | Atom identifier |
| RUI | R + digits | R12345 | Relationship ID |
| ATUI | AT + digits | AT12345 | Attribute ID |
| NDC | 10-11 digits | 00363-0109-01 | National Drug Code |

---

## Enumerations

### Term Types (TTY) - Drug Hierarchy

| TTY | Description | Example |
|-----|-------------|---------|
| IN | Ingredient | Acetaminophen |
| MIN | Multiple Ingredients | Acetaminophen / Codeine |
| PIN | Precise Ingredient | Acetaminophen (substance) |
| BN | Brand Name | Tylenol |
| DF | Dose Form | Oral Tablet |
| DFG | Dose Form Group | Pills |
| SCDC | Semantic Clinical Drug Component | Acetaminophen 500 MG |
| SCDF | Semantic Clinical Drug Form | Acetaminophen Oral Tablet |
| SCDG | Semantic Clinical Drug Group | Acetaminophen Oral Product |
| SCD | Semantic Clinical Drug | Acetaminophen 500 MG Oral Tablet |
| SBDC | Semantic Branded Drug Component | Tylenol 500 MG |
| SBDF | Semantic Branded Drug Form | Tylenol Oral Tablet |
| SBDG | Semantic Branded Drug Group | Tylenol Oral Product |
| SBD | Semantic Branded Drug | Tylenol 500 MG Oral Tablet |
| GPCK | Generic Pack | Acetaminophen Pack |
| BPCK | Branded Pack | Tylenol Pack |

### Relationship Attributes (RELA)

| RELA | Description |
|------|-------------|
| consists_of | Pack consists of components |
| contains | Drug contains ingredient |
| dose_form_of | Is dose form of |
| has_dose_form | Has specific dose form |
| has_ingredient | Has active ingredient |
| has_quantified_form | Has quantified version |
| ingredient_of | Is ingredient of |
| isa | Is-a hierarchical relationship |
| tradename_of | Is brand name of |
| has_tradename | Has brand name |

### Relationship Types (REL)

| REL | Description |
|-----|-------------|
| RN | Narrower relationship |
| RO | Other related |
| RB | Broader relationship |
| PAR | Parent |
| CHD | Child |
| SIB | Sibling |
| SY | Synonym |

### Source Abbreviations (SAB)

| SAB | Description |
|-----|-------------|
| RXNORM | RxNorm vocabulary |
| MMSL | Multum MediSource Lexicon |
| NDDF | First Databank NDDF Plus |
| VANDF | VA National Drug File |
| MTHSPL | DailyMed SPL |
| MMX | Micromedex |
| GS | Gold Standard Drug Database |
| DRUGBANK | DrugBank |

---

## Entity Relationships

### Ingredient to Clinical Drug (SCD)
- **Cardinality:** 1:N
- **Description:** Ingredients used in multiple clinical drugs
- **Key Fields:** RXCUI (IN), RXCUI (SCD), RELA=has_ingredient

### SCD to Branded Drug (SBD)
- **Cardinality:** 1:N
- **Description:** Generic drugs with multiple brand versions
- **Key Fields:** RXCUI (SCD), RXCUI (SBD), RELA=tradename_of

### Drug to NDC
- **Cardinality:** 1:N
- **Description:** Drug concepts mapped to NDC codes
- **Key Fields:** RXCUI, ATN=NDC

### Ingredient Hierarchy
- **Cardinality:** N:1
- **Description:** Multiple ingredients under one category
- **Key Fields:** RXCUI (MIN), RXCUI (IN), RELA=has_ingredient

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| RXCUI | RxNorm Concept Unique Identifier | Primary identifier |
| SCD | Semantic Clinical Drug | Generic drug representation |
| SBD | Semantic Branded Drug | Brand drug representation |
| IN | Ingredient | Active substance |
| BN | Brand Name | Tradename |
| DF | Dose Form | Formulation type |
| TTY | Term Type | Concept category |
| SAB | Source Abbreviation | Vocabulary source |
| RELA | Relationship Attribute | Relationship type |
| RRF | Rich Release Format | UMLS file format |
| UMLS | Unified Medical Language System | NLM terminology system |
| NDC | National Drug Code | Package identifier |
| NLM | National Library of Medicine | RxNorm maintainer |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| DailyMed | Set ID | Drug labels |
| Orange Book | Application # | Approval status |
| DrugBank | DrugBank ID | Drug information |
| FDA NDC | NDC | Package codes |
| SNOMED CT | SNOMED ID | Clinical terminology |
| ATC | ATC Code | Drug classification |
| MeSH | MeSH ID | Medical terminology |
| VANDF | VUID | VA drug vocabulary |

---

## Data Quality Notes

1. **Normalization:** All drug names normalized to standard forms
2. **Monthly Updates:** New RxNorm release each month
3. **Hierarchy:** Semantic types enable navigation from ingredient to branded pack
4. **Multi-Source:** Integrates 12+ source vocabularies
5. **NDC Mapping:** Links to historical and current NDC codes
6. **UMLS Format:** Uses Rich Release Format (RRF) files
7. **API Access:** RESTful API via RxNav for real-time queries
8. **License:** UMLS license required for bulk download
