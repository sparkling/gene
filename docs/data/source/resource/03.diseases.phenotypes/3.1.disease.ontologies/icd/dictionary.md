# ICD - Data Dictionary

## Overview

This data dictionary documents the schema for ICD (International Classification of Diseases).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | icd |
| **Name** | ICD |
| **Parent** | 3.1.disease.ontologies |
| **Total Fields** | 15+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### ICD Entity

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| code | string | 1:1 | Yes | ICD code | E11.9 |
| title | string | 1:1 | Yes | Disease name | Type 2 diabetes mellitus |
| definition | string | 1:1 | No | Clinical definition | A metabolic disorder... |
| chapter | string | 1:1 | Yes | Chapter number | 5 |
| block | string | 1:1 | Yes | Block range | E10-E14 |
| parent_code | string | 1:1 | No | Parent code | E11 |
| includes | array | 1:N | No | Inclusion terms | - |
| excludes | array | 1:N | No | Exclusion terms | - |
| coding_hint | string | 1:1 | No | Coding guidance | Use additional code... |

### ICD-11 Extension Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| entity_id | string | 1:1 | Yes | ICD-11 entity ID | 5A11 |
| foundation_uri | string | 1:1 | Yes | Foundation URI | http://id.who.int/... |
| linearization | string | 1:1 | Yes | Linearization type | MMS |
| postcoordination | object | 1:1 | No | Extension axes | severity, laterality |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| ICD-10-CM | [A-Z][0-9]{2}(\.[0-9]{1,4})? | E11.9 | US clinical modification |
| ICD-10 | [A-Z][0-9]{2}(\.[0-9])? | I21.0 | WHO version |
| ICD-11 | [A-Z0-9]{4,6} | 5A11 | New WHO version |
| Block Code | [A-Z][0-9]{2}-[A-Z][0-9]{2} | E10-E14 | Code range |
| Chapter | [0-9]{1,2} | 5 | Chapter number |

---

## Enumerations

### ICD-10 Chapters

| Chapter | Code Range | Description |
|---------|------------|-------------|
| 1 | A00-B99 | Infectious diseases |
| 2 | C00-D48 | Neoplasms |
| 3 | D50-D89 | Blood diseases |
| 4 | E00-E90 | Endocrine disorders |
| 5 | F00-F99 | Mental disorders |
| 6 | G00-G99 | Nervous system |
| 7 | H00-H59 | Eye diseases |
| 8 | H60-H95 | Ear diseases |
| 9 | I00-I99 | Circulatory system |
| 10 | J00-J99 | Respiratory system |
| 11 | K00-K93 | Digestive system |
| 12 | L00-L99 | Skin diseases |
| 13 | M00-M99 | Musculoskeletal |
| 14 | N00-N99 | Genitourinary |
| 15 | O00-O99 | Pregnancy |
| 16 | P00-P96 | Perinatal conditions |
| 17 | Q00-Q99 | Congenital anomalies |
| 18 | R00-R99 | Symptoms, signs |
| 19 | S00-T98 | Injury, poisoning |
| 20 | V01-Y98 | External causes |
| 21 | Z00-Z99 | Health services |
| 22 | U00-U99 | Special purposes |

### Code Specificity Levels

| Level | Description | Example |
|-------|-------------|---------|
| Category | 3-character | E11 |
| Subcategory | 4-character | E11.6 |
| Code | 5-7 character | E11.65 |

### ICD-11 Extension Axes

| Axis | Description | Examples |
|------|-------------|----------|
| Severity | Disease severity | Mild, Moderate, Severe |
| Laterality | Body side | Left, Right, Bilateral |
| Temporality | Time aspects | Acute, Chronic |
| Etiology | Causal factors | Infectious, Neoplastic |
| Anatomy | Body location | Specific sites |
| Histopathology | Tissue type | Specific histology |

### Coding Conventions

| Convention | Description |
|------------|-------------|
| Use additional code | Add secondary diagnosis |
| Code first | Primary condition |
| Code also | Related condition |
| Excludes1 | Never coded together |
| Excludes2 | May code both |
| Includes | Inclusion terms |

---

## Entity Relationships

### Code to Parent
- **Cardinality:** N:1
- **Description:** Hierarchical code structure
- **Key Fields:** code, parent_code

### Code to Block
- **Cardinality:** N:1
- **Description:** Codes belong to blocks
- **Key Fields:** code, block

### Code to Chapter
- **Cardinality:** N:1
- **Description:** Blocks belong to chapters
- **Key Fields:** block, chapter

### ICD-10 to ICD-11
- **Cardinality:** N:M
- **Description:** Version mappings
- **Key Fields:** icd10_code, icd11_entity_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ICD | International Classification of Diseases | WHO standard |
| ICD-10 | ICD 10th Revision | Current version |
| ICD-10-CM | ICD-10 Clinical Modification | US adaptation |
| ICD-11 | ICD 11th Revision | New version (2022) |
| WHO | World Health Organization | Publisher |
| MMS | Mortality and Morbidity Statistics | ICD-11 linearization |
| NOS | Not Otherwise Specified | Coding term |
| NEC | Not Elsewhere Classified | Coding term |
| CM | Clinical Modification | US adaptation |
| PCS | Procedure Coding System | Procedure codes |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| SNOMED CT | SCTID | Clinical terms |
| MeSH | MeSH ID | Medical vocabulary |
| MONDO | MONDO ID | Disease mapping |
| Orphanet | ORPHA ID | Rare diseases |
| UMLS | CUI | Unified vocabulary |

---

## Data Quality Notes

1. **ICD-10 Coverage:** ~70,000 codes
2. **ICD-11 Coverage:** ~80,000 entities
3. **Global Adoption:** 194 countries
4. **Clinical Coding:** Primary use for billing/statistics
5. **WHO Standard:** International reporting
6. **Version Transition:** ICD-11 adoption ongoing
7. **API Access:** ICD API with registration
8. **Annual Updates:** Regular code additions

