---
id: schema-rxnorm
title: "RxNorm Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, terminology, drug-names, normalization, umls]
---

# RxNorm - Normalized Drug Nomenclature Schema

**Document ID:** SCHEMA-RXNORM
**Version:** 1.0
**Source Version:** Current (monthly updates)

---

## TL;DR

RxNorm provides standardized drug names linking clinical drugs to a common identifier (RXCUI). The schema organizes drug concepts hierarchically from ingredients to branded products with relationships between concepts, enabling medication interoperability across healthcare systems.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Drug Concepts (RXCUIs) | 300,000+ | Concept table |
| Source Vocabularies | 12+ | SAB mappings |
| Ingredient Concepts | 14,000+ | IN term type |
| Branded Products | 45,000+ | SBD term type |
| Relationships | 2,000,000+ | RXNREL table |

---

## Entity Relationship Overview

```
Ingredients (IN) → Semantic Clinical Drug (SCD) → Branded Drug (SBD)
       ↓                    ↓                          ↓
    RXCUI               RXCUI + DF                RXCUI + Brand

                    ↓           ↓
             Clinical Pack    Branded Pack
               (GPCK)          (BPCK)

All concepts linked via:
RXNREL (relationships)
RXNSAT (attributes)
RXNSTY (semantic types)
```

---

## Core Tables/Entities

### RXNCONSO

**Description:** RxNorm concept names and sources.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| RXCUI | string | Yes | RxNorm Concept Unique Identifier |
| LAT | string | Yes | Language (ENG) |
| TS | string | No | Term status |
| LUI | string | Yes | Lexical Unique Identifier |
| STT | string | No | String type |
| SUI | string | Yes | String Unique Identifier |
| ISPREF | string | No | Preferred atom indicator |
| RXAUI | string | Yes | RxNorm Atom Unique Identifier |
| SAUI | string | No | Source Atom Unique Identifier |
| SCUI | string | No | Source Concept Unique Identifier |
| SDUI | string | No | Source Descriptor Unique Identifier |
| SAB | string | Yes | Source Abbreviation |
| TTY | string | Yes | Term Type (IN, SCD, SBD, etc.) |
| CODE | string | No | Source code |
| STR | string | Yes | String (drug name) |
| SRL | string | No | Source restriction level |
| SUPPRESS | string | No | Suppressible flag |
| CVF | string | No | Content view flag |

### RXNREL

**Description:** Relationships between RxNorm concepts.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| RXCUI1 | string | Yes | First concept RXCUI |
| RXAUI1 | string | No | First atom identifier |
| STYPE1 | string | Yes | Source type 1 |
| REL | string | Yes | Relationship type |
| RXCUI2 | string | Yes | Second concept RXCUI |
| RXAUI2 | string | No | Second atom identifier |
| STYPE2 | string | Yes | Source type 2 |
| RELA | string | No | Relationship attribute |
| RUI | string | Yes | Relationship identifier |
| SRUI | string | No | Source relationship ID |
| SAB | string | Yes | Source abbreviation |
| SL | string | No | Source of relationship |
| RG | string | No | Relationship group |
| DIR | string | No | Direction flag |
| SUPPRESS | string | No | Suppressible flag |
| CVF | string | No | Content view flag |

### RXNSAT

**Description:** Attributes of RxNorm concepts.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| RXCUI | string | Yes | Concept identifier |
| LUI | string | No | Lexical identifier |
| SUI | string | No | String identifier |
| RXAUI | string | Yes | Atom identifier |
| STYPE | string | Yes | Source type |
| CODE | string | No | Attribute code |
| ATUI | string | Yes | Attribute identifier |
| SATUI | string | No | Source attribute ID |
| ATN | string | Yes | Attribute name |
| SAB | string | Yes | Source abbreviation |
| ATV | string | Yes | Attribute value |
| SUPPRESS | string | No | Suppressible flag |
| CVF | string | No | Content view flag |

---

## Term Types (TTY)

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

---

## Relationship Types (REL/RELA)

| RELA | Description |
|------|-------------|
| consists_of | Pack consists of components |
| contains | Drug contains ingredient |
| dose_form_of | Is dose form of |
| has_dose_form | Has specific dose form |
| has_ingredient | Has active ingredient |
| has_quantified_form | Has quantified version |
| ingredient_of | Is ingredient of |
| isa | Is-a relationship |
| tradename_of | Is brand name of |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /rxcui/{rxcui} | GET | Get concept by RXCUI |
| /drugs | GET | Search drugs by name |
| /rxcui/{rxcui}/related | GET | Get related concepts |
| /rxcui/{rxcui}/allrelated | GET | Get all relationships |
| /rxcui/{rxcui}/ndcs | GET | Get NDC codes |
| /approximateTerm | GET | Approximate string match |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | Rich Release Format (RRF) |
| Alternative | XML via API |
| Encoding | UTF-8 |
| Delimiter | Pipe (|) separated |

---

## Sample Record

```
# RXNCONSO sample
2670|ENG|S|L0014479|PF|S0042379|Y|A0043984||2670||RXNORM|SCD|2670|Acetaminophen 500 MG Oral Tablet|||

# RXNREL sample
2670|A0043984|SCUI|RN|161|A0000089|SCUI|has_ingredient|R12345|||RXNORM|||N||
```

```json
{
  "rxcui": "2670",
  "name": "Acetaminophen 500 MG Oral Tablet",
  "tty": "SCD",
  "ingredients": [
    {
      "rxcui": "161",
      "name": "Acetaminophen"
    }
  ],
  "dose_form": "Oral Tablet",
  "strength": "500 MG"
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| RXCUI | RxNorm Concept Unique Identifier |
| SCD | Semantic Clinical Drug (generic) |
| SBD | Semantic Branded Drug (brand) |
| TTY | Term Type |
| SAB | Source Abbreviation |
| RRF | Rich Release Format (UMLS standard) |

---

## References

1. RxNav: https://rxnav.nlm.nih.gov/
2. RxNorm Technical Documentation: https://www.nlm.nih.gov/research/umls/rxnorm/docs/
3. API Documentation: https://rxnav.nlm.nih.gov/RxNormAPIs.html
