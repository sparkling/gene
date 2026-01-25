---
id: schema-phenol.explorer
title: "Phenol-Explorer Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, polyphenols, flavonoids, food-composition, bioavailability]
---

# Phenol-Explorer - Polyphenol Database Schema

**Document ID:** SCHEMA-PHENOL-EXPLORER
**Version:** 3.6
**Source Version:** 2024

---

## TL;DR

Phenol-Explorer provides comprehensive data on polyphenol content in foods and their metabolism. The schema links compounds to food sources with quantitative composition data and analytical methods, plus metabolite information for bioavailability studies. Essential for dietary polyphenol research.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Polyphenol Compounds | 500+ | Compound database |
| Foods Covered | 450+ | Food database |
| Composition Data Points | 35,000+ | Food content data |
| Metabolites | 380+ | Metabolite database |
| Pharmacokinetic Studies | 250+ | PK data |
| Publications | 1,500+ | Literature references |

---

## Entity Relationship Overview

```
Compounds (1) ←→ (many) Food_Contents (many) ←→ (1) Foods
     ↓                        ↓                      ↓
 PubChem ID            mg/100g value            Food name

Compounds (1) ←→ (many) Metabolites
                            ↓
                    Phase I/II, Microbial

Compounds (1) ←→ (many) Pharmacokinetics
                            ↓
                    Cmax, Tmax, AUC values

Compounds (1) ←→ (1) Classification
                       ↓
                 Class → Subclass
```

---

## Core Tables/Entities

### compounds

**Description:** Polyphenol compounds with structural information.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| compound_id | integer | Yes | Primary identifier |
| name | string | Yes | Compound name |
| class | string | Yes | Polyphenol class |
| subclass | string | No | Subclass |
| pubchem_cid | integer | No | PubChem compound ID |
| cas_number | string | No | CAS registry number |
| molecular_formula | string | No | Chemical formula |
| molecular_weight | decimal | No | MW in Daltons |
| smiles | string | No | SMILES structure |

### foods

**Description:** Food items with composition data.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| food_id | integer | Yes | Primary identifier |
| food_name | string | Yes | Food name |
| food_group | string | No | Food category |
| scientific_name | string | No | Plant species |
| langual_code | string | No | LanguaL food code |

### food_contents

**Description:** Polyphenol content in foods.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| content_id | integer | Yes | Primary identifier |
| compound_id | integer | Yes | Foreign key to compounds |
| food_id | integer | Yes | Foreign key to foods |
| content_value | decimal | Yes | Amount in food |
| content_unit | string | Yes | mg/100g, mg/100ml |
| min_value | decimal | No | Minimum reported |
| max_value | decimal | No | Maximum reported |
| data_points | integer | No | Number of samples |
| analytical_method | string | No | HPLC, LC-MS, etc. |
| reference_id | integer | No | Literature source |

### metabolites

**Description:** Polyphenol metabolites in humans.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| metabolite_id | integer | Yes | Primary identifier |
| parent_compound_id | integer | Yes | Foreign key to compounds |
| metabolite_name | string | Yes | Metabolite name |
| metabolite_type | string | Yes | Phase I, Phase II, Microbial |
| transformation | string | No | Glucuronidation, etc. |
| biofluid | string | No | plasma, urine, etc. |

### pharmacokinetics

**Description:** Pharmacokinetic parameters from studies.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| pk_id | integer | Yes | Primary identifier |
| compound_id | integer | Yes | Foreign key to compounds |
| metabolite_id | integer | No | Specific metabolite |
| dose | decimal | No | Administered dose |
| dose_unit | string | No | mg, g |
| food_matrix | string | No | Source food/pure |
| cmax | decimal | No | Maximum concentration |
| cmax_unit | string | No | umol/L, ng/mL |
| tmax | decimal | No | Time to Cmax |
| tmax_unit | string | No | hours, minutes |
| auc | decimal | No | Area under curve |
| half_life | decimal | No | Elimination half-life |
| reference_id | integer | No | Study reference |

---

## Polyphenol Classes

| Class | Subclasses |
|-------|------------|
| Flavonoids | Flavonols, Flavones, Flavan-3-ols, Anthocyanins, Flavanones, Isoflavones |
| Phenolic acids | Hydroxybenzoic acids, Hydroxycinnamic acids |
| Stilbenes | Resveratrol and derivatives |
| Lignans | Secoisolariciresinol, Matairesinol |
| Other polyphenols | Curcuminoids, Ellagitannins |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | Web interface with export |
| Alternative | CSV/Excel downloads |
| Encoding | UTF-8 |
| Units | mg/100g fresh weight (standard) |

---

## Sample Record

```json
{
  "compound": {
    "id": 1,
    "name": "Quercetin",
    "class": "Flavonoids",
    "subclass": "Flavonols",
    "pubchem_cid": 5280343,
    "molecular_weight": 302.24
  },
  "food_content": {
    "food": "Onion, red, raw",
    "content": 39.21,
    "unit": "mg/100g",
    "min": 20.3,
    "max": 50.6,
    "method": "HPLC-DAD"
  },
  "metabolites": [
    {
      "name": "Quercetin-3-O-glucuronide",
      "type": "Phase II",
      "transformation": "Glucuronidation"
    }
  ]
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| Polyphenol | Plant secondary metabolite with phenol groups |
| Flavonoid | Six-carbon ring compounds |
| Cmax | Maximum plasma concentration |
| Tmax | Time to reach Cmax |
| AUC | Area Under the Curve (bioavailability) |
| Phase II metabolism | Conjugation reactions (glucuronidation, sulfation) |

---

## References

1. Phenol-Explorer: http://phenol-explorer.eu
2. Rothwell JA, et al. (2013) Database (Oxford). 2013:bat070
3. Neveu V, et al. (2010) Database (Oxford). 2010:bap024
