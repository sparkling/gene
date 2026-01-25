# FooDB - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | foodb |
| **Name** | FooDB |
| **Total Fields** | 35+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| food_id | integer | 1:1 | Yes | Internal numeric identifier | `234` |
| public_id | string | 1:1 | Yes | Public identifier with FOOD prefix | `FOOD00234` |
| name | string | 1:1 | Yes | Common food name | `Apple` |
| name_scientific | string | 1:1 | No | Latin binomial nomenclature | `Malus domestica` |

### Classification

| Field Name | Data Type | Cardinality | Required | Description | Allowed Values |
|------------|-----------|-------------|----------|-------------|----------------|
| food_group | string | 1:1 | No | High-level category | Herbs and Spices, Fruits, Vegetables, etc. |
| food_subgroup | string | 1:1 | No | Secondary category | Herbs, Berries, Root vegetables |
| food_type | string | 1:1 | No | Classification type | - |

### Compound Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| compound_id | string | 1:N | No | FooDB compound ID | `FDB000001` |
| moldb_formula | string | 1:1 | No | Molecular formula | `C21H21O11+` |
| moldb_smiles | string | 1:1 | No | SMILES notation | `CC(=O)Oc1ccccc1C(=O)O` |
| moldb_inchikey | string | 1:1 | No | InChI key | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` |
| moldb_average_mass | number | 1:1 | No | Average molecular mass (Da) | `449.386` |
| cas_number | string | 1:1 | No | CAS registry number | `7084-24-4` |

### Content Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| orig_content | number | 1:1 | No | Average concentration | `25.3` |
| orig_min | number | 1:1 | No | Minimum concentration | `15.2` |
| orig_max | number | 1:1 | No | Maximum concentration | `42.8` |
| orig_unit | string | 1:1 | No | Unit of measurement | `mg/100g` |
| source_type | string | 1:1 | Yes | Content source type | `Compound`, `Nutrient` |

---

## Enumerations

### Compound Status

| Value | Code | Description |
|-------|------|-------------|
| Quantified | 0 | Detected and quantified |
| Detected | 1 | Detected but not quantified |
| Expected | 2 | Expected but not quantified |
| Predicted | 3 | Predicted to exist |

### Export Status

| Value | Description |
|-------|-------------|
| 1 | Active record |
| 2 | Inactive record |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| FooDB | Food Database | World's largest food constituent database |
| SMILES | Simplified Molecular-Input Line-Entry System | Chemical structure notation |
| InChI | International Chemical Identifier | IUPAC structure identifier |
| CAS | Chemical Abstracts Service | Chemical registry system |
| ITIS | Integrated Taxonomic Information System | Species classification |
| IUPAC | International Union of Pure and Applied Chemistry | Nomenclature authority |
| ppm | Parts Per Million | Concentration unit |

---

## Data Quality Notes

1. **Concentration data** varies in quality based on original literature source
2. **Unit standardization**: Most values normalized to mg/100g or mcg/100g
3. **N:M Relationships**: Foods have many compounds; compounds found in many foods
4. **Null handling**: Empty concentration fields indicate data not available

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
