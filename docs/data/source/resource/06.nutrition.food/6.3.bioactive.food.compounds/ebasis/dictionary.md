# eBASIS - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | ebasis |
| **Name** | Bioactive Substances in Food Information System |
| **Total Fields** | 40+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| composition_id | integer | 1:1 | Yes | Record identifier | `12345` |
| food_id | integer | 1:1 | Yes | Food item ID | `101` |
| component_id | integer | 1:1 | Yes | Compound ID | `1234` |
| food_name | string | 1:1 | Yes | Food name | `Apple, raw, with skin` |
| component_name | string | 1:1 | Yes | Compound name | `Quercetin` |

### Classification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| component_group | string | 1:1 | Yes | Major category | `Flavonoids` |
| component_subgroup | string | 1:1 | No | Subcategory | `Flavonols` |
| food_group | string | 1:1 | No | Food category | `Fruits` |

### Composition Values

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| value | number | 1:1 | Yes | Concentration | `4.01` |
| unit | string | 1:1 | Yes | Measurement unit | `mg/100g` |
| value_type | string | 1:1 | No | Value type | `mean` |
| min_value | number | 1:1 | No | Minimum observed | `2.10` |
| max_value | number | 1:1 | No | Maximum observed | `6.45` |
| n_samples | integer | 1:1 | No | Sample count | `12` |

### Quality Assessment

| Field Name | Data Type | Cardinality | Required | Description | Allowed Values |
|------------|-----------|-------------|----------|-------------|----------------|
| overall_score | integer | 1:1 | No | Quality rating | 1-5 |
| analytical_method_score | integer | 1:1 | No | Method quality | 1-5 |

---

## Enumerations

### Component Groups

| Value | Examples |
|-------|----------|
| Flavonoids | Quercetin, Kaempferol, Catechin |
| Carotenoids | Beta-carotene, Lutein, Lycopene |
| Phenolic Acids | Gallic acid, Caffeic acid |
| Stilbenes | Resveratrol |
| Lignans | Secoisolariciresinol |
| Phytosterols | Beta-sitosterol |
| Glucosinolates | Glucoraphanin, Sinigrin |
| Alkaloids | Caffeine, Theobromine |

### Quality Score Scale

| Score | Description |
|-------|-------------|
| 5 | Excellent - Full documentation, validated method |
| 4 | Good - Good documentation, appropriate method |
| 3 | Acceptable - Adequate documentation |
| 2 | Poor - Limited documentation |
| 1 | Very Poor - Minimal information |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| eBASIS | Bioactive Substances in Food Information System | Database name |
| EuroFIR | European Food Information Resource | Network/standards |
| CAS | Chemical Abstracts Service | Identifier system |
| LanguaL | Langua aLimentaria | Food thesaurus |
| HPLC | High-Performance Liquid Chromatography | Analytical method |
| DAD | Diode Array Detector | Detection method |
| GC-MS | Gas Chromatography-Mass Spectrometry | Analytical method |

---

## Data Quality Notes

1. **Quality-evaluated data**: All records scored using EuroFIR criteria
2. **Literature-sourced**: Values from peer-reviewed publications
3. **Multiple sources**: Mean values aggregate multiple studies
4. **Access restricted**: Registration required; commercial use needs license

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
