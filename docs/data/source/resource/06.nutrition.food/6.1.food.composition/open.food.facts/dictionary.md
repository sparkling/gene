# Open Food Facts - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | open-food-facts |
| **Name** | Open Food Facts |
| **Total Fields** | 100+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| code | string | 1:1 | Yes | Product barcode (EAN-13) | `3017620422003` |
| product_name | string | 1:1 | No | Product name | `Nutella` |
| brands | string | 1:1 | No | Brand name(s), comma-separated | `Ferrero` |
| generic_name | string | 1:1 | No | Generic product description | `Hazelnut spread` |

### Nutrition

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| energy-kcal_100g | number | 1:1 | No | Energy per 100g | `539` |
| fat_100g | number | 1:1 | No | Total fat (g) | `30.7` |
| saturated-fat_100g | number | 1:1 | No | Saturated fat (g) | `10.6` |
| carbohydrates_100g | number | 1:1 | No | Total carbs (g) | `57.5` |
| sugars_100g | number | 1:1 | No | Total sugars (g) | `56.3` |
| fiber_100g | number | 1:1 | No | Dietary fiber (g) | `3.4` |
| proteins_100g | number | 1:1 | No | Protein (g) | `6.3` |
| salt_100g | number | 1:1 | No | Salt (g) | `0.107` |

### Scores

| Field Name | Data Type | Cardinality | Required | Description | Allowed Values |
|------------|-----------|-------------|----------|-------------|----------------|
| nutriscore_grade | string | 1:1 | No | Nutri-Score grade | a, b, c, d, e |
| nova_group | integer | 1:1 | No | NOVA classification | 1, 2, 3, 4 |
| ecoscore_grade | string | 1:1 | No | Eco-Score grade | a, b, c, d, e |

---

## Enumerations

### NOVA Classification

| Value | Name | Description |
|-------|------|-------------|
| 1 | Unprocessed | Unprocessed or minimally processed foods |
| 2 | Processed Ingredients | Processed culinary ingredients |
| 3 | Processed Foods | Processed foods |
| 4 | Ultra-processed | Ultra-processed food and drink products |

### Nutri-Score Grades

| Grade | Score Range | Description |
|-------|-------------|-------------|
| a | -15 to -1 | Best nutritional quality |
| b | 0 to 2 | Good nutritional quality |
| c | 3 to 10 | Moderate nutritional quality |
| d | 11 to 18 | Poor nutritional quality |
| e | 19 to 40 | Worst nutritional quality |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| OFF | Open Food Facts | Crowdsourced food database |
| EAN | European Article Number | Barcode standard (EAN-13) |
| ODbL | Open Database License | Database license |
| NOVA | - | Not an acronym, food classification |
| JSONL | JSON Lines | Newline-delimited JSON |
| SFA | Saturated Fatty Acids | Fat classification |
| MUFA | Monounsaturated Fatty Acids | Fat classification |
| PUFA | Polyunsaturated Fatty Acids | Fat classification |

---

## Data Quality Notes

1. **Crowdsourced data**: Quality varies; some fields may be incomplete
2. **Nutrient basis**: All values per 100g or 100ml
3. **Tag format**: Tags use `language:value` format (e.g., `en:organic`)
4. **Images**: Multiple images available (front, nutrition, ingredients)

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
