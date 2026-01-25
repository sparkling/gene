# USDA FoodData Central - Data Dictionary

## Overview

This data dictionary documents the schema for USDA FoodData Central.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | usda.fooddata |
| **Name** | USDA FoodData Central |
| **Parent** | 2.4.food.compounds.nutrients |
| **Total Fields** | 35+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Food Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| fdc_id | integer | 1:1 | Yes | FoodData Central ID | 167512 |
| data_type | string | 1:1 | Yes | Data source type | Foundation, SR Legacy, Branded |
| description | string | 1:1 | Yes | Food description | Apple, raw |
| food_class | string | 1:1 | No | Food classification | Branded, Survey |
| publication_date | date | 1:1 | No | Data publication date | 2024-01-01 |
| scientific_name | string | 1:1 | No | Species name | Malus domestica |

### Nutrient Definition

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| nutrient_id | integer | 1:1 | Yes | Nutrient identifier | 1008 |
| name | string | 1:1 | Yes | Nutrient name | Energy |
| unit_name | string | 1:1 | Yes | Unit of measurement | kcal, g, mg, mcg |
| nutrient_nbr | string | 1:1 | No | USDA nutrient number | 208 |
| rank | integer | 1:1 | No | Display order | 300 |

### Food Nutrient Values

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | Primary identifier | 12345 |
| amount | decimal | 1:1 | Yes | Nutrient amount per 100g | 52.0 |
| data_points | integer | 1:1 | No | Number of samples | 15 |
| derivation_id | integer | 1:1 | No | Calculation method | 1 |
| min | decimal | 1:1 | No | Minimum value | 48.0 |
| max | decimal | 1:1 | No | Maximum value | 56.0 |
| median | decimal | 1:1 | No | Median value | 52.0 |
| loq | decimal | 1:1 | No | Limit of quantitation | 0.01 |
| footnote | string | 1:1 | No | Additional notes | Calculated value |

### Food Portions

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | Primary identifier | 101 |
| sequence_number | integer | 1:1 | No | Portion order | 1 |
| amount | decimal | 1:1 | Yes | Portion amount | 1.0 |
| portion_description | string | 1:1 | No | Portion description | 1 medium (3" dia) |
| modifier | string | 1:1 | No | Additional description | without skin |
| gram_weight | decimal | 1:1 | Yes | Weight in grams | 182.0 |

### Branded Foods

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| brand_owner | string | 1:1 | No | Brand owner | General Mills |
| brand_name | string | 1:1 | No | Product brand | Cheerios |
| gtin_upc | string | 1:1 | No | GTIN/UPC barcode | 016000275287 |
| ingredients | text | 1:1 | No | Ingredient list | Whole grain oats, sugar... |
| serving_size | decimal | 1:1 | No | Serving size | 28.0 |
| serving_size_unit | string | 1:1 | No | Serving unit | g |
| household_serving | string | 1:1 | No | Household measure | 1 cup |
| branded_food_category | string | 1:1 | No | Product category | Breakfast cereals |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| FDC ID | Integer | 167512 | Primary food identifier |
| NDB Number | 5 digits | 09003 | Legacy nutrient database |
| Nutrient ID | Integer | 1008 | Nutrient identifier |
| GTIN/UPC | 12-14 digits | 016000275287 | Product barcode |
| Derivation ID | Integer | 1 | Calculation method |

---

## Enumerations

### Data Types

| Type | Description | Records |
|------|-------------|---------|
| Foundation | Comprehensive analytical profiles | 2,100+ |
| SR Legacy | Standard Reference (historical) | 7,800+ |
| Branded | Industry-submitted products | 360,000+ |
| FNDDS | Survey foods (NHANES) | 7,000+ |
| Experimental | Research data | Varies |

### Nutrient Categories

| Category | Examples |
|----------|----------|
| Proximates | Energy, protein, fat, carbohydrates, water |
| Minerals | Calcium, iron, magnesium, zinc, potassium |
| Vitamins | A, C, D, E, K, B-complex |
| Lipids | Fatty acids, cholesterol, trans fats |
| Amino acids | Essential and non-essential |
| Other | Fiber, sugars, alcohol, caffeine |

### Key Nutrients (by ID)

| ID | Name | Unit |
|----|------|------|
| 1008 | Energy | kcal |
| 1003 | Protein | g |
| 1004 | Total lipid (fat) | g |
| 1005 | Carbohydrate | g |
| 1079 | Fiber, total dietary | g |
| 1087 | Calcium, Ca | mg |
| 1089 | Iron, Fe | mg |
| 1092 | Potassium, K | mg |
| 1093 | Sodium, Na | mg |
| 1106 | Vitamin A, RAE | mcg |
| 1162 | Vitamin C | mg |

### Derivation Methods

| Method | Description |
|--------|-------------|
| Analytical | Direct laboratory analysis |
| Calculated | Computed from other values |
| Imputed | Estimated from similar foods |
| Assumed | Standard assumption applied |
| Manufacturer | Manufacturer-provided data |

### Food Classes

| Class | Description |
|-------|-------------|
| FinalFood | Consumable food item |
| SubSample | Analytical subsample |
| Composite | Blended sample |
| Branded | Commercial product |
| Survey | NHANES dietary survey |

---

## Entity Relationships

### Food to Nutrients
- **Cardinality:** 1:N
- **Description:** Each food has multiple nutrient values
- **Key Fields:** fdc_id, nutrient_id

### Food to Portions
- **Cardinality:** 1:N
- **Description:** Each food has multiple serving sizes
- **Key Fields:** fdc_id, portion_id

### Branded Food to Food
- **Cardinality:** 1:1
- **Description:** Branded products linked to base food entry
- **Key Fields:** fdc_id

### Foundation Food to Food
- **Cardinality:** 1:1
- **Description:** Foundation entries linked to base food
- **Key Fields:** fdc_id, ndb_number

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| FDC | FoodData Central | Database name |
| NDB | Nutrient Data Bank | Legacy identifier |
| SR | Standard Reference | Historical data type |
| FNDDS | Food and Nutrient Database for Dietary Studies | Survey data |
| GTIN | Global Trade Item Number | Product identifier |
| UPC | Universal Product Code | Barcode |
| USDA | United States Department of Agriculture | Agency |
| ARS | Agricultural Research Service | USDA division |
| RAE | Retinol Activity Equivalents | Vitamin A unit |
| DFE | Dietary Folate Equivalents | Folate unit |
| NHANES | National Health and Nutrition Examination Survey | Survey |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| LanguaL | Food codes | Food description |
| FooDB | FooDB ID | Food composition |
| Phenol-Explorer | Internal | Polyphenol content |
| HMDB | HMDB ID | Metabolites |
| FDA | GTIN/UPC | Product registration |

---

## Data Quality Notes

1. **Authoritative Source:** Official USDA food composition data
2. **Multiple Data Types:** Foundation (analytical), SR Legacy, Branded, Survey
3. **Comprehensive Nutrients:** 150+ nutrients tracked per food
4. **Portion Data:** Standardized serving sizes with gram weights
5. **Branded Coverage:** 360,000+ commercial products
6. **API Access:** RESTful API with JSON responses
7. **Update Frequency:** Continuous rolling updates
8. **Public Domain:** Data freely available for any use
