---
id: schema-usda.fooddata
title: "USDA FoodData Central Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, food-composition, nutrients, usda, nutrition]
---

# USDA FoodData Central Schema

**Document ID:** SCHEMA-USDA-FOODDATA
**Version:** 2.0
**Source Version:** Current (continuously updated)

---

## TL;DR

USDA FoodData Central provides authoritative food composition data integrating multiple USDA data types (Foundation, SR Legacy, Branded, Survey). The schema organizes foods with detailed nutrient profiles, portions, and metadata, supporting dietary assessment, food labeling, and nutrition research.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Total Foods | 380,000+ | All data types |
| Foundation Foods | 2,100+ | Comprehensive profiles |
| SR Legacy Foods | 7,800+ | Historical data |
| Branded Foods | 360,000+ | Industry data |
| Nutrients Tracked | 150+ | Nutrient database |
| Update Frequency | Continuous | Rolling updates |

---

## Entity Relationship Overview

```
Foods (1) ←→ (many) Food_Nutrients (many) ←→ (1) Nutrients
  ↓                       ↓                        ↓
FDC ID               Amount/100g              Nutrient ID

Foods (1) ←→ (many) Food_Portions
                         ↓
                  Serving sizes

Foods (1) ←→ (many) Food_Components
                         ↓
                  Ingredients, additives

Branded_Foods (1) ←→ (1) Foods
      ↓
  GTIN/UPC, Brand
```

---

## Core Tables/Entities

### foods

**Description:** Core food items with identification.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| fdc_id | integer | Yes | FoodData Central ID |
| data_type | string | Yes | Foundation, SR Legacy, etc. |
| description | string | Yes | Food description |
| food_class | string | No | Branded, Survey, etc. |
| publication_date | date | No | Data publication date |
| scientific_name | string | No | Species name (if applicable) |

### nutrients

**Description:** Nutrient definitions.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| nutrient_id | integer | Yes | Nutrient identifier |
| name | string | Yes | Nutrient name |
| unit_name | string | Yes | g, mg, mcg, kcal, IU |
| nutrient_nbr | string | No | USDA nutrient number |
| rank | integer | No | Display order |

### food_nutrients

**Description:** Nutrient values for foods.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| id | integer | Yes | Primary identifier |
| fdc_id | integer | Yes | Foreign key to foods |
| nutrient_id | integer | Yes | Foreign key to nutrients |
| amount | decimal | Yes | Nutrient amount |
| data_points | integer | No | Number of samples |
| derivation_id | integer | No | Calculation method |
| min | decimal | No | Minimum value |
| max | decimal | No | Maximum value |
| median | decimal | No | Median value |
| loq | decimal | No | Limit of quantitation |
| footnote | string | No | Additional notes |

### food_portions

**Description:** Serving size information.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| id | integer | Yes | Primary identifier |
| fdc_id | integer | Yes | Foreign key to foods |
| sequence_number | integer | No | Portion order |
| amount | decimal | Yes | Portion amount |
| measure_unit_id | integer | Yes | Unit type |
| portion_description | string | No | "1 cup", "1 medium", etc. |
| modifier | string | No | Additional description |
| gram_weight | decimal | Yes | Weight in grams |

### branded_foods

**Description:** Industry-submitted branded food products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| fdc_id | integer | Yes | Foreign key to foods |
| brand_owner | string | No | Brand name |
| brand_name | string | No | Product brand |
| gtin_upc | string | No | GTIN/UPC barcode |
| ingredients | text | No | Ingredient list |
| serving_size | decimal | No | Serving size |
| serving_size_unit | string | No | g, ml, etc. |
| household_serving | string | No | Household measure |
| branded_food_category | string | No | Product category |

### foundation_foods

**Description:** Comprehensive analytical data.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| fdc_id | integer | Yes | Foreign key to foods |
| ndb_number | string | No | Legacy NDB number |
| footnote | text | No | Additional notes |

---

## Data Types

| Type | Description | Records |
|------|-------------|---------|
| Foundation | Comprehensive analytical profiles | 2,100+ |
| SR Legacy | Standard Reference (historical) | 7,800+ |
| Branded | Industry-submitted products | 360,000+ |
| FNDDS | Survey foods (NHANES) | 7,000+ |
| Experimental | Research data | Varies |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /v1/foods | POST | Search foods |
| /v1/food/{fdcId} | GET | Get food by ID |
| /v1/foods/list | GET | List foods |
| /v1/foods/search | GET | Search by query |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | JSON (API) |
| Alternative | CSV downloads |
| Encoding | UTF-8 |
| Units | Per 100g edible portion |

---

## Sample Record

```json
{
  "fdcId": 167512,
  "dataType": "Foundation",
  "description": "Apple, raw",
  "publicationDate": "2024-01-01",
  "foodNutrients": [
    {
      "nutrient": {
        "id": 1008,
        "name": "Energy",
        "unitName": "kcal"
      },
      "amount": 52.0
    },
    {
      "nutrient": {
        "id": 1003,
        "name": "Protein",
        "unitName": "g"
      },
      "amount": 0.26
    },
    {
      "nutrient": {
        "id": 1004,
        "name": "Total lipid (fat)",
        "unitName": "g"
      },
      "amount": 0.17
    },
    {
      "nutrient": {
        "id": 1005,
        "name": "Carbohydrate, by difference",
        "unitName": "g"
      },
      "amount": 13.81
    }
  ],
  "foodPortions": [
    {
      "portionDescription": "1 medium (3\" dia)",
      "gramWeight": 182.0
    }
  ]
}
```

---

## Nutrient Categories

| Category | Examples |
|----------|----------|
| Proximates | Energy, protein, fat, carbs, water |
| Minerals | Calcium, iron, magnesium, zinc |
| Vitamins | A, C, D, E, K, B-vitamins |
| Lipids | Fatty acids, cholesterol |
| Amino acids | Essential and non-essential |
| Other | Fiber, sugars, alcohol |

---

## Glossary

| Term | Definition |
|------|------------|
| FDC ID | FoodData Central identifier |
| NDB Number | Legacy USDA nutrient database number |
| Foundation Food | Comprehensive analytical profile |
| SR Legacy | Standard Reference historical data |
| GTIN/UPC | Product barcode identifier |
| Derivation | Method used to calculate nutrient value |

---

## References

1. USDA FoodData Central: https://fdc.nal.usda.gov
2. API Documentation: https://fdc.nal.usda.gov/api-guide.html
3. Data Downloads: https://fdc.nal.usda.gov/download-datasets.html
