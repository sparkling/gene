---
id: schemas-usda-fooddata-central
title: "USDA FoodData Central Schema Documentation"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, usda, fooddata-central, nutrition, food-composition, api]
---

**Parent:** [Schema Documentation](./_index.md)

# USDA FoodData Central - Schema Documentation

**Document ID:** SCHEMA-USDA-FDC
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://fdc.nal.usda.gov/

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "12345" |
| `name` | string | Entity name | "Apple" |
| `type` | string | Record type | "food" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Overview

FoodData Central is the USDA's integrated food and nutrient data system, providing expanded nutrient profile data for foods. It merges data from five distinct data sources: Foundation Foods, SR Legacy, Survey Foods (FNDDS), Branded Foods, and Experimental Foods.

---

## API Specification

### Base URL
```
https://api.nal.usda.gov/fdc/v1/
```

### Authentication
- **Type:** API Key (query parameter)
- **Parameter:** `api_key`
- **Obtain Key:** https://api.data.gov/signup/
- **Demo Key:** `DEMO_KEY` (for testing, rate limited)

### Rate Limits
- Default: 1,000 requests per hour per IP address
- Exceeding limit: Temporary 1-hour block
- Higher limits: Contact FoodData Central

---

## Endpoints

### 1. Single Food Retrieval
```
GET /v1/food/{fdcId}
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| fdcId | integer | Yes | FoodData Central ID |
| format | string | No | `abridged` or `full` (default: full) |
| nutrients | array | No | Up to 25 nutrient numbers to filter |

### 2. Multiple Foods Retrieval
```
GET/POST /v1/foods
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| fdcIds | array | Yes | Array of FDC IDs (max 20) |
| format | string | No | `abridged` or `full` |
| nutrients | array | No | Nutrient numbers to filter |

### 3. Foods List (Paginated)
```
GET/POST /v1/foods/list
```

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| dataType | array | No | All | `Branded`, `Foundation`, `Survey`, `SR Legacy` |
| pageSize | integer | No | 50 | Results per page (1-200) |
| pageNumber | integer | No | 1 | Page to retrieve |
| sortBy | string | No | - | `dataType`, `description`, `fdcId`, `publishedDate` |
| sortOrder | string | No | asc | `asc` or `desc` |

### 4. Foods Search
```
GET/POST /v1/foods/search
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| query | string | Yes | Search terms (supports operators) |
| dataType | array | No | Filter by data types |
| brandOwner | string | No | Filter by brand owner |
| pageSize | integer | No | Results per page (1-200) |
| pageNumber | integer | No | Page number |
| sortBy | string | No | Sort field |
| sortOrder | string | No | `asc` or `desc` |

---

## Data Types

### Foundation Foods
- Nutrient values from chemical analysis
- Sample-level data with analytical methods
- Scientific names, food components
- Published biannually (April, December)

### SR Legacy
- USDA Standard Reference database (final 2018 release)
- Static dataset, no longer updated
- Historical reference data

### Survey Foods (FNDDS)
- Food and Nutrient Database for Dietary Studies
- Updated biennially with NHANES
- Portion weights and descriptions

### Branded Foods
- Labeled nutrient data from manufacturers
- GS1 GTIN identifiers
- Updated monthly

### Experimental Foods
- Limited release research data
- Special analytical studies

---

## Core Data Schemas

### Food Schema

```json
{
  "fdcId": 534358,
  "dataType": "Foundation",
  "description": "Apples, raw, with skin",
  "foodClass": "FinalFood",
  "publicationDate": "2019-04-01",
  "scientificName": "Malus domestica",
  "foodCategory": {
    "id": 9,
    "code": "0900",
    "description": "Fruits and Fruit Juices"
  },
  "foodNutrients": [
    {
      "nutrient": {
        "id": 1003,
        "number": "203",
        "name": "Protein",
        "unitName": "g"
      },
      "amount": 0.26,
      "dataPoints": 8,
      "min": 0.19,
      "max": 0.33,
      "median": 0.25
    }
  ],
  "ndbNumber": "09003",
  "isHistoricalReference": false
}
```

---

## Key Nutrient Numbers

| Number | Name | Unit | Description |
|--------|------|------|-------------|
| 203 | Protein | g | Total protein |
| 204 | Total lipid (fat) | g | Total fat |
| 205 | Carbohydrate, by difference | g | Total carbohydrates |
| 208 | Energy | kcal | Calories |
| 269 | Sugars, total | g | Total sugars |
| 291 | Fiber, total dietary | g | Dietary fiber |
| 301 | Calcium, Ca | mg | Calcium |
| 303 | Iron, Fe | mg | Iron |
| 306 | Potassium, K | mg | Potassium |
| 307 | Sodium, Na | mg | Sodium |
| 318 | Vitamin A, IU | IU | Vitamin A |
| 401 | Vitamin C | mg | Ascorbic acid |
| 601 | Cholesterol | mg | Cholesterol |
| 606 | Fatty acids, saturated | g | Saturated fat |
| 645 | Fatty acids, monounsaturated | g | MUFA |
| 646 | Fatty acids, polyunsaturated | g | PUFA |

---

## License

**License:** CC0 (Public Domain)
**Attribution:** Not required, but appreciated
**Commercial Use:** Allowed without restriction

---

## Data Set Size

| Metric | Value |
|--------|-------|
| **Total Foods** | 20,900+ documented foods |
| **Branded Foods** | 400,000+ (Supplemental data) |
| **Foundation Foods** | 1,500+ with analytical data |
| **SR Legacy Foods** | 8,700+ (archived 2018) |
| **Survey Foods (FNDDS)** | 6,700+ (dietary survey data) |
| **Experimental Foods** | Limited, varies by release |
| **Nutrient Types** | 150+ measured nutrients |
| **Food Categories** | 20+ major groups |
| **Countries Represented** | United States primarily |
| **Data Points per Nutrient** | 1-50+ samples (Foundation Foods) |
| **API Requests Daily** | 100,000+ |
| **CSV Export Size** | ~500 MB (full dataset) |
| **Branded Foods Dataset** | 400,000+ products in supplemental files |
| **Batch Download Frequency** | Updated monthly (Branded), biannually (Foundation) |
| **Update Cycle** | April and December for Foundation; Monthly for Branded |

---

## Sample Data

### Example Record
```json
{
  "fdc_id": 173944,
  "description": "Broccoli, raw",
  "data_type": "SR Legacy",
  "nutrient_id": 1104,
  "nutrient_name": "Fiber, total dietary",
  "amount": 2.4,
  "unit": "g"
}
```

### Sample Query Result
| fdc_id | description | nutrient_name | amount | unit |
|--------|-----------|---------------|--------|------|
| 173944 | Broccoli, raw | Fiber, total dietary | 2.4 | g |
| 169145 | Banana, raw | Potassium, K | 358 | mg |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 400,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | JSON (API) |
| Alternative | CSV |
| Encoding | UTF-8 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| USDA FoodData Central | HTTP | https://fdc.nal.usda.gov/download-datasets.html |
| FDC API | REST | https://fdc.nal.usda.gov/api |
| Bulk Data | HTTP/FTP | See USDA FDC resources |

**Access Requirements:** Open access, registration required for API key (free)

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| USDA FoodData Central | Public Domain (USDA) | Yes |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `fdcId` | FoodData Central unique food identifier | `534358` |
| `dataType` | Category of food data source | `Foundation` |
| `description` | Human-readable food name | `Apples, raw, with skin` |
| `foodCategory` | Hierarchical classification with code and description | `Fruits and Fruit Juices` |
| `nutrient.number` | USDA nutrient identification number | `203` (Protein) |
| `amount` | Measured nutrient quantity per 100g | `0.26` |
| `ndbNumber` | Legacy USDA National Nutrient Database number | `09003` |
| `scientificName` | Taxonomic name for food source | `Malus domestica` |
| `brandOwner` | Company owning branded food product | `General Mills` |
| `gtinUpc` | GS1 barcode identifier for branded foods | `016000275287` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Foundation Foods | Analytically derived nutrient values with sample data | dataType |
| SR Legacy | Final 2018 USDA Standard Reference database | dataType |
| Survey Foods (FNDDS) | Data supporting dietary surveys like NHANES | dataType |
| Branded Foods | Label-derived data from food manufacturers | dataType |
| Nutrient Profile | Complete set of nutrient values for a food | foodNutrients |
| Data Points | Number of samples analyzed for a nutrient value | Statistical precision |
| Abridged Format | Minimal response with essential fields only | API format parameter |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| FDC | FoodData Central | USDA integrated food data system |
| USDA | United States Department of Agriculture | Database maintainer |
| ARS | Agricultural Research Service | USDA research division |
| FNDDS | Food and Nutrient Database for Dietary Studies | Survey food data source |
| NHANES | National Health and Nutrition Examination Survey | Uses FNDDS data |
| SR | Standard Reference | Legacy nutrient database |
| NDB | National Nutrient Database | Former name of SR |
| GTIN | Global Trade Item Number | Barcode standard |
| GS1 | Global Standards 1 | Barcode organization |
| MUFA | Monounsaturated Fatty Acids | Nutrient type |
| PUFA | Polyunsaturated Fatty Acids | Nutrient type |
| IU | International Unit | Vitamin measurement unit |
| CC0 | Creative Commons Zero | Public domain dedication |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
