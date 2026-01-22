# USDA FoodData Central - Schema Documentation

**Document ID:** SCHEMA-USDA-FDC
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://fdc.nal.usda.gov/

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
- Nutrient conversion factors

### Survey Foods (FNDDS)
- Food and Nutrient Database for Dietary Studies
- WWEIA food codes and categories
- Input food mappings
- Updated annually

### Branded Foods
- Data from food industry submissions
- GTIN/UPC codes, brand information
- Ingredient lists, serving sizes
- Published biannually

---

## JSON Schemas

### Search Response
```json
{
  "totalHits": 20939,
  "currentPage": 1,
  "totalPages": 10470,
  "pageList": [1, 2, 3, 4, 5],
  "foodSearchCriteria": {
    "query": "chicken",
    "pageNumber": 1,
    "pageSize": 2,
    "numberOfResultsPerPage": 2,
    "requireAllWords": false
  },
  "foods": [
    {
      "fdcId": 2091022,
      "description": "CHICKEN",
      "dataType": "Branded",
      "gtinUpc": "855203003043",
      "publishedDate": "2021-03-19",
      "brandOwner": "Keystone",
      "brandName": "KEYSTONE",
      "ingredients": "CHICKEN",
      "marketCountry": "United States",
      "foodCategory": "Canned Meat",
      "modifiedDate": "2021-01-29",
      "dataSource": "LI",
      "packageWeight": "28 oz/794g",
      "servingSizeUnit": "g",
      "servingSize": 56,
      "tradeChannels": ["NO_TRADE_CHANNEL"],
      "allHighlightFields": "<b>Chicken</b>",
      "score": 866.0668,
      "foodNutrients": [
        {
          "nutrientId": 1003,
          "nutrientName": "Protein",
          "nutrientNumber": "203",
          "unitName": "G",
          "derivationCode": "LCCD",
          "derivationDescription": "Calculated from a daily value percentage per serving size measure",
          "derivationId": 75,
          "value": 21.4,
          "foodNutrientSourceId": 9,
          "foodNutrientSourceCode": "12",
          "foodNutrientSourceDescription": "Manufacturer's analytical; partial documentation",
          "rank": 600,
          "indentLevel": 1,
          "foodNutrientId": 22988265,
          "percentDailyValue": 43
        }
      ]
    }
  ],
  "aggregations": {
    "dataType": {
      "Branded": 20093,
      "Survey (FNDDS)": 374,
      "SR Legacy": 267,
      "Foundation": 205
    }
  }
}
```

### Foundation Food Item
```json
{
  "fdcId": 2346407,
  "description": "Cabbage, raw",
  "dataType": "Foundation",
  "publicationDate": "2022-10-28",
  "foodClass": "FinalFood",
  "scientificName": "Brassica oleracea var. capitata",
  "foodCategory": {
    "id": 11,
    "code": "1100",
    "description": "Vegetables and Vegetable Products"
  },
  "foodNutrients": [
    {
      "id": 30198742,
      "nutrient": {
        "id": 1051,
        "number": "255",
        "name": "Water",
        "rank": 100,
        "unitName": "g"
      },
      "type": "FoodNutrient",
      "amount": 91.87,
      "dataPoints": 8,
      "min": 90.08,
      "max": 93.1,
      "median": 91.9,
      "foodNutrientDerivation": {
        "id": 1,
        "code": "A",
        "description": "Analytical",
        "foodNutrientSource": {
          "id": 1,
          "code": "1",
          "description": "Analytical or derived from analytical"
        }
      },
      "nutrientAnalysisDetails": {
        "subSampleId": 1,
        "amount": 91.87,
        "nutrientId": 1051,
        "labMethodDescription": "Gravimetric",
        "labMethodLink": "AOAC 934.01",
        "nutrientAcquisitionDetails": [
          {
            "sampleUnitId": 123,
            "purchaseDate": "2022-04-15",
            "storeCity": "Blacksburg",
            "storeState": "VA"
          }
        ]
      }
    }
  ],
  "foodComponents": [
    {
      "id": 1234,
      "name": "Seeds",
      "gramWeight": 2.5,
      "percentWeight": 1.2,
      "isRefuse": false
    }
  ],
  "foodPortions": [
    {
      "id": 567,
      "amount": 1,
      "gramWeight": 70,
      "portionDescription": "cup, chopped",
      "measureUnit": {
        "id": 9999,
        "name": "cup",
        "abbreviation": "cup"
      }
    }
  ],
  "inputFoods": [
    {
      "id": 789,
      "foodDescription": "Raw cabbage, green",
      "inputFood": {
        "fdcId": 170397,
        "description": "Cabbage, raw"
      }
    }
  ]
}
```

### Branded Food Item
```json
{
  "fdcId": 534358,
  "description": "NUT 'N BERRY MIX",
  "dataType": "Branded",
  "publicationDate": "2019-04-01",
  "brandOwner": "Kar Nut Products Company",
  "brandName": "KAR'S",
  "gtinUpc": "077034085228",
  "ingredients": "PEANUTS, RAISINS, DRIED CRANBERRIES (CRANBERRIES, SUGAR, SUNFLOWER OIL), SUNFLOWER KERNELS, ALMONDS",
  "marketCountry": "United States",
  "servingSize": 28,
  "servingSizeUnit": "g",
  "householdServingFullText": "1 oz",
  "foodCategory": "Snacks",
  "modifiedDate": "2019-03-15",
  "foodNutrients": [
    {
      "id": 4544789,
      "nutrient": {
        "id": 1008,
        "number": "208",
        "name": "Energy",
        "rank": 300,
        "unitName": "kcal"
      },
      "type": "FoodNutrient",
      "amount": 500,
      "foodNutrientDerivation": {
        "code": "LCCS",
        "description": "Calculated from value per serving size measure"
      }
    }
  ],
  "labelNutrients": {
    "fat": {"value": 9.0},
    "saturatedFat": {"value": 1.0},
    "transFat": {"value": 0.0},
    "cholesterol": {"value": 0.0},
    "sodium": {"value": 5.0},
    "carbohydrates": {"value": 12.0},
    "fiber": {"value": 2.0},
    "sugars": {"value": 8.0},
    "protein": {"value": 4.0},
    "calcium": {"value": 20.0},
    "iron": {"value": 0.72},
    "calories": {"value": 140.0}
  },
  "foodUpdateLog": [
    {
      "fdcId": 534358,
      "availableDate": "2019-04-01",
      "brandOwner": "Kar Nut Products Company",
      "dataSource": "LI",
      "dataType": "Branded",
      "description": "NUT 'N BERRY MIX",
      "gtinUpc": "077034085228",
      "householdServingFullText": "1 oz",
      "ingredients": "PEANUTS, RAISINS...",
      "modifiedDate": "2019-03-15",
      "publicationDate": "2019-04-01",
      "servingSize": 28,
      "servingSizeUnit": "g"
    }
  ]
}
```

### SR Legacy Food Item
```json
{
  "fdcId": 170567,
  "description": "Nuts, almonds",
  "dataType": "SR Legacy",
  "publicationDate": "2019-04-01",
  "foodClass": "FinalFood",
  "scientificName": "Prunus dulcis",
  "ndbNumber": "12061",
  "nutrientConversionFactors": [
    {
      "type": "CalorieConversionFactor",
      "proteinValue": 3.47,
      "fatValue": 8.37,
      "carbohydrateValue": 4.07
    },
    {
      "type": "ProteinConversionFactor",
      "value": 5.18
    }
  ],
  "foodCategory": {
    "id": 12,
    "code": "1200",
    "description": "Nut and Seed Products"
  },
  "foodNutrients": [
    {
      "id": 1783958,
      "nutrient": {
        "id": 1003,
        "number": "203",
        "name": "Protein",
        "rank": 600,
        "unitName": "g"
      },
      "type": "FoodNutrient",
      "amount": 21.15,
      "dataPoints": 90
    }
  ],
  "foodPortions": [
    {
      "id": 89234,
      "amount": 1,
      "gramWeight": 1.2,
      "portionDescription": "almond",
      "sequenceNumber": 1
    },
    {
      "id": 89235,
      "amount": 1,
      "gramWeight": 143,
      "portionDescription": "cup whole kernels",
      "sequenceNumber": 6
    }
  ]
}
```

### Abridged Food Item (List Response)
```json
{
  "fdcId": 2262074,
  "description": "Almond butter, creamy",
  "dataType": "Foundation",
  "publicationDate": "2022-04-28",
  "ndbNumber": "12195",
  "foodNutrients": [
    {
      "number": "717",
      "name": "Daidzin",
      "amount": 0.0281,
      "unitName": "MG",
      "derivationCode": "A",
      "derivationDescription": "Analytical"
    },
    {
      "number": "208",
      "name": "Energy",
      "amount": 614,
      "unitName": "KCAL",
      "derivationCode": "NC",
      "derivationDescription": ""
    }
  ]
}
```

---

## Nutrient Fields Reference

### Core Nutrient IDs

| ID | Number | Name | Unit | Category |
|----|--------|------|------|----------|
| 1003 | 203 | Protein | g | Proximates |
| 1004 | 204 | Total lipid (fat) | g | Proximates |
| 1005 | 205 | Carbohydrate, by difference | g | Proximates |
| 1008 | 208 | Energy | kcal | Energy |
| 1051 | 255 | Water | g | Proximates |
| 1079 | 291 | Fiber, total dietary | g | Proximates |
| 1087 | 301 | Calcium, Ca | mg | Minerals |
| 1089 | 303 | Iron, Fe | mg | Minerals |
| 1090 | 304 | Magnesium, Mg | mg | Minerals |
| 1091 | 305 | Phosphorus, P | mg | Minerals |
| 1092 | 306 | Potassium, K | mg | Minerals |
| 1093 | 307 | Sodium, Na | mg | Minerals |
| 1095 | 309 | Zinc, Zn | mg | Minerals |
| 1098 | 312 | Copper, Cu | mg | Minerals |
| 1101 | 315 | Manganese, Mn | mg | Minerals |
| 1103 | 317 | Selenium, Se | mcg | Minerals |
| 1104 | 318 | Vitamin A, IU | IU | Vitamins |
| 1106 | 320 | Vitamin A, RAE | mcg | Vitamins |
| 1109 | 323 | Vitamin E | mg | Vitamins |
| 1162 | 401 | Vitamin C | mg | Vitamins |
| 1165 | 404 | Thiamin | mg | Vitamins |
| 1166 | 405 | Riboflavin | mg | Vitamins |
| 1167 | 406 | Niacin | mg | Vitamins |
| 1170 | 410 | Pantothenic acid | mg | Vitamins |
| 1175 | 415 | Vitamin B-6 | mg | Vitamins |
| 1177 | 417 | Folate, total | mcg | Vitamins |
| 1178 | 418 | Vitamin B-12 | mcg | Vitamins |
| 1180 | 421 | Choline, total | mg | Vitamins |
| 1185 | 430 | Vitamin K | mcg | Vitamins |

### Derivation Codes

| Code | Description |
|------|-------------|
| A | Analytical |
| AR | Analytical, averaged from multiple samples |
| AS | Analytical, averaged, composite sample |
| BFNN | Based on food and nutrient source |
| BFZN | Based on another food, zero nutrient |
| BP | Based on label, protein calculated |
| LCCD | Calculated from daily value percentage |
| LCCS | Calculated from value per serving size |
| NC | Calculated from nutrient profile |
| NR | Calculated, assumed zero |

---

## Download Files

### December 2025 Release

| Data Type | CSV | JSON |
|-----------|-----|------|
| Foundation Foods | 3.4 MB (29 MB uncompressed) | 467 KB (6.5 MB) |
| Branded Foods | 427 MB (2.9 GB) | 195 MB (3.1 GB) |
| SR Legacy | Static (2018) | Static (2018) |
| FNDDS (Survey) | Annual release | Annual release |
| Full Database | 458 MB (3.1 GB) | Combined |

### Download URL
```
https://fdc.nal.usda.gov/download-datasets/
```

---

## Food Counts by Data Type

| Data Type | Approximate Count | Update Frequency |
|-----------|-------------------|------------------|
| Foundation Foods | 205+ | Biannual (Apr, Dec) |
| SR Legacy | 267+ | Static (final 2018) |
| Survey (FNDDS) | 374+ | Annual |
| Branded Foods | 20,000+ | Biannual (Apr, Dec) |
| **Total** | ~20,900+ | Varies |

---

## Response Status Codes

| Code | Description |
|------|-------------|
| 200 | Success |
| 400 | Bad input parameter |
| 404 | No results found |

---

## Example API Calls

### Search for Foods
```bash
curl "https://api.nal.usda.gov/fdc/v1/foods/search?api_key=DEMO_KEY&query=Cheddar%20Cheese"
```

### Get Single Food by ID
```bash
curl "https://api.nal.usda.gov/fdc/v1/food/2346407?api_key=DEMO_KEY"
```

### List Foundation Foods
```bash
curl "https://api.nal.usda.gov/fdc/v1/foods/list?api_key=DEMO_KEY&dataType=Foundation&pageSize=25"
```

### POST Search with Filters
```bash
curl -X POST "https://api.nal.usda.gov/fdc/v1/foods/search?api_key=DEMO_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "query": "chicken",
    "dataType": ["Foundation", "SR Legacy"],
    "pageSize": 10,
    "pageNumber": 1,
    "sortBy": "dataType.keyword",
    "sortOrder": "asc"
  }'
```

---

## License

- **Data License:** CC0 1.0 Universal (Public Domain)
- **Suggested Citation:** U.S. Department of Agriculture, Agricultural Research Service. FoodData Central, [Year]. fdc.nal.usda.gov.

---

## Related Documents

- [../diseases/SLEEP-LONGEVITY-NUTRI.md](../../diseases/SLEEP-LONGEVITY-NUTRI.md) - References USDA FoodData
- [../INDEX.md](../../INDEX.md) - Data sources index
