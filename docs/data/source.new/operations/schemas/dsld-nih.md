---
id: schema-dsld-nih
title: "DSLD (Dietary Supplement Label Database) - Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./_index.md)

# DSLD (Dietary Supplement Label Database) - Schema Documentation

**Document ID:** SCHEMA-DSLD-NIH
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://dsld.od.nih.gov/
**API Base URL:** https://api.ods.od.nih.gov/dsld/v9/

---

## Overview

The Dietary Supplement Label Database (DSLD) is maintained by the NIH Office of Dietary Supplements (ODS). It contains full label information from dietary supplement products marketed in the United States, including both currently marketed and historical (off-market) products.

---

## Database Statistics

| Metric | Value |
|--------|-------|
| Total Labels | 200,000+ |
| As of 2021 | 125,565 product labels |
| Initial Launch (2013) | 16,712 labels |
| Coverage | On-market and off-market products |
| Update Frequency | Continuous |

---

## API Specification

### Base URL
```
https://api.ods.od.nih.gov/dsld/v9/
```

### Authentication
- **Type:** None required (public API)
- **Rate Limits:** Not documented; recommend using `size=10` for testing

---

## Endpoints

### 1. Browse Products
```
GET /v9/browse-products
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| letter | string | No | Filter by starting letter (A-Z) |
| keyword | string | No | Search within product names |
| size | integer | No | Results per page (default: 25) |
| from | integer | No | Offset for pagination |

### 2. Browse Brands
```
GET /v9/browse-brands
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| letter | string | No | Filter by starting letter |
| keyword | string | No | Search within brand names |
| size | integer | No | Results per page |
| from | integer | No | Offset for pagination |

### 3. Brand Products
```
GET /v9/brand-products
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| brand | string | Yes | Brand name to filter by |
| size | integer | No | Results per page |
| from | integer | No | Offset for pagination |

### 4. Get Label by ID
```
GET /v9/label/{id}
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| id | integer | Yes | DSLD label ID |

### 5. Search Filter (Advanced)
```
GET /v9/search-filter
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| q | string | No | Search query |
| product_type | string | No | Filter by product type code |
| ingredient_group | string | No | Filter by ingredient category |
| size | integer | No | Results per page |
| from | integer | No | Offset for pagination |

### 6. Search Filter Histogram
```
GET /v9/search-filter-histogram
```

Returns time-based distribution of search results for trend analysis.

### 7. Ingredient Groups
```
GET /v9/ingredient-groups
```

Returns list of all ingredient categories with synonyms.

### 8. Version
```
GET /v9/version
```

Returns API version information.

---

## Ingredient Categories

The DSLD classifies ingredients into 18 primary categories:

| Category Code | Category Name |
|---------------|---------------|
| amino_acid | Amino Acid |
| animal_part_or_source | Animal Part or Source |
| blend | Blend |
| botanical | Botanical |
| sugars | Sugars |
| complex_carbohydrate | Complex Carbohydrate |
| enzyme | Enzyme |
| fat | Fat |
| fatty_acid | Fatty Acid |
| fiber | Fiber |
| hormone | Hormone |
| chemical | Chemical |
| bacteria | Bacteria (Probiotic) |
| protein | Protein |
| vitamin | Vitamin |
| mineral | Mineral |
| other | Other |
| tbd | To Be Determined |

---

## JSON Schemas

### Label Object (Full)
```json
{
  "id": 3000,
  "fullName": "American Ginseng",
  "brandName": "Root To Health",
  "upcSku": "0 53181 01000 9",
  "manufacturer": {
    "name": "Hsu's Ginseng Enterprises, Inc.",
    "address": {
      "street": "Street Address",
      "city": "Wausau",
      "state": "WI",
      "country": "USA"
    },
    "phone": "1-800-826-1577"
  },
  "netContents": {
    "value": 50,
    "unit": "capsules",
    "weight": {
      "value": 25,
      "unit": "g"
    }
  },
  "servingsPerContainer": 50,
  "servingSize": {
    "value": 1,
    "unit": "capsule"
  },
  "physicalState": "Capsule",
  "targetGroups": [
    {
      "group": "Adults",
      "ageRange": "18-50 years"
    }
  ],
  "ingredientRows": [
    {
      "name": "American Ginseng",
      "quantity": 500,
      "unit": "mg",
      "dailyValue": null,
      "dailyValueUnit": null,
      "dailyValueTargetGroup": "Adults",
      "uniiCode": null,
      "ingredientGroup": "botanical",
      "synonyms": ["Panax quinquefolius"]
    }
  ],
  "otherIngredients": [
    {
      "name": "Gelatin",
      "purpose": "capsule material"
    },
    {
      "name": "Water",
      "purpose": null
    }
  ],
  "statements": [
    {
      "type": "disclaimer",
      "text": "*These statements have not been evaluated by the Food and Drug Administration. This product is not intended to diagnose, treat, cure, or prevent disease."
    }
  ],
  "claims": [
    {
      "type": "structure_function",
      "text": "Stimulates Energy and Helps Reduce Stress"
    }
  ],
  "entryDate": "2011-12-23",
  "marketStatus": {
    "status": "off_market",
    "statusDate": "2023-11-17"
  }
}
```

### Label Object (GNC Coral Calcium Example)
```json
{
  "id": 1000,
  "fullName": "Coral Calcium",
  "brandName": "GNC",
  "upcSku": "0 48107 05843 2",
  "manufacturer": {
    "name": "General Nutrition Corporation",
    "address": {
      "city": "Pittsburgh",
      "state": "PA",
      "zip": "15222"
    },
    "website": "www.gnc.com"
  },
  "netContents": {
    "value": 60,
    "unit": "capsules"
  },
  "servingsPerContainer": 30,
  "servingSize": {
    "value": 2,
    "unit": "capsules"
  },
  "physicalState": "Capsule",
  "targetGroups": [
    {
      "group": "Adults",
      "ageRange": "18-50 years"
    }
  ],
  "ingredientRows": [
    {
      "name": "Vitamin D",
      "form": "Cholecalciferol",
      "quantity": 200,
      "unit": "IU",
      "dailyValue": 50,
      "dailyValueUnit": "%",
      "ingredientGroup": "vitamin"
    },
    {
      "name": "Calcium",
      "form": "Coral Calcium",
      "quantity": 400,
      "unit": "mg",
      "dailyValue": 40,
      "dailyValueUnit": "%",
      "ingredientGroup": "mineral"
    },
    {
      "name": "Magnesium",
      "form": "Magnesium Oxide",
      "quantity": 200,
      "unit": "mg",
      "dailyValue": 50,
      "dailyValueUnit": "%",
      "ingredientGroup": "mineral"
    }
  ],
  "otherIngredients": [
    {"name": "Gelatin"},
    {"name": "Titanium Dioxide"}
  ],
  "marketStatus": {
    "status": "off_market",
    "statusDate": "2022-12-13"
  },
  "productOrigin": "Made in the USA; sourced from fossilized coral above the Okinawan sea"
}
```

### Search Response
```json
{
  "hits": [
    {
      "id": 12345,
      "fullName": "Product Name",
      "brandName": "Brand",
      "ingredientRows": [...],
      "marketStatus": {...}
    }
  ],
  "stats": {
    "count": 150,
    "pct": 100
  },
  "aggregations": {
    "product_types": [...],
    "ingredient_groups": [...]
  }
}
```

### Histogram Response
```json
[
  {
    "key_as_string": "2020-01-01",
    "key": 1577836800000,
    "doc_count": 245
  },
  {
    "key_as_string": "2020-02-01",
    "key": 1580515200000,
    "doc_count": 312
  }
]
```

### Ingredient Row Object
```json
{
  "name": "Vitamin D",
  "form": "Cholecalciferol",
  "quantity": 200,
  "unit": "IU",
  "dailyValue": 50,
  "dailyValueUnit": "%",
  "dailyValueTargetGroup": "Adults",
  "uniiCode": "1C6V77QF41",
  "ingredientGroup": "vitamin",
  "synonyms": ["Vitamin D3", "Cholecalciferol"],
  "partUsed": null,
  "extractRatio": null
}
```

---

## UNII Codes

UNII (Unique Ingredient Identifier) codes are standardized identifiers assigned by the FDA's Substance Registration System.

### Format
- 10-character alphanumeric string
- Example: `19F5HK2737` (Vitamin C/Ascorbic Acid)
- Example: `1C6V77QF41` (Vitamin D3/Cholecalciferol)

### Usage in DSLD
- Stored in `ingredientRows[].uniiCode` field
- Not all ingredients have UNII codes assigned
- Botanical ingredients often lack UNII codes
- Used for cross-referencing with FDA databases

### UNII Lookup
- FDA SRS: https://precision.fda.gov/uniisearch
- DailyMed: https://dailymed.nlm.nih.gov/dailymed/

---

## Data Fields Reference

### Product Identification

| Field | Type | Description |
|-------|------|-------------|
| id | integer | DSLD internal ID |
| fullName | string | Full product name |
| brandName | string | Brand/manufacturer name |
| upcSku | string | UPC or SKU code |

### Physical Characteristics

| Field | Type | Description |
|-------|------|-------------|
| physicalState | string | Form (Capsule, Tablet, Liquid, Powder, etc.) |
| servingsPerContainer | number | Number of servings |
| servingSize | object | Amount and unit per serving |
| netContents | object | Package contents (count, weight) |

### Market Information

| Field | Type | Description |
|-------|------|-------------|
| entryDate | string (date) | Date added to database |
| marketStatus.status | string | `on_market` or `off_market` |
| marketStatus.statusDate | string (date) | Date of status change |

### Regulatory Fields

| Field | Type | Description |
|-------|------|-------------|
| statements | array | FDA disclaimers, warnings |
| claims | array | Structure/function claims |
| events | array | Adverse event reports (if any) |

---

## Physical States (Product Forms)

| Physical State | Description |
|----------------|-------------|
| Capsule | Hard or soft gel capsule |
| Tablet | Compressed tablet |
| Softgel | Soft gelatin capsule |
| Liquid | Liquid solution or suspension |
| Powder | Loose or encapsulated powder |
| Gummy | Gummy/chewable form |
| Lozenge | Dissolving lozenge |
| Spray | Oral spray |
| Drops | Liquid drops |
| Bar | Nutrition/protein bar |
| Wafer | Chewable wafer |

---

## Example API Calls

### Browse Products Starting with "V"
```bash
curl "https://api.ods.od.nih.gov/dsld/v9/browse-products?letter=V&size=10"
```

### Get Specific Label
```bash
curl "https://api.ods.od.nih.gov/dsld/v9/label/3000"
```

### Search for Vitamin D Products
```bash
curl "https://api.ods.od.nih.gov/dsld/v9/search-filter?q=vitamin%20D&size=10"
```

### Get Products by Brand
```bash
curl "https://api.ods.od.nih.gov/dsld/v9/brand-products?brand=GNC&size=10"
```

### Filter by Ingredient Group
```bash
curl "https://api.ods.od.nih.gov/dsld/v9/search-filter?ingredient_group=botanical&size=10"
```

### Get Historical Trend
```bash
curl "https://api.ods.od.nih.gov/dsld/v9/search-filter-histogram?q=melatonin"
```

---

## Data Quality Notes

### Completeness
- Not all labels have complete data
- UNII codes may be missing for some ingredients
- Daily values not established for all ingredients
- Some historical products lack detailed information

### Updates
- Database updated regularly with new products
- Formulation changes tracked
- Off-market status updated when products discontinued

### Limitations
- API documentation is sparse
- Some endpoints return 500 errors under load
- Empty result sets return minimal structure

---

## Timeline

| Date | Milestone |
|------|-----------|
| June 2013 | Initial launch (16,712 labels) |
| 2017 | Mobile redesign |
| 2020 | Platform modernization |
| 2021 | 125,565 labels |
| 2026 | 200,000+ labels |

---

## License

- **Access:** Public, free access
- **Usage:** Research and educational purposes
- **Attribution:** NIH Office of Dietary Supplements

---

## Citation

```
Office of Dietary Supplements. Dietary Supplement Label Database (DSLD).
National Institutes of Health.
https://dsld.od.nih.gov/
```

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | JSON (REST API) |
| Alternative | None (web only) |
| Compression | None |
| Encoding | UTF-8 |
| API Response | JSON |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 30,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| DSLD | HTTP | https://dsld.nlm.nih.gov/dsld.html |
| DSLD Data | FTP | See main database |

**Access Requirements:** Open access via NIH, no registration required

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| DSLD | Public Domain (NIH) | Yes |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "3000" |
| `name` | string | Entity name | "American Ginseng" |
| `type` | string | Record type | "product" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Sample Data

### Example Record
```json
{
  "id": 3000,
  "fullName": "American Ginseng",
  "brandName": "Root To Health",
  "ingredientGroup": "botanical"
}
```

### Sample Query Result
| id | fullName | brandName | ingredientGroup |
|----|----------|-----------|-----------------|
| 3000 | American Ginseng | Root To Health | botanical |
| 1000 | Coral Calcium | GNC | mineral |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `id` | DSLD internal identifier for a product label | `3000` |
| `fullName` | Complete product name as shown on label | `American Ginseng` |
| `brandName` | Manufacturer or brand name | `GNC` |
| `upcSku` | Universal Product Code or Stock Keeping Unit | `0 48107 05843 2` |
| `physicalState` | Form of the supplement product | `Capsule` |
| `servingSize` | Amount constituting one serving with unit | `2 capsules` |
| `dailyValue` | Percentage of recommended daily intake | `50` |
| `uniiCode` | FDA Unique Ingredient Identifier | `1C6V77QF41` |
| `ingredientGroup` | Category classification for ingredients | `botanical` |
| `marketStatus` | Whether product is currently sold or discontinued | `on_market` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Daily Value (DV) | FDA reference amount for nutrients based on 2000-calorie diet | dailyValue |
| Supplement Facts | FDA-required nutrition label panel for dietary supplements | ingredientRows |
| Structure/Function Claim | Statement about supplement's effect on body structure or function | claims |
| Other Ingredients | Non-active ingredients in supplement formulation | otherIngredients |
| Net Contents | Total quantity of product in package by count or weight | servingsPerContainer |
| Target Group | Population for which daily values are calculated | `Adults` |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DSLD | Dietary Supplement Label Database | NIH product label database |
| NIH | National Institutes of Health | US health research agency |
| ODS | Office of Dietary Supplements | NIH division maintaining DSLD |
| UNII | Unique Ingredient Identifier | FDA substance registration |
| FDA | Food and Drug Administration | US regulatory agency |
| UPC | Universal Product Code | Barcode standard |
| SKU | Stock Keeping Unit | Inventory identifier |
| DV | Daily Value | Nutrient reference amount |
| IU | International Unit | Vitamin potency measure |
| SRS | Substance Registration System | FDA UNII database |

---

## Related Documents

- [sleep-longevity-nutri.md](../diseases/sleep-longevity-nutri.md) - Nutrigenomics context
- [pharmaceuticals.md](../compounds/pharmaceuticals.md) - Related drug databases
- [usda-fooddata-central.md](./usda-fooddata-central.md) - Food composition data
- [foodb.md](./foodb.md) - Food compound data
