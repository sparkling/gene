---
id: schema-dsld
title: "DSLD - Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, supplements, nih, labels, ingredients]
---

# DSLD - Dietary Supplement Label Database Schema Documentation

**Document ID:** SCHEMA-DSLD
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://dsld.od.nih.gov/

---

## TL;DR

The Dietary Supplement Label Database (DSLD) is maintained by the NIH Office of Dietary Supplements. It contains full label information from dietary supplements marketed in the United States, including supplement facts, ingredient lists, suggested use, and manufacturer information.

---

## Database Statistics

| Metric | Count |
|--------|-------|
| Products | 150,000+ |
| Unique Ingredients | 13,000+ |
| Brands | 4,000+ |
| Manufacturers | 2,500+ |
| Product Categories | 15+ |

---

## Database Schema

### Core Tables

#### 1. Products Table
```sql
CREATE TABLE products (
    dsld_id INTEGER PRIMARY KEY,      -- DSLD internal ID
    product_name VARCHAR(500),        -- Full product name
    brand_name VARCHAR(255),          -- Brand name
    net_contents VARCHAR(100),        -- Net quantity/count
    serving_size VARCHAR(100),        -- Serving size description
    servings_per_container VARCHAR(50), -- Number of servings
    product_type VARCHAR(100),        -- Form (capsule, tablet, etc.)
    dietary_claims TEXT,              -- Health/dietary claims
    upc VARCHAR(20),                  -- UPC barcode
    nhp_id VARCHAR(20),               -- Canadian NHP ID (if applicable)
    date_entered DATE,                -- Database entry date
    last_updated DATE,                -- Last modification date
    INDEX idx_brand (brand_name),
    INDEX idx_upc (upc)
);
```

#### 2. Ingredients Table
```sql
CREATE TABLE ingredients (
    ingredient_id INTEGER PRIMARY KEY,
    dsld_id INTEGER REFERENCES products(dsld_id),
    ingredient_name VARCHAR(500),     -- Ingredient name as on label
    amount VARCHAR(50),               -- Amount per serving
    unit VARCHAR(20),                 -- Unit (mg, mcg, IU, etc.)
    daily_value_percent VARCHAR(10),  -- % Daily Value
    ingredient_group VARCHAR(100),    -- Category/group
    source_name VARCHAR(255),         -- Source organism/material
    is_blend_component BOOLEAN,       -- Part of proprietary blend
    blend_name VARCHAR(255),          -- Proprietary blend name
    INDEX idx_ingredient (ingredient_name),
    INDEX idx_dsld (dsld_id)
);
```

#### 3. Manufacturers Table
```sql
CREATE TABLE manufacturers (
    manufacturer_id INTEGER PRIMARY KEY,
    company_name VARCHAR(255),        -- Company name
    address VARCHAR(500),             -- Full address
    city VARCHAR(100),
    state VARCHAR(50),
    country VARCHAR(100),
    contact_info VARCHAR(255),        -- Phone/website
    INDEX idx_company (company_name)
);
```

#### 4. Product-Manufacturer Junction
```sql
CREATE TABLE product_manufacturers (
    dsld_id INTEGER REFERENCES products(dsld_id),
    manufacturer_id INTEGER REFERENCES manufacturers(manufacturer_id),
    relationship_type VARCHAR(50),    -- "Manufacturer", "Distributor", etc.
    PRIMARY KEY (dsld_id, manufacturer_id)
);
```

---

## API Specification

### Base URL
```
https://api.ods.od.nih.gov/dsld/v8/
```

### Authentication
- **Type:** API Key (header)
- **Header:** `X-Api-Key`
- **Registration:** https://api.ods.od.nih.gov/signup

### Endpoints

#### 1. Search Products
```
GET /browse
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| q | string | No | Search query |
| brandName | string | No | Filter by brand |
| ingredient | string | No | Filter by ingredient |
| pagesize | integer | No | Results per page (default: 25) |
| offset | integer | No | Pagination offset |

#### 2. Get Product Details
```
GET /label/{dsld_id}
```

Returns complete label information for a specific product.

#### 3. Get Ingredient List
```
GET /ingredients/{dsld_id}
```

Returns all ingredients for a specific product.

#### 4. Search Ingredients
```
GET /ingredient/browse
```

| Parameter | Type | Description |
|-----------|------|-------------|
| name | string | Ingredient name search |
| category | string | Ingredient category |

---

## JSON Schemas

### Product Object
```json
{
  "dsldId": 123456,
  "productName": "Vitamin D3 5000 IU",
  "brandName": "Nature's Way",
  "netContents": "120 Softgels",
  "servingSize": "1 softgel",
  "servingsPerContainer": "120",
  "productType": "Softgel",
  "upc": "033674100233",
  "dateEntered": "2023-05-15",
  "manufacturer": {
    "name": "Nature's Way Products",
    "address": "825 Challenger Drive",
    "city": "Green Bay",
    "state": "WI",
    "country": "USA"
  },
  "dietaryClaims": [
    "Gluten Free",
    "No artificial colors"
  ]
}
```

### Ingredient Object
```json
{
  "ingredientId": 789012,
  "dsldId": 123456,
  "ingredientName": "Vitamin D3 (as Cholecalciferol)",
  "amount": "5000",
  "unit": "IU",
  "dailyValuePercent": "1250",
  "ingredientGroup": "Vitamins",
  "sourceName": "Lanolin",
  "isBlendComponent": false,
  "blendName": null
}
```

### Search Results Object
```json
{
  "totalCount": 15234,
  "pageSize": 25,
  "offset": 0,
  "results": [
    {
      "dsldId": 123456,
      "productName": "...",
      "brandName": "...",
      "...": "..."
    }
  ]
}
```

---

## Data Categories

### Product Types

| Type | Description |
|------|-------------|
| Capsule | Hard or soft capsule |
| Tablet | Compressed tablet |
| Softgel | Soft gelatin capsule |
| Liquid | Liquid form |
| Powder | Powder form |
| Gummy | Gummy/chewable |
| Lozenge | Dissolving lozenge |
| Spray | Oral spray |

### Ingredient Groups

| Group | Examples |
|-------|----------|
| Vitamins | Vitamin A, B12, C, D, E, K |
| Minerals | Calcium, Iron, Magnesium, Zinc |
| Herbs/Botanicals | Echinacea, Ginkgo, Turmeric |
| Amino Acids | L-Arginine, L-Glutamine, BCAA |
| Enzymes | Lipase, Protease, Amylase |
| Probiotics | Lactobacillus, Bifidobacterium |
| Omega Fatty Acids | Fish oil, Flaxseed oil |
| Specialty | CoQ10, Glucosamine, Melatonin |

---

## Data Access

### API Access

```bash
# Get API key from: https://api.ods.od.nih.gov/signup

# Search for vitamin D products
curl -H "X-Api-Key: YOUR_API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/browse?q=vitamin%20d&pagesize=10"

# Get specific product
curl -H "X-Api-Key: YOUR_API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/label/123456"
```

### Bulk Downloads

| Format | URL | Description |
|--------|-----|-------------|
| Full Database | https://dsld.od.nih.gov/download | Complete export |
| API Documentation | https://api.ods.od.nih.gov/docs | Swagger/OpenAPI |

---

## Use Cases

### 1. Find Products by Ingredient
```bash
curl -H "X-Api-Key: $API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/browse?ingredient=ashwagandha"
```

### 2. Get Label Information
```bash
curl -H "X-Api-Key: $API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/label/123456" | jq
```

### 3. Search by Brand
```bash
curl -H "X-Api-Key: $API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/browse?brandName=Nature%27s%20Way"
```

---

## Relationships

### Entity Relationship Diagram
```
products (1) ----< ingredients (N)
    |
    +---- (N) manufacturers (M)
    |
    +---- (1) labels (1)
```

---

## License

- **License Type:** Public Domain (US Government)
- **Commercial Use:** Yes
- **Attribution:** Recommended
- **API Terms:** Requires registration for API key

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| DSLD ID | Unique product identifier in database | 123456 |
| Daily Value | FDA recommended daily intake percentage | 1250% |
| Net Contents | Total amount in package | 120 Softgels |
| Serving Size | Amount per single serving | 1 softgel |
| Proprietary Blend | Combined ingredients without individual amounts | "Energy Blend" |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DSLD | Dietary Supplement Label Database | NIH database |
| ODS | Office of Dietary Supplements | NIH office |
| NIH | National Institutes of Health | US agency |
| NHP | Natural Health Product | Canadian identifier |
| UPC | Universal Product Code | Barcode |
| IU | International Unit | Vitamin measurement |
| DV | Daily Value | Nutrient reference |

---

## Related Documents

- [Download Instructions](./download.md)
- [ConsumerLab](../consumerlab/_index.md) - Product testing
- [Natural Medicines](../natural.medicines/_index.md) - Evidence database
