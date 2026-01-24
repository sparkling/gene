---
id: schemas-open-food-facts
title: Open Food Facts Schema Documentation
category: schemas
subcategory: nutrition
tier: 2
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, nutrition, food, crowdsourced, open-data]
---

**Parent:** [Schema Documentation](./_index.md)

# Open Food Facts - Schema Documentation

**Document ID:** SCHEMA-OPEN-FOOD-FACTS
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://world.openfoodfacts.org/
**API Documentation:** https://wiki.openfoodfacts.org/API

---

## TL;DR

Open Food Facts is the world's largest open food products database with 3M+ products, providing detailed nutrition facts, ingredients, allergens, and environmental scores. Data is crowdsourced and available under ODbL license. API supports JSON responses with comprehensive product data including Nutri-Score, Eco-Score, and NOVA classification.

---

## Database Statistics

| Metric | Value |
|--------|-------|
| Total Products | 3,000,000+ |
| Countries | 180+ |
| Contributors | 30,000+ |
| Languages | 50+ |
| Daily Additions | ~10,000 products |
| Data Formats | JSON, JSONL, Parquet, CSV, RDF |

---

## License

| Attribute | Value |
|-----------|-------|
| Database License | ODbL (Open Database License) 1.0 |
| Content License | Database Contents License (DbCL) 1.0 |
| Images License | CC BY-SA 3.0 |
| Commercial Use | Permitted with attribution |

---

## Data Access

### REST API (v2)

**Base URL:** `https://world.openfoodfacts.org/api/v2/`

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/product/{barcode}` | GET | Get single product by barcode |
| `/search` | GET | Search products with filters |

### Bulk Downloads

| Format | URL | Use Case |
|--------|-----|----------|
| CSV | https://static.openfoodfacts.org/data/en.openfoodfacts.org.products.csv.gz | Spreadsheet analysis |
| JSONL | https://static.openfoodfacts.org/data/openfoodfacts-products.jsonl.gz | Streaming processing |
| Parquet | https://huggingface.co/datasets/openfoodfacts/product-database | ML/analytics |

---

## Product Schema

### Core Identification Fields

| Field | Type | Description |
|-------|------|-------------|
| code | string | Product barcode (EAN-13 or internal) |
| product_name | string | Product name |
| brands | string | Brand name(s), comma-separated |
| categories_tags | array | Hierarchical category tags |
| quantity | string | Package quantity and unit |

### Nutrition Facts Fields

| Nutrient Key | Unit | Description |
|--------------|------|-------------|
| energy-kcal | kcal | Energy in calories |
| fat | g | Total fat |
| saturated-fat | g | Saturated fatty acids |
| carbohydrates | g | Total carbohydrates |
| sugars | g | Total sugars |
| fiber | g | Dietary fiber |
| proteins | g | Protein |
| salt | g | Salt |

### Nutrition Scores

| Field | Type | Description |
|-------|------|-------------|
| nutriscore_grade | string | Nutri-Score grade (a-e) |
| nova_group | integer | NOVA group (1-4) |
| ecoscore_grade | string | Eco-Score grade (a-e) |

---

## NOVA Classification

| Group | Description |
|-------|-------------|
| 1 | Unprocessed or minimally processed foods |
| 2 | Processed culinary ingredients |
| 3 | Processed foods |
| 4 | Ultra-processed foods |

---

## Download

### Bulk Data Downloads

| Source | Format | Size | URL |
|--------|--------|------|-----|
| CSV Export | CSV | ~800 MB (compressed) | https://static.openfoodfacts.org/data/en.openfoodfacts.org.products.csv.gz |
| JSONL Export | JSONL | ~1.2 GB (compressed) | https://static.openfoodfacts.org/data/openfoodfacts-products.jsonl.gz |
| Parquet (Hugging Face) | Parquet | ~600 MB | https://huggingface.co/datasets/openfoodfacts/product-database |
| RDF Dump | RDF/Turtle | ~5 GB | https://static.openfoodfacts.org/data/openfoodfacts-products.rdf.gz |

### Update Frequency

- New products: Added continuously (average 10,000+ per day)
- Bulk exports: Updated daily
- Database snapshots: Available monthly

### API Rate Limits

- No authentication required
- Rate limit: ~1 request per second recommended
- Bulk operations: Use download files instead

---

## Data Format

| Aspect | Details |
|--------|---------|
| **Primary Format** | JSON (API), CSV (bulk exports) |
| **Alternative Formats** | JSONL, RDF, Parquet |
| **Compression** | gzip (.gz) for bulk downloads |
| **Encoding** | UTF-8 |
| **API Response Format** | JSON |
| **Bulk Export Encoding** | UTF-8 with BOM for Excel compatibility |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| **Total Products** | 3,000,000+ (updated continuously) |
| **Countries Represented** | 180+ |
| **Active Contributors** | 30,000+ |
| **Languages** | 50+ |
| **Daily New Entries** | ~10,000 |
| **CSV Export Size** | ~800 MB (gzip compressed) |
| **JSONL Export Size** | ~1.2 GB (gzip compressed) |
| **Parquet Dataset** | ~600 MB |
| **Update Frequency** | Continuous (real-time API), daily bulk export |
| **Last Full Snapshot** | Updated daily at 02:00 UTC |

---

## Sample Data

### Example Record
```json
{
  "barcode": "3017620425035",
  "product_name": "Nutella hazelnut spread",
  "brand": "Ferrero",
  "categories": "Spreads, Cocoa and hazelnuts",
  "energy_kcal": 539,
  "protein_g": 6.3,
  "fat_g": 30.7,
  "carbs_g": 57.5
}
```

### Sample Query Result
| barcode | product_name | brand | energy_kcal | protein_g | fat_g |
|---------|-------------|-------|------------|-----------|--------|
| 3017620425035 | Nutella hazelnut spread | Ferrero | 539 | 6.3 | 30.7 |
| 8711200762022 | Coca-Cola Original | Coca-Cola | 42 | 0.0 | 0.0 |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "5000112118957" |
| `name` | string | Entity name | "Product Name" |
| `type` | string | Record type | "food_product" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `contains` | Ingredient | N:M |
| `has_label` | Label / Allergen | N:M |
| `belongs_to` | Category | N:M |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `code` | Product barcode identifier (EAN-13 or internal) | `3017620422003` |
| `product_name` | Name of the food product | `Nutella` |
| `brands` | Brand name(s), comma-separated if multiple | `Ferrero` |
| `categories_tags` | Hierarchical category classification tags | `en:breakfast-cereals` |
| `nutriscore_grade` | Nutri-Score nutritional grade from A (best) to E (worst) | `c` |
| `nova_group` | NOVA food processing classification (1-4) | `4` |
| `ecoscore_grade` | Environmental impact grade from A (best) to E (worst) | `b` |
| `energy-kcal` | Energy content in kilocalories per 100g | `539` |
| `quantity` | Package quantity with unit | `400 g` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Nutri-Score | Front-of-pack nutrition label using colors A-E | nutriscore_grade |
| NOVA Classification | System categorizing foods by processing level | nova_group |
| Eco-Score | Environmental impact score based on life cycle analysis | ecoscore_grade |
| Ultra-processed | NOVA group 4 foods with industrial additives | nova_group |
| Crowdsourced Data | User-contributed product information | Data collection method |
| Per 100g | Standard reference amount for nutrition comparison | Nutrient values |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| OFF | Open Food Facts | Crowdsourced food database |
| EAN | European Article Number | Barcode standard (EAN-13) |
| ODbL | Open Database License | Database license type |
| DbCL | Database Contents License | Content license type |
| NOVA | - | Food processing classification (not an acronym) |
| JSONL | JSON Lines | Newline-delimited JSON format |
| RDF | Resource Description Framework | Linked data format |
| ML | Machine Learning | Analysis/prediction applications |
| MUFA | Monounsaturated Fatty Acids | Healthy fat type |
| PUFA | Polyunsaturated Fatty Acids | Essential fatty acids |
| SFA | Saturated Fatty Acids | Fat type tracked in nutrition |

---

## Related Documents

- [foodb.md](./foodb.md) - Food composition database
- [usda-fooddata-central.md](./usda-fooddata-central.md) - USDA nutrition database
