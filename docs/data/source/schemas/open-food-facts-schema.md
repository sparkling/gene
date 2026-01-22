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

## Overview

Open Food Facts is a free, open, collaborative database of food products from around the world. The database includes:

- Product identification (barcodes, names, brands)
- Ingredient lists with analysis
- Nutrition facts per 100g and per serving
- Allergen and trace information
- Environmental impact scores (Eco-Score)
- Nutrition quality grades (Nutri-Score)
- Food processing levels (NOVA classification)
- Packaging information

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
| Share-Alike | Required for derivative databases |
| Attribution | "Source: Open Food Facts" required |

**Attribution Example:**
```
Data source: Open Food Facts (https://openfoodfacts.org)
License: Open Database License (ODbL)
```

---

## Data Access

### REST API (v2)

**Base URL:** `https://world.openfoodfacts.org/api/v2/`

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/product/{barcode}` | GET | Get single product by barcode |
| `/search` | GET | Search products with filters |
| `/cgi/product_jqm2.pl` | POST | Add/edit product |

**Example Requests:**
```bash
# Get product by barcode
curl "https://world.openfoodfacts.org/api/v2/product/737628064502.json"

# Search products
curl "https://world.openfoodfacts.org/api/v2/search?categories_tags=breakfast-cereals&page_size=20"

# Search with nutrition filter
curl "https://world.openfoodfacts.org/api/v2/search?nutrition_grades=a&page_size=50"
```

### Bulk Downloads

| Format | URL | Use Case |
|--------|-----|----------|
| CSV | https://static.openfoodfacts.org/data/en.openfoodfacts.org.products.csv.gz | Spreadsheet analysis |
| JSONL | https://static.openfoodfacts.org/data/openfoodfacts-products.jsonl.gz | Streaming processing |
| Parquet | https://huggingface.co/datasets/openfoodfacts/product-database | ML/analytics |
| MongoDB | https://static.openfoodfacts.org/data/openfoodfacts-mongodbdump.gz | Full database |
| Delta | Daily incremental updates | Sync updates |

### Country-Specific Instances

| Country | URL |
|---------|-----|
| World | https://world.openfoodfacts.org |
| USA | https://us.openfoodfacts.org |
| France | https://fr.openfoodfacts.org |
| UK | https://uk.openfoodfacts.org |
| Germany | https://de.openfoodfacts.org |

---

## Product Schema

### Core Identification Fields

| Field | Type | Description |
|-------|------|-------------|
| code | string | Product barcode (EAN-13 or internal) |
| product_name | string | Product name |
| product_name_{lc} | string | Localized product name |
| generic_name | string | Generic product category |
| brands | string | Brand name(s), comma-separated |
| brands_tags | array | Normalized brand tags |
| categories | string | Product categories |
| categories_tags | array | Hierarchical category tags |
| quantity | string | Package quantity and unit |

**Sample JSON:**
```json
{
  "code": "737628064502",
  "product_name": "Banana Chips",
  "product_name_en": "Banana Chips",
  "generic_name": "Dried banana slices",
  "brands": "Trader Joe's",
  "brands_tags": ["trader-joe-s"],
  "categories": "Snacks, Dried fruits, Banana chips",
  "categories_tags": ["en:snacks", "en:dried-fruits", "en:banana-chips"],
  "quantity": "170 g"
}
```

### Ingredients Fields

| Field | Type | Description |
|-------|------|-------------|
| ingredients_text | string | Raw ingredients text |
| ingredients_text_{lc} | string | Localized ingredients |
| ingredients | array | Parsed ingredient objects |
| ingredients_n | integer | Number of ingredients |
| additives_n | integer | Number of additives |
| additives_tags | array | Additive identifiers |
| allergens | string | Allergen information |
| allergens_tags | array | Normalized allergen tags |
| traces | string | May contain traces |
| traces_tags | array | Normalized trace tags |

**Ingredients Array Structure:**
```json
{
  "ingredients": [
    {
      "id": "en:banana",
      "text": "Bananas",
      "percent": 65,
      "percent_estimate": 65.5,
      "percent_min": 60,
      "percent_max": 70,
      "vegan": "yes",
      "vegetarian": "yes"
    },
    {
      "id": "en:coconut-oil",
      "text": "Coconut oil",
      "percent_estimate": 20,
      "from_palm_oil": "no"
    },
    {
      "id": "en:sugar",
      "text": "Sugar",
      "percent_estimate": 10
    }
  ]
}
```

### Nutrition Facts Fields

| Field | Type | Description |
|-------|------|-------------|
| nutriments | object | All nutrition data |
| nutrition_data | string | Nutrition source (product, calculated) |
| nutrition_data_per | string | "100g" or "serving" |
| serving_size | string | Serving size description |
| serving_quantity | number | Serving size in grams |

**Nutriments Object Structure:**

Each nutrient has multiple field variants:
- `{nutrient}` - Value as entered
- `{nutrient}_100g` - Per 100g/100ml
- `{nutrient}_serving` - Per serving
- `{nutrient}_unit` - Unit (g, mg, kcal, kJ)
- `{nutrient}_value` - Raw numeric value

**Core Nutriments:**

| Nutrient Key | Unit | Description |
|--------------|------|-------------|
| energy | kJ | Energy in kilojoules |
| energy-kcal | kcal | Energy in calories |
| fat | g | Total fat |
| saturated-fat | g | Saturated fatty acids |
| monounsaturated-fat | g | Monounsaturated fatty acids |
| polyunsaturated-fat | g | Polyunsaturated fatty acids |
| trans-fat | g | Trans fatty acids |
| cholesterol | mg | Cholesterol |
| carbohydrates | g | Total carbohydrates |
| sugars | g | Total sugars |
| fiber | g | Dietary fiber |
| proteins | g | Protein |
| salt | g | Salt |
| sodium | mg | Sodium |

**Extended Nutriments:**

| Nutrient Key | Unit | Description |
|--------------|------|-------------|
| vitamin-a | mcg | Vitamin A |
| vitamin-c | mg | Vitamin C |
| vitamin-d | mcg | Vitamin D |
| vitamin-e | mg | Vitamin E |
| vitamin-k | mcg | Vitamin K |
| vitamin-b1 | mg | Thiamin |
| vitamin-b2 | mg | Riboflavin |
| vitamin-b6 | mg | Vitamin B6 |
| vitamin-b12 | mcg | Vitamin B12 |
| vitamin-b9 | mcg | Folate |
| calcium | mg | Calcium |
| iron | mg | Iron |
| magnesium | mg | Magnesium |
| zinc | mg | Zinc |
| potassium | mg | Potassium |
| phosphorus | mg | Phosphorus |

**Sample Nutriments JSON:**
```json
{
  "nutriments": {
    "energy": 2092,
    "energy_100g": 2092,
    "energy_unit": "kJ",
    "energy-kcal": 500,
    "energy-kcal_100g": 500,
    "fat": 25,
    "fat_100g": 25,
    "fat_unit": "g",
    "saturated-fat": 22,
    "saturated-fat_100g": 22,
    "carbohydrates": 65,
    "carbohydrates_100g": 65,
    "sugars": 35,
    "sugars_100g": 35,
    "fiber": 3,
    "fiber_100g": 3,
    "proteins": 2,
    "proteins_100g": 2,
    "salt": 0.01,
    "salt_100g": 0.01,
    "sodium": 0.004,
    "sodium_100g": 0.004
  }
}
```

### Nutrition Scores

| Field | Type | Description |
|-------|------|-------------|
| nutriscore_grade | string | Nutri-Score grade (a-e) |
| nutriscore_score | integer | Nutri-Score numeric value |
| nutriscore_data | object | Detailed scoring components |
| nutrition_grade_fr | string | French nutrition grade |
| nutrition_grades | string | Nutrition grade value |
| nutrition_grades_tags | array | Grade tags |

**Nutri-Score Object (2023 version):**
```json
{
  "nutriscore": {
    "2023": {
      "grade": "c",
      "score": 8,
      "category_available": true,
      "data": {
        "negative_points": 15,
        "positive_points": 7,
        "energy": 5,
        "sugars": 6,
        "saturated_fat": 4,
        "sodium": 0,
        "fiber": 2,
        "proteins": 0,
        "fruits_vegetables_legumes": 5
      }
    }
  }
}
```

### Environmental Impact

| Field | Type | Description |
|-------|------|-------------|
| ecoscore_grade | string | Eco-Score grade (a-e) |
| ecoscore_score | integer | Eco-Score numeric value |
| ecoscore_data | object | Detailed environmental data |
| carbon_footprint | number | Carbon footprint (g CO2e) |
| carbon_footprint_100g | number | Per 100g |

**Eco-Score Data Structure:**
```json
{
  "ecoscore_data": {
    "grade": "c",
    "score": 45,
    "adjustments": {
      "origins_of_ingredients": {
        "value": -5,
        "origins_from_categories": ["en:france"]
      },
      "packaging": {
        "value": -10,
        "packagings": [
          {"material": "en:plastic", "shape": "en:bag"}
        ]
      },
      "production_system": {
        "value": 5,
        "labels": ["en:organic"]
      }
    },
    "agribalyse": {
      "code": "31032",
      "ef_total": 0.5
    }
  }
}
```

### NOVA Classification

| Field | Type | Description |
|-------|------|-------------|
| nova_group | integer | NOVA group (1-4) |
| nova_groups | string | NOVA group value |
| nova_groups_tags | array | NOVA group tags |

**NOVA Groups:**

| Group | Description |
|-------|-------------|
| 1 | Unprocessed or minimally processed foods |
| 2 | Processed culinary ingredients |
| 3 | Processed foods |
| 4 | Ultra-processed foods |

### Packaging Information

| Field | Type | Description |
|-------|------|-------------|
| packaging | string | Packaging description |
| packaging_tags | array | Normalized packaging tags |
| packagings | array | Detailed packaging objects |
| packaging_materials_tags | array | Material tags |
| packaging_recycling_tags | array | Recycling info |

**Packagings Array:**
```json
{
  "packagings": [
    {
      "material": "en:plastic",
      "shape": "en:bag",
      "recycling": "en:recycle",
      "number_of_units": 1,
      "quantity_per_unit": "170 g",
      "weight_measured": 5
    }
  ]
}
```

### Geographic and Origin Fields

| Field | Type | Description |
|-------|------|-------------|
| origins | string | Ingredient origins |
| origins_tags | array | Normalized origin tags |
| manufacturing_places | string | Manufacturing locations |
| manufacturing_places_tags | array | Normalized place tags |
| countries | string | Countries where sold |
| countries_tags | array | Normalized country tags |
| stores | string | Stores selling product |
| purchase_places | string | Purchase locations |

### Labels and Certifications

| Field | Type | Description |
|-------|------|-------------|
| labels | string | Product labels/certifications |
| labels_tags | array | Normalized label tags |
| labels_hierarchy | array | Hierarchical labels |
| emb_codes | string | Packaging/traceability codes |
| emb_codes_tags | array | Normalized EMB codes |

**Common Labels:**
```json
{
  "labels_tags": [
    "en:organic",
    "en:fair-trade",
    "en:vegan",
    "en:vegetarian",
    "en:gluten-free",
    "en:no-palm-oil",
    "en:nutriscore-grade-a"
  ]
}
```

### Image Fields

| Field | Type | Description |
|-------|------|-------------|
| image_url | string | Main product image |
| image_small_url | string | Thumbnail image |
| image_front_url | string | Front of package |
| image_ingredients_url | string | Ingredients image |
| image_nutrition_url | string | Nutrition facts image |
| images | object | All image metadata |

**Images Object:**
```json
{
  "images": {
    "front_en": {
      "imgid": "1",
      "sizes": {
        "100": {"h": 100, "w": 75},
        "200": {"h": 200, "w": 150},
        "400": {"h": 400, "w": 300},
        "full": {"h": 1200, "w": 900}
      },
      "uploaded_t": 1609459200,
      "uploader": "contributor123"
    }
  }
}
```

### Metadata Fields

| Field | Type | Description |
|-------|------|-------------|
| created_t | integer | Creation timestamp (UNIX) |
| created_datetime | string | Creation date (ISO8601) |
| last_modified_t | integer | Last modification timestamp |
| last_modified_datetime | string | Last modification date |
| creator | string | Original contributor |
| last_editor | string | Last editor username |
| editors_tags | array | All editors |
| completeness | number | Data completeness (0-1) |
| states | string | Product states |
| states_tags | array | State tags |

### Quality and Completeness

| Field | Type | Description |
|-------|------|-------------|
| completeness | number | Completeness score (0-1) |
| data_quality_tags | array | Data quality issues |
| data_quality_warnings_tags | array | Quality warnings |
| data_quality_bugs_tags | array | Known data issues |
| data_quality_errors_tags | array | Data errors |
| misc_tags | array | Miscellaneous tags |

**States Tags (completeness indicators):**
```json
{
  "states_tags": [
    "en:complete",
    "en:product-name-completed",
    "en:brands-completed",
    "en:categories-completed",
    "en:ingredients-completed",
    "en:nutrition-facts-completed",
    "en:photos-validated"
  ]
}
```

---

## Complete Product Example

```json
{
  "code": "737628064502",
  "product_name": "Organic Banana Chips",
  "product_name_en": "Organic Banana Chips",
  "generic_name": "Dried banana slices",
  "brands": "Trader Joe's",
  "brands_tags": ["trader-joe-s"],
  "categories": "Snacks, Dried fruits, Banana chips",
  "categories_tags": ["en:snacks", "en:dried-fruits", "en:banana-chips"],
  "categories_hierarchy": [
    "en:plant-based-foods-and-beverages",
    "en:plant-based-foods",
    "en:snacks",
    "en:dried-fruits",
    "en:banana-chips"
  ],
  "quantity": "170 g",

  "ingredients_text": "Organic bananas, organic coconut oil, organic sugar",
  "ingredients_text_en": "Organic bananas, organic coconut oil, organic sugar",
  "ingredients": [
    {
      "id": "en:banana",
      "text": "Organic bananas",
      "percent_estimate": 65,
      "vegan": "yes",
      "vegetarian": "yes"
    },
    {
      "id": "en:coconut-oil",
      "text": "Organic coconut oil",
      "percent_estimate": 25,
      "from_palm_oil": "no"
    },
    {
      "id": "en:sugar",
      "text": "Organic sugar",
      "percent_estimate": 10
    }
  ],
  "ingredients_n": 3,
  "additives_n": 0,
  "allergens": "",
  "allergens_tags": [],
  "traces": "tree nuts",
  "traces_tags": ["en:tree-nuts"],

  "nutriments": {
    "energy": 2092,
    "energy_100g": 2092,
    "energy-kcal": 500,
    "energy-kcal_100g": 500,
    "fat": 25,
    "fat_100g": 25,
    "saturated-fat": 22,
    "saturated-fat_100g": 22,
    "carbohydrates": 65,
    "carbohydrates_100g": 65,
    "sugars": 35,
    "sugars_100g": 35,
    "fiber": 3,
    "fiber_100g": 3,
    "proteins": 2,
    "proteins_100g": 2,
    "salt": 0.01,
    "salt_100g": 0.01
  },
  "nutrition_data_per": "100g",
  "serving_size": "30 g",
  "serving_quantity": 30,

  "nutriscore_grade": "d",
  "nutriscore_score": 14,
  "nutrition_grades": "d",
  "nova_group": 3,
  "nova_groups": "3",

  "ecoscore_grade": "b",
  "ecoscore_score": 65,

  "labels": "Organic, USDA Organic",
  "labels_tags": ["en:organic", "en:usda-organic"],

  "origins": "Philippines",
  "origins_tags": ["en:philippines"],
  "manufacturing_places": "USA",
  "countries": "United States",
  "countries_tags": ["en:united-states"],
  "stores": "Trader Joe's",

  "packaging": "Plastic bag",
  "packaging_tags": ["en:plastic", "en:bag"],

  "image_url": "https://images.openfoodfacts.org/images/products/073/762/806/4502/front_en.jpg",
  "image_small_url": "https://images.openfoodfacts.org/images/products/073/762/806/4502/front_en.100.jpg",

  "created_t": 1609459200,
  "created_datetime": "2021-01-01T00:00:00Z",
  "last_modified_t": 1704067200,
  "last_modified_datetime": "2024-01-01T00:00:00Z",
  "creator": "openfoodfacts-contributors",
  "completeness": 0.85
}
```

---

## Search API Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| search_terms | string | Full-text search query |
| categories_tags | string | Filter by category |
| brands_tags | string | Filter by brand |
| labels_tags | string | Filter by label |
| nutrition_grades | string | Filter by Nutri-Score (a-e) |
| nova_groups | integer | Filter by NOVA group (1-4) |
| ecoscore_grades | string | Filter by Eco-Score (a-e) |
| countries_tags | string | Filter by country |
| allergens_tags | string | Filter by allergen |
| traces_tags | string | Filter by traces |
| page | integer | Page number |
| page_size | integer | Results per page (max 100) |
| fields | string | Comma-separated fields to return |
| sort_by | string | Sort field |

**Example Search:**
```bash
curl "https://world.openfoodfacts.org/api/v2/search?categories_tags=breakfast-cereals&nutrition_grades=a&page_size=10&fields=code,product_name,nutriscore_grade,nutriments"
```

---

## Tag Normalization

Open Food Facts uses normalized tags for consistent querying:

| Raw Value | Normalized Tag |
|-----------|----------------|
| "Trader Joe's" | `trader-joe-s` |
| "United States" | `en:united-states` |
| "Organic" | `en:organic` |
| "Breakfast cereals" | `en:breakfast-cereals` |

**Tag Prefixes:**
- `en:` - English canonical tag
- `fr:` - French tag
- No prefix - User-contributed (may not be normalized)

---

## Data Quality Indicators

| Tag | Description |
|-----|-------------|
| `en:nutrition-not-filled-in` | Missing nutrition data |
| `en:ingredients-not-filled-in` | Missing ingredients |
| `en:quantity-not-filled-in` | Missing quantity |
| `en:categories-not-filled-in` | Missing categories |
| `en:product-name-not-filled-in` | Missing product name |
| `en:photos-to-be-validated` | Photos need review |
| `en:nutriscore-not-computed` | Nutri-Score not calculated |

---

## Rate Limits and Best Practices

| Aspect | Guideline |
|--------|-----------|
| Rate Limit | ~100 requests/minute |
| User-Agent | Required, identify your app |
| Caching | Cache responses when possible |
| Bulk Data | Use dumps for large-scale analysis |
| Updates | Use delta files for daily sync |

**Recommended User-Agent:**
```
MyApp/1.0 (myapp.com; contact@myapp.com)
```

---

## Integration Notes

### Mapping to Other Databases

| Open Food Facts | USDA FoodData | FooDB |
|-----------------|---------------|-------|
| code | fdcId | - |
| product_name | description | name |
| nutriments.energy-kcal | foodNutrients[energy] | - |
| nutriments.proteins | foodNutrients[protein] | - |
| ingredients | ingredientsList | - |

### Common Transformations

```python
# Convert nutriments to standard format
def normalize_nutriments(nutriments):
    return {
        'energy_kcal': nutriments.get('energy-kcal_100g', 0),
        'protein_g': nutriments.get('proteins_100g', 0),
        'fat_g': nutriments.get('fat_100g', 0),
        'carbs_g': nutriments.get('carbohydrates_100g', 0),
        'fiber_g': nutriments.get('fiber_100g', 0),
        'sugar_g': nutriments.get('sugars_100g', 0),
        'salt_g': nutriments.get('salt_100g', 0),
        'sodium_mg': nutriments.get('sodium_100g', 0) * 1000
    }
```

---

## Citation

```
Open Food Facts Contributors. Open Food Facts.
https://world.openfoodfacts.org

License: Open Database License (ODbL) v1.0
https://opendatacommons.org/licenses/odbl/1.0/
```

---

## Related Documents

- [foodb.md](./foodb.md) - Food composition database
- [usda-fooddata-central.md](./usda-fooddata-central.md) - USDA nutrition database
- [unified-schema-analysis.md](./unified-schema-analysis.md) - Cross-database schema mapping
