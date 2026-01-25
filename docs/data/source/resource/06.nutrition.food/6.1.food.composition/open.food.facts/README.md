---
id: open.food.facts
title: "Open Food Facts"
type: data-source
category: nutrition
subcategory: food.composition
parent: ../README.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [food, products, labels, nutrition, crowdsourced, open-data]
---

# Open Food Facts

**Category:** [Nutrition](../../../README.md) > [Food Composition](../README.md)

## Overview

Open Food Facts is a free, open, collaborative database of food products from around the world. It collects ingredients, allergens, nutrition facts, and all information found on product labels through crowdsourced contributions.

The database covers packaged food products globally, with extensive coverage in Europe and growing coverage worldwide. Each product entry includes barcode information, brand, categories, ingredients lists, Nutri-Score grades, NOVA processing classifications, and environmental impact scores.

Open Food Facts is particularly valuable for food product analysis, nutritional epidemiology, and consumer applications. Its open nature makes it ideal for research and app development.

## Key Statistics

| Metric | Value |
|--------|-------|
| Products | >3,000,000 |
| Countries | 180+ |
| Contributors | 30,000+ |
| Last Update | Continuous |
| Coverage | Packaged foods worldwide |

## Primary Use Cases

1. Product barcode lookup and identification
2. Nutritional analysis of commercial products
3. Allergen and ingredient tracking
4. Food product categorization
5. Environmental impact assessment (Eco-Score)

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Product Barcode | EAN-13/UPC | 3017620422003 |
| OFF ID | URL-based | nutella-ferrero-400g |
| Brand | Text | Ferrero |
| Categories | Hierarchical | en:spreads |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://world.openfoodfacts.org | Browse and search |
| API | https://world.openfoodfacts.org/api/v2 | REST API |
| Downloads | https://world.openfoodfacts.org/data | CSV, MongoDB dumps |
| Mobile Apps | iOS/Android | Barcode scanning |

## License

| Aspect | Value |
|--------|-------|
| License | Open Database License (ODbL) |
| Commercial Use | Yes |
| Attribution | Required |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [FooDB](../foodb/README.md) - Food compound database
- [USDA FoodData](../usda.fooddata/README.md) - US nutrient database
