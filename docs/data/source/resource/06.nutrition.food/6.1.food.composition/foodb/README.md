---
id: foodb
title: "FooDB - The Food Database"
type: data-source
category: nutrition
subcategory: food.composition
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [food, compounds, chemistry, nutrients, phytochemicals]
---

# FooDB - The Food Database

**Category:** [Nutrition](../../../_index.md) > [Food Composition](../_index.md)

## Overview

FooDB is the world's largest and most comprehensive resource on food constituents, chemistry, and biology. It provides detailed information about the chemical compounds found in common foods, including their structures, concentrations, and physiological effects.

The database contains information on both macronutrients and micronutrients, as well as thousands of non-nutrient compounds including flavor compounds, colors, additives, and contaminants. Each compound is linked to detailed descriptions, concentrations in specific foods, and external database cross-references.

FooDB is particularly valuable for nutritional research, food chemistry, metabolomics studies, and understanding the complexity of human dietary exposure.

## Key Statistics

| Metric | Value |
|--------|-------|
| Foods | ~1,000 |
| Compounds | ~70,000 |
| Food-Compound Links | ~800,000 |
| Last Update | 2023 |
| Coverage | Common foods worldwide |

## Primary Use Cases

1. Identifying chemical compounds in specific foods
2. Understanding food chemistry and flavor profiles
3. Nutritional composition analysis
4. Metabolomics study design
5. Food safety and contaminant assessment

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| FooDB Compound ID | `FDB[0-9]{6}` | FDB000123 |
| FooDB Food ID | `FOOD[0-9]{5}` | FOOD00123 |
| PubChem CID | Numeric | 5280343 |
| HMDB ID | `HMDB[0-9]{7}` | HMDB0000001 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://foodb.ca | Browse and search |
| Downloads | https://foodb.ca/downloads | Full database dumps |
| API | N/A | No public API |

## License

| Aspect | Value |
|--------|-------|
| License | Creative Commons Attribution-NonCommercial 4.0 |
| Commercial Use | Requires license |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [HMDB](../../6.4.metabolomics/hmdb/_index.md) - Human Metabolome Database
- [USDA FoodData](../usda.fooddata/) - US nutrient database
