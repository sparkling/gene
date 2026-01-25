---
id: foodb
title: "FooDB - The Food Database"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: food.compounds.nutrients
tags:
  - food
  - compounds
  - chemistry
  - nutrients
  - phytochemicals
---

# FooDB - The Food Database

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

## Limitations

- No public API available for programmatic access
- Commercial use requires license agreement
- Updates less frequent than some competitors
- Concentration data quality varies by source

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [USDA FoodData](../usda.fooddata/README.md) - US nutrient database
