---
id: usda.fooddata
title: "USDA FoodData Central"
type: data-source
category: compounds
subcategory: food-compounds
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [food-composition, nutrients, usda, dietary-reference, nutrition]
---

# USDA FoodData Central

**Category:** [Compounds & Molecules](../../_index.md) > [Food Compounds & Nutrients](../_index.md)

## Overview

USDA FoodData Central is the authoritative source for food composition data in the United States. It integrates data from multiple USDA data types including Foundation Foods (comprehensive nutrient profiles), SR Legacy (historical data), Branded Foods (industry-submitted data), and Survey Foods (NHANES data). The database provides nutrient data essential for dietary assessment, food labeling, and nutrition research.

FoodData Central contains nutrient profiles for thousands of foods, covering macronutrients, vitamins, minerals, and other food components. The data supports dietary intake assessment, food labeling compliance, nutrient database development, and nutrition research worldwide.

The database provides extensive API access and bulk download options, making it the foundation for many nutrition software applications and research databases.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Foods | 380,000+ |
| Foundation Foods | 2,100+ |
| SR Legacy Foods | 7,800+ |
| Branded Foods | 360,000+ |
| Nutrients Tracked | 150+ |
| Updates | Continuous |

## Primary Use Cases

1. **Dietary Assessment** - Calculate nutrient intake from food records
2. **Food Labeling** - Reference values for Nutrition Facts labels
3. **Research** - Source data for nutrition studies
4. **Software Development** - Foundation for nutrition apps
5. **Food Industry** - Product formulation and reformulation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| FDC ID | Integer | 167512 |
| NDB Number | 5 digits | 01001 (legacy) |
| GTIN/UPC | 12-14 digits | Branded foods |
| Food Description | Text | "Cheese, cheddar" |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://fdc.nal.usda.gov | Search and browse |
| REST API | https://api.nal.usda.gov/fdc | Requires API key (free) |
| Bulk Download | https://fdc.nal.usda.gov/download-datasets.html | CSV, JSON formats |
| Documentation | https://fdc.nal.usda.gov/api-guide.html | API reference |

## Data Types

| Type | Description |
|------|-------------|
| Foundation Foods | Comprehensive analytical data |
| SR Legacy | Standard Reference historical |
| Branded Foods | Industry product data |
| FNDDS | Survey foods for NHANES |
| Experimental Foods | Research data |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution Required | No |

## Related Resources

- [Phenol-Explorer](../phenol.explorer/_index.md) - Polyphenols
- [PhytoHub](../phytohub/_index.md) - Phytochemicals
- [PubChem](../../2.6.chemical.ontology.classification/pubchem/_index.md) - Chemical data

## References

1. U.S. Department of Agriculture, Agricultural Research Service. FoodData Central. https://fdc.nal.usda.gov/
