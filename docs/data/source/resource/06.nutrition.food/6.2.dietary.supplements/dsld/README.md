---
id: dsld
title: "DSLD - Dietary Supplement Label Database"
type: data-source
category: nutrition
subcategory: dietary.supplements
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [supplements, labels, nih, ods, ingredients]
---

# DSLD - Dietary Supplement Label Database

**Category:** [Nutrition](../../../_index.md) > [Dietary Supplements](../_index.md)

## Overview

The Dietary Supplement Label Database (DSLD) is a project of the NIH Office of Dietary Supplements (ODS). It contains the full label contents of dietary supplement products marketed in the United States, including vitamins, minerals, herbs, botanicals, amino acids, and specialty supplements.

DSLD provides searchable access to supplement facts, ingredient lists, suggested use, and manufacturer information exactly as they appear on product labels. The database is updated regularly and represents one of the most comprehensive collections of US dietary supplement product data.

DSLD is essential for supplement research, regulatory analysis, ingredient surveillance, and understanding the US dietary supplement marketplace.

## Key Statistics

| Metric | Value |
|--------|-------|
| Products | 150,000+ |
| Unique Ingredients | 13,000+ |
| Brands | 4,000+ |
| Last Update | Monthly |
| Coverage | US marketed supplements |

## Primary Use Cases

1. Supplement ingredient surveillance
2. Market research and analysis
3. Regulatory compliance checking
4. Ingredient dosage assessment
5. Product formulation research

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| DSLD ID | Numeric | 123456 |
| UPC Barcode | Standard | 012345678901 |
| NHP ID | Canada | 80012345 |
| Brand Name | Text | Nature's Way |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://dsld.od.nih.gov | Search and browse |
| API | https://api.ods.od.nih.gov/dsld | REST API |
| Downloads | https://dsld.od.nih.gov/download | Database dumps |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [Natural Medicines](../natural.medicines/README.md) - Evidence-based database
- [ConsumerLab](../consumerlab/README.md) - Product testing
