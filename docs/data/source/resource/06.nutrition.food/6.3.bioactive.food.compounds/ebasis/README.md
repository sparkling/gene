---
id: ebasis
title: "eBASIS - Bioactive Substances in Food Information System"
type: data-source
category: nutrition
subcategory: bioactive.food.compounds
parent: ../README.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [bioactive, phytochemicals, food, europe, composition]
---

# eBASIS - Bioactive Substances in Food Information System

**Category:** [Nutrition](../../../README.md) > [Bioactive Food Compounds](../README.md)

## Overview

eBASIS (Bioactive Substances in Food Information System) is a European Food Information Resource (EuroFIR) database that compiles quality-evaluated composition data on bioactive compounds in plant-based foods. It focuses on compounds with potential health benefits beyond basic nutrition.

The database contains systematic data on bioactive compound concentrations in foods, including polyphenols, carotenoids, phytosterols, glucosinolates, and other phytochemicals. Each data point includes quality evaluation using standardized EuroFIR criteria.

eBASIS is particularly valuable for exposure assessment, dietary intake calculations, and research on bioactive compounds and their health effects.

## Key Statistics

| Metric | Value |
|--------|-------|
| Compounds | 300+ bioactives |
| Foods | 500+ plant foods |
| Data Points | 25,000+ |
| Last Update | 2022 |
| Coverage | European plant foods |

## Primary Use Cases

1. Dietary bioactive exposure assessment
2. Food composition research
3. Health claims substantiation
4. Functional food development
5. Epidemiological study support

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| eBASIS Compound ID | Internal | EB001234 |
| EuroFIR Food Code | LanguaL-based | A0148 |
| Component Group | Text | Flavonoids |
| PubChem CID | Numeric | 5280343 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://ebasis.eurofir.org | Registration required |
| API | N/A | No public API |
| Downloads | Contact EuroFIR | Restricted access |

## License

| Aspect | Value |
|--------|-------|
| License | EuroFIR Agreement |
| Commercial Use | Requires license |
| Attribution | Required |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [FooDB](../../6.1.food.composition/foodb/README.md) - Food compound database
- [Phenol-Explorer](../phenol.explorer/README.md) - Polyphenol database
