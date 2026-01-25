---
id: exposome.explorer
title: "Exposome-Explorer"
type: data-source
category: nutrition
subcategory: metabolomics
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [exposome, biomarkers, diet, environment, metabolomics]
---

# Exposome-Explorer

**Category:** [Nutrition](../../../_index.md) > [Metabolomics](../_index.md)

## Overview

Exposome-Explorer is a manually curated database dedicated to biomarkers of dietary and environmental exposures. It systematically catalogs published associations between biomarkers measured in human samples and dietary or environmental exposures.

The database includes detailed information on biomarker concentrations, specimen types, analytical methods, and statistical associations with specific exposures. It covers dietary biomarkers (nutrients, phytochemicals, food additives) and environmental biomarkers (pollutants, toxins, occupational exposures).

Exposome-Explorer is essential for exposure assessment research, biomarker discovery, and epidemiological studies investigating diet-disease relationships.

## Key Statistics

| Metric | Value |
|--------|-------|
| Biomarkers | 900+ |
| Dietary Exposures | 600+ |
| Publications | 1,500+ |
| Last Update | 2023 |
| Coverage | Global studies |

## Primary Use Cases

1. Dietary biomarker identification
2. Exposure assessment methodology
3. Epidemiological study design
4. Metabolomics biomarker validation
5. Diet-disease association research

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Biomarker ID | Internal | EXP000123 |
| Exposure ID | Internal | DIET00456 |
| PubChem CID | Numeric | 5280343 |
| HMDB ID | `HMDB[0-9]{7}` | HMDB0000001 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://exposome-explorer.iarc.fr | Free access |
| Downloads | http://exposome-explorer.iarc.fr/downloads | TSV files |
| API | N/A | No public API |

## License

| Aspect | Value |
|--------|-------|
| License | Creative Commons Attribution 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [HMDB](../hmdb/README.md) - Human Metabolome Database
- [FooDB](../../6.1.food.composition/foodb/README.md) - Food compound database
