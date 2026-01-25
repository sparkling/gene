---
id: hmdb
title: "HMDB - Human Metabolome Database"
type: data-source
category: nutrition
subcategory: metabolomics
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [metabolome, metabolites, human, biomarkers, spectroscopy]
---

# HMDB - Human Metabolome Database

**Category:** [Nutrition](../../../README.md) > [Metabolomics](../README.md)

## Overview

The Human Metabolome Database (HMDB) is the most comprehensive database of human metabolites and their associated biological, chemical, clinical, and molecular biology data. It contains detailed information about small molecule metabolites found in the human body.

HMDB includes metabolite concentration data for biofluid types (blood, urine, CSF, etc.), reference spectra (NMR, MS/MS), disease associations, pathway links, and extensive cross-references. It is designed to support metabolomics research, clinical chemistry, and biomarker discovery.

HMDB is a cornerstone resource for metabolomics research, providing gold-standard reference data for metabolite identification and quantification.

## Key Statistics

| Metric | Value |
|--------|-------|
| Metabolites | 220,000+ |
| With Concentrations | 5,000+ |
| MS/MS Spectra | 650,000+ |
| NMR Spectra | 5,000+ |
| Last Update | 2024 |
| Coverage | Human metabolome |

## Primary Use Cases

1. Metabolite identification from mass spectrometry
2. Biomarker discovery and validation
3. Metabolic pathway analysis
4. Clinical chemistry reference
5. Disease-metabolite association studies

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| HMDB ID | `HMDB[0-9]{7}` | HMDB0000001 |
| KEGG ID | `C[0-9]{5}` | C00031 |
| PubChem CID | Numeric | 5793 |
| ChEBI ID | `CHEBI:[0-9]+` | CHEBI:17234 |
| InChIKey | Standard | WQZGKKKJIJFFOK-... |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://hmdb.ca | Browse and search |
| API | https://hmdb.ca/api | REST API |
| Downloads | https://hmdb.ca/downloads | XML, SDF, CSV |
| MetaboAnalyst | https://www.metaboanalyst.ca | Analysis platform |

## License

| Aspect | Value |
|--------|-------|
| License | Creative Commons Attribution-NonCommercial 4.0 |
| Commercial Use | Requires license |
| Attribution | Required |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [Exposome-Explorer](../exposome.explorer/README.md) - Exposure biomarkers
- [FooDB](../../6.1.food.composition/foodb/README.md) - Food compound database
