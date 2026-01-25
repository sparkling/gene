---
id: tcmbank
title: "TCMBank"
type: source
parent: ../README.md
tier: 2
status: active
category: traditional.medicine
subcategory: traditional.chinese.medicine
tags:
  - tcm
  - traditional-chinese-medicine
  - compound-library
  - drug-discovery
  - database
---

# TCMBank

**Category:** [Traditional Medicine](../../README.md) > [Traditional Chinese Medicine](../README.md)

## Overview

TCMBank is a comprehensive database designed to support drug discovery from Traditional Chinese Medicine. It provides a systematic compilation of TCM herbs, formulas, compounds, and their pharmacological properties, with a focus on structural data and drug-likeness assessment for identifying drug leads from natural sources.

The database emphasizes compound-centric information including chemical structures, physicochemical properties, ADMET predictions, and target predictions. This makes TCMBank particularly valuable for medicinal chemists and drug discovery researchers interested in natural product-derived lead compounds.

TCMBank integrates data from multiple TCM databases and literature sources, providing a unified resource for accessing TCM compound information with standardized chemical structures.

## Key Statistics

| Metric | Value |
|--------|-------|
| TCM Herbs | 9,000+ |
| Compounds | 60,000+ |
| Prescriptions/Formulas | 75,000+ |
| Target Proteins | 15,000+ |
| Diseases | 8,000+ |
| Compound-Target Pairs | 250,000+ |

## Primary Use Cases

1. Virtual screening of TCM compound libraries
2. Identifying drug-like natural products from TCM
3. Structure-based drug discovery from herbs
4. Lead compound optimization starting from TCM scaffolds
5. Target fishing for active TCM ingredients

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Compound ID | TCMBank internal | TCM001234 |
| Herb ID | TCMBank internal | HERB001234 |
| PubChem CID | Integer | 5280343 |
| SMILES | Standard notation | CC(=O)OC1=CC=CC=C1 |
| InChIKey | 27-character string | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://tcmbank.cn/ | Interactive search |
| Download | Compound SDF files | Structure download |
| Search | Multi-parameter | Structure, property, target |

## Data Model

```
Herbs (9,000+) -----> Formulas (75,000+)
     |
     v
Compounds (60,000+) -----> Drug-likeness Assessment
     |
     v
Targets (15,000+) -----> Diseases (8,000+)
```

## Compound Properties

| Property Category | Details |
|-------------------|---------|
| Structural | SMILES, InChI, SDF, MOL |
| Physicochemical | MW, logP, TPSA, HBD, HBA |
| Drug-likeness | Lipinski, Veber, PAINS |
| ADMET | Absorption, metabolism, toxicity predictions |
| Target | Predicted and known targets |

## Drug-Likeness Filters

| Filter | Criteria |
|--------|----------|
| Lipinski RO5 | MW <= 500, logP <= 5, HBD <= 5, HBA <= 10 |
| Veber | RotBonds <= 10, TPSA <= 140 |
| PAINS | No pan-assay interference compounds |
| Lead-like | MW 250-350, logP <= 3.5 |

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

## Limitations

- Focus on Chinese publications may limit coverage
- ADMET predictions are computational
- Regular updates needed for new literature

## See Also

- [BATMAN-TCM](../batman.tcm/README.md)
- [TCMSID](../tcmsid/README.md)
- [ETCM](../etcm/README.md)
