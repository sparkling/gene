---
id: etcm
title: "ETCM - Encyclopedia of Traditional Chinese Medicine"
type: source
parent: ../README.md
tier: 2
status: active
category: traditional.medicine
subcategory: traditional.chinese.medicine
tags:
  - tcm
  - traditional-chinese-medicine
  - herb
  - formula
  - pharmacology
---

# ETCM - Encyclopedia of Traditional Chinese Medicine

**Category:** [Traditional Medicine](../../README.md) > [Traditional Chinese Medicine](../README.md)

## Overview

ETCM (Encyclopedia of Traditional Chinese Medicine) is a comprehensive database integrating diverse TCM information including herbs, formulas, and their associated compounds, targets, and diseases. The database provides standardized, structured data following Traditional Chinese Medicine theory while bridging to modern pharmacological understanding.

ETCM emphasizes the systematic organization of TCM knowledge according to classical theory, including herb properties (nature, flavor, meridian tropism), formula composition principles, and therapeutic indications. This makes it particularly valuable for researchers seeking to understand TCM from both traditional and modern scientific perspectives.

The database includes quantitative data on herb-compound relationships, predicted drug-target interactions, and disease associations derived from text mining and computational predictions.

## Key Statistics

| Metric | Value |
|--------|-------|
| TCM Herbs | 403 |
| TCM Formulas | 3,677 |
| Compounds | 7,274 |
| Targets | 3,000+ |
| Diseases | 500+ |
| Herb-Compound Links | 23,000+ |

## Primary Use Cases

1. Understanding TCM formulas and their composition principles
2. Exploring herb properties according to TCM theory
3. Identifying active compounds in traditional herbs
4. Target prediction for TCM compounds
5. Disease-TCM association discovery

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Herb ID | ETCM internal | H001 |
| Formula ID | ETCM internal | F001 |
| Compound ID | ETCM internal / PubChem CID | C001 / 5280343 |
| Target ID | UniProt accession | P12345 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://www.tcmip.cn/ETCM/ | Interactive browsing |
| Download | Via web interface | Excel/CSV exports |
| Search | Multi-field search | Herb, formula, compound |

## Data Model

```
TCM Theory Properties
        |
    Herbs (403) -----> Formulas (3,677)
        |                    |
        v                    v
  Compounds (7,274) <--------+
        |
        v
   Targets (3,000+)
        |
        v
   Diseases (500+)
```

## TCM Properties Captured

| Property Type | Description |
|---------------|-------------|
| Nature (Xing) | Cold, Cool, Neutral, Warm, Hot |
| Flavor (Wei) | Sour, Bitter, Sweet, Pungent, Salty |
| Meridian Tropism | Target organ systems (Lung, Heart, Liver, etc.) |
| Actions | Traditional therapeutic effects |
| Indications | Disease patterns treated |

## License

| Aspect | Value |
|--------|-------|
| License | Academic use |
| Commercial Use | Contact maintainers |
| Attribution | Required for publications |

## Limitations

- No public REST API
- Manual download required for bulk data
- Less comprehensive than BATMAN-TCM for predicted interactions

## See Also

- [BATMAN-TCM](../batman.tcm/README.md)
- [TCMBank](../tcmbank/README.md)
- [HERB Database](../herb/README.md)
