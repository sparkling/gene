---
id: npact
title: "NPACT - Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target"
type: source
parent: ../README.md
tier: 2
status: active
category: traditional.medicine
subcategory: south.east.asian.systems
tags:
  - natural-products
  - anti-cancer
  - indian-plants
  - bioactivity
  - drug-discovery
---

# NPACT - Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target

**Category:** [Traditional Medicine](../../README.md) > [South East Asian Systems](../README.md)

## Overview

NPACT (Naturally occurring Plant-based Anti-cancer Compound-activity-Target database) is a curated resource focusing specifically on natural compounds with demonstrated anti-cancer activity. It compiles experimentally validated anti-cancer compounds from plant sources, with emphasis on species used in Indian traditional medicine.

Unlike prediction-based databases, NPACT prioritizes experimentally confirmed bioactivity data including IC50, EC50, and other quantitative measures of anti-cancer effects. Each entry includes the cancer cell lines tested, activity values, and literature references for validation.

The database serves as a critical resource for cancer drug discovery from natural sources, providing validated starting points for lead optimization.

## Key Statistics

| Metric | Value |
|--------|-------|
| Anti-cancer Compounds | 1,574 |
| Plant Sources | 1,000+ |
| Cancer Cell Lines | 150+ |
| Target Proteins | 200+ |
| Literature References | 3,000+ |
| Activity Data Points | 5,000+ |

## Primary Use Cases

1. Identifying validated anti-cancer natural products
2. Finding lead compounds for cancer drug discovery
3. Exploring plant sources for bioactive molecules
4. Correlating structural features with anti-cancer activity
5. Target identification for natural anti-cancer agents

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| NPACT ID | Internal | NPACT001234 |
| Compound Name | Text | Curcumin |
| PubChem CID | Integer | 969516 |
| Cell Line | Standard name | MCF-7, HeLa |
| Activity Type | IC50/EC50/GI50 | IC50 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://crdd.osdd.net/raghava/npact/ | Interactive search |
| Download | Via web interface | Excel/CSV |
| Browse | By compound, plant, cancer type | Structured navigation |

## Data Model

```
Plant Sources (1,000+)
        |
        v
Anti-cancer Compounds (1,574) ---> Chemical Properties
        |                          - SMILES, MW, logP
        |                          - Drug-likeness
        v
Cancer Cell Lines (150+) <-----> Activity Data (IC50, etc.)
        |
        v
Target Proteins (200+)
```

## Activity Data Fields

| Field | Description | Example |
|-------|-------------|---------|
| Activity Type | Measurement type | IC50, EC50, GI50 |
| Activity Value | Numeric value | 2.5 |
| Activity Unit | Concentration unit | uM, nM |
| Cell Line | Cancer cell tested | MCF-7 (breast) |
| Reference | PubMed ID | PMID:12345678 |

## Cancer Type Coverage

| Cancer Type | Compound Count |
|-------------|----------------|
| Breast | 400+ |
| Lung | 300+ |
| Colon | 250+ |
| Leukemia | 200+ |
| Liver | 150+ |
| Prostate | 150+ |

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

## Key Strengths

- **Experimental Validation**: All compounds have tested activity
- **Quantitative Data**: IC50/EC50 values, not just predictions
- **Cancer Focus**: Specialized for oncology research
- **Indian Plants**: Strong coverage of Ayurvedic sources

## Limitations

- Focus on anti-cancer limits scope
- Activity data from diverse assays (heterogeneous)
- No standardized target prediction
- Updates may lag behind literature

## See Also

- [IMPPAT](../imppat/README.md)
- [TM-MC](../tm.mc/README.md)
- [CMAUP](../../5.4.multi.system.integration/{hit.2.0}/README.md)
