---
id: download-npact
title: "NPACT Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# NPACT Download Instructions

## Quick Start

```bash
# Visit the web interface for data access
# URL: http://crdd.osdd.net/raghava/npact/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ for data processing (optional)
- Spreadsheet software for Excel/CSV analysis

## Download Methods

### Primary: Web Interface

Access NPACT through the main portal:
- **URL**: http://crdd.osdd.net/raghava/npact/
- **Data Export**: Excel/CSV via web interface
- **Browse Options**: By compound, plant, cancer type

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| Anti-cancer Compounds | 1,574 | Experimentally validated |
| Plant Sources | 1,000+ | Traditional plant sources |
| Cancer Cell Lines | 150+ | Tested cell lines |
| Target Proteins | 200+ | Molecular targets |
| Literature References | 3,000+ | PubMed citations |
| Activity Data Points | 5,000+ | IC50/EC50 values |

## Data Structure

### Compound Entry Fields
- NPACT ID (internal)
- Compound name
- PubChem CID
- SMILES notation
- Molecular formula
- Molecular weight
- LogP
- Drug-likeness scores
- Source plants
- Activity data

### Activity Data Fields
- Activity type (IC50, EC50, GI50)
- Activity value (numeric)
- Activity unit (uM, nM)
- Cell line tested
- Cancer type
- PubMed reference

## Cancer Type Coverage

| Cancer Type | Compound Count |
|-------------|----------------|
| Breast | 400+ |
| Lung | 300+ |
| Colon | 250+ |
| Leukemia | 200+ |
| Liver | 150+ |
| Prostate | 150+ |

## Download Steps

1. Navigate to http://crdd.osdd.net/raghava/npact/
2. Use browse function (by compound, plant, or cancer)
3. Search for specific compounds or cancer types
4. Export search results as Excel/CSV
5. Compile comprehensive dataset

## Data Processing Example

```python
import pandas as pd

# Load NPACT compound data
compounds = pd.read_csv('npact_compounds.csv')

# Filter by cancer type
breast_cancer = compounds[compounds['cancer_type'] == 'Breast']
print(f"Breast cancer compounds: {len(breast_cancer)}")

# Filter by activity threshold
potent = compounds[
    (compounds['activity_type'] == 'IC50') &
    (compounds['activity_value'] <= 10) &
    (compounds['activity_unit'] == 'uM')
]
print(f"Potent compounds (IC50 <= 10 uM): {len(potent)}")

# Get compounds from specific plant
ashwagandha = compounds[
    compounds['plant_source'].str.contains('Withania', na=False)
]
print(f"Withania somnifera compounds: {len(ashwagandha)}")
```

## Key Strengths

- **Experimental Validation**: All compounds have tested activity
- **Quantitative Data**: IC50/EC50 values, not just predictions
- **Cancer Focus**: Specialized for oncology research
- **Indian Plants**: Strong coverage of Ayurvedic sources

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Anti-cancer Compounds | 1,574 |
| Plant Sources | 1,000+ |
| Cell Lines | 150+ |
| Targets | 200+ |
| References | 3,000+ |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Periodic |
| Notification | Check website |
| License | Academic use free |

## Notes

- Focus on anti-cancer activity limits scope
- Activity data from diverse assays (heterogeneous)
- No standardized target prediction method
- Updates may lag behind literature
- Contact maintainers for commercial licensing
- Complements IMPPAT for validated bioactivity
