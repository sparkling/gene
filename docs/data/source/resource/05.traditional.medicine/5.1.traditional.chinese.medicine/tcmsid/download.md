---
id: download-tcmsid
title: "TCMSID Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# TCMSID Download Instructions

## Quick Start

```bash
# Visit the web interface for data access
# URL: http://lsp.nwu.edu.cn/tcmsid.php
```

## Prerequisites

- Web browser for data access
- Python 3.8+ for data processing (optional)
- Spreadsheet software for data analysis

## Download Methods

### Primary: Web Interface

Access TCMSID through the analysis platform:
- **URL**: http://lsp.nwu.edu.cn/tcmsid.php
- **Data Access**: Bulk data available via download section

### Available Downloads

TCMSID provides comprehensive ADME property data:

| File Type | Format | Description |
|-----------|--------|-------------|
| Compound Data | CSV/Excel | All compounds with ADME properties |
| Herb-Compound Links | CSV/Excel | Associations |
| Target Data | CSV/Excel | Protein targets |
| Network Files | Various | For pathway analysis |

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| TCM Herbs | 499 | Herb metadata |
| Compounds | 29,384 | With full ADME properties |
| Target Proteins | 3,311 | Predicted targets |
| Diseases | 837 | Disease associations |
| Compound-Target Pairs | 98,215 | Interaction predictions |

## ADME Property Fields

Each compound includes:

| Property | Description | Threshold |
|----------|-------------|-----------|
| OB (Oral Bioavailability) | Fraction reaching circulation | >= 30% favorable |
| DL (Drug-Likeness) | Similarity to known drugs | >= 0.18 favorable |
| Caco-2 | Intestinal permeability | >= -0.4 favorable |
| BBB | Blood-brain barrier penetration | >= -0.3 favorable |
| HL (Half-Life) | Pharmacokinetic duration | Long/Short |
| Lipinski | Rule of five compliance | 0-4 violations |

## Download Steps

1. Navigate to http://lsp.nwu.edu.cn/tcmsid.php
2. Access the download section
3. Select data type (compounds, herbs, targets)
4. Download CSV/Excel files
5. Filter by ADME thresholds as needed

## Filtering Workflow

Recommended compound prioritization:

```python
import pandas as pd

# Load TCMSID compound data
compounds = pd.read_csv('tcmsid_compounds.csv')

# Apply standard ADME filters
filtered = compounds[
    (compounds['OB'] >= 30) &      # Oral bioavailability
    (compounds['DL'] >= 0.18) &     # Drug-likeness
    (compounds['Caco2'] >= -0.4)    # Intestinal permeability
]

print(f"Drug-like compounds: {len(filtered)} of {len(compounds)}")
```

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Herbs | 499 |
| Compounds | 29,384 |
| Targets | 3,311 |
| Diseases | 837 |
| Compound-Target Pairs | 98,215 |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Periodic |
| Notification | Check website |

---

## Dataset Versions

### Current Release: TCMSID

| Property | Value |
|----------|-------|
| Version | 1.0 |
| Release Date | 2018-01-01 |
| Total Size | ~150 MB |
| Focus | ADME Properties |

### Version Contents

| Component | Records | Description |
|-----------|---------|-------------|
| TCM Herbs | 499 | Herb metadata |
| Compounds | 29,384 | With full ADME properties |
| Target Proteins | 3,311 | Predicted targets |
| Diseases | 837 | Disease associations |
| Compound-Target Pairs | 98,215 | Interaction predictions |

### ADME Properties Included

| Property | Threshold | Description |
|----------|-----------|-------------|
| OB | >= 30% | Oral Bioavailability |
| DL | >= 0.18 | Drug-Likeness |
| Caco-2 | >= -0.4 | Intestinal permeability |
| BBB | >= -0.3 | Blood-brain barrier |
| HL | Long/Short | Half-life |

---

## Notes

- All compounds have OB and DL scores
- ADME predictions are computational models
- Standard thresholds (OB>=30%, DL>=0.18) widely used
- Network visualization tools integrated
- Complements BATMAN-TCM with pharmacokinetic data
