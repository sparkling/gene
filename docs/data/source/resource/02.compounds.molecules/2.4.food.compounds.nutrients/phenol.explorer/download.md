---
id: download-phenol.explorer
title: "Phenol-Explorer Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Phenol-Explorer Download Instructions

## Quick Start

```bash
# Access via web interface
# http://phenol-explorer.eu/contents

# Data can be exported from search results
# Contact database maintainers for bulk access
```

## Prerequisites

- Web browser for interface access
- Spreadsheet software for exported data
- Contact maintainers for programmatic access

## Registration

Free for academic use. Contact for commercial access.

## Download Methods

### Method 1: Web Interface Export (Recommended)

```bash
# 1. Navigate to Phenol-Explorer
#    http://phenol-explorer.eu

# 2. Use Composition Explorer
#    http://phenol-explorer.eu/contents
#    - Select compound class or food group
#    - Browse food content data
#    - Export search results to Excel/CSV

# 3. Use Metabolite Explorer
#    http://phenol-explorer.eu/metabolites
#    - Browse metabolite data
#    - Export selected records
```

### Method 2: Food-Based Search

```bash
# 1. Go to: http://phenol-explorer.eu/contents/food/list
# 2. Select food of interest (e.g., "Apple")
# 3. View polyphenol content by compound class
# 4. Export data via "Download" button
```

### Method 3: Compound-Based Search

```bash
# 1. Go to: http://phenol-explorer.eu/compounds/list
# 2. Select compound class (e.g., Flavonoids > Flavonols)
# 3. Click specific compound (e.g., Quercetin)
# 4. View food sources and content values
# 5. Export food content data
```

### Method 4: Bulk Data Request

```bash
# For complete database access:
# 1. Contact: phenol-explorer@inra.fr
# 2. Describe intended use
# 3. Request appropriate data files

# Bulk files may include:
# - Complete composition database
# - Metabolite database
# - Pharmacokinetic data
```

## Data Categories

| Category | Content | Access |
|----------|---------|--------|
| Composition | Polyphenol content in foods | Web export |
| Metabolites | Human/microbial metabolites | Web export |
| Pharmacokinetics | Cmax, Tmax, AUC values | Web export |
| References | Literature citations | Web browse |

## File Formats

| Format | Description |
|--------|-------------|
| Excel (.xlsx) | Primary export format |
| CSV | Tab or comma-separated |
| PDF | Summary reports |

## Data Fields in Exports

### Composition Data

| Field | Description |
|-------|-------------|
| Food | Food name |
| Compound | Polyphenol name |
| Class | Polyphenol class |
| Mean | Mean content value |
| Min/Max | Range of values |
| Unit | mg/100g or mg/100mL |
| n | Number of data points |
| Method | Analytical method |
| Reference | Publication source |

### Metabolite Data

| Field | Description |
|-------|-------------|
| Parent | Parent compound |
| Metabolite | Metabolite name |
| Type | Phase I, II, Microbial |
| Biofluid | Detection matrix |
| Reference | Publication source |

## Post-Download Processing

```bash
# Convert Excel to CSV (if needed)
# Using LibreOffice:
libreoffice --headless --convert-to csv phenol_data.xlsx

# Or using Python:
python3 << 'EOF'
import pandas as pd
df = pd.read_excel('phenol_data.xlsx')
df.to_csv('phenol_data.csv', index=False)
EOF

# Basic analysis
# Count foods by compound class
awk -F',' 'NR>1 {print $3}' phenol_data.csv | sort | uniq -c | sort -rn

# Filter for specific compound
grep -i "quercetin" phenol_data.csv > quercetin_foods.csv

# Calculate average content
awk -F',' 'NR>1 && $4!="" {sum+=$4; n++} END {print sum/n}' phenol_data.csv
```

## Data Quality Notes

- Content values are from peer-reviewed publications
- Analytical methods documented for each value
- Min/max ranges show natural variation
- Data points (n) indicate sample size

## Verification

```bash
# Check exported file integrity
file phenol_data.xlsx
# Should be: Microsoft Excel 2007+

# Check CSV structure
head -5 phenol_data.csv

# Verify field count
head -1 phenol_data.csv | tr ',' '\n' | wc -l
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| v3.6 | 2015-06 | ~50 MB | Current |
| v3.0 | 2013 | ~40 MB | Archived |
| v2.0 | 2012 | ~30 MB | Archived |

### Version Notes

Phenol-Explorer 3.6 contains:
- 35,000+ content values for 500 polyphenols
- 400+ foods with polyphenol data
- Metabolism data for 375 metabolites
- Effects of food processing on polyphenol content

## API Access

| Property | Value |
|----------|-------|
| Base URL | `http://phenol-explorer.eu` |
| Rate Limit | N/A (web interface) |
| Auth Required | No |
| Documentation | Contact phenol-explorer@inra.fr |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | Periodic |
| Major versions | Every few years |
| Last major update | Version 3.6 (June 2015) |

## Common Issues

- **Export limits**: Web exports may be limited; request bulk for large datasets
- **Unit variations**: Ensure consistent units (mg/100g fresh weight standard)
- **Missing values**: Not all compounds measured in all foods
- **Method differences**: Consider analytical method when comparing values

## Integration Notes

Cross-reference with:
- **PubChem**: Via compound CID
- **USDA FoodData**: Via food descriptions
- **ChEBI**: Via compound structure
- **PhytoHub**: For metabolite mapping

## Citation Requirements

When using data, cite:
```
Rothwell JA, et al. (2013) "Phenol-Explorer 3.0: a major update of the
Phenol-Explorer database to incorporate data on the effects of food
processing on polyphenol content." Database (Oxford). 2013:bat070.
```

## License

Free for academic use. Contact for commercial licensing.
