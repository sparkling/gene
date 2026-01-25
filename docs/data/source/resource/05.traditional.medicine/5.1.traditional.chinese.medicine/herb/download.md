---
id: download-herb
title: "HERB Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# HERB Download Instructions

## Quick Start

```bash
# Visit the web interface for bulk downloads
# URL: http://herb.ac.cn/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ for data processing (optional)
- Sufficient storage for expression data (~500MB)

## Download Methods

### Primary: Web Interface

Access the HERB database through the main portal:
- **URL**: http://herb.ac.cn/
- **Data Access**: Download section available after browsing

### Available Data Files

HERB provides bulk download of herb-ingredient-target associations.

| File Type | Format | Description |
|-----------|--------|-------------|
| Herb Data | TSV/Excel | 1,037 herbs with metadata |
| Ingredients | TSV/Excel | 12,933 ingredients with structures |
| Targets | TSV/Excel | 2,064 gene targets |
| Herb-Ingredient Links | TSV | 49,258 associations |
| Ingredient-Target Links | TSV | 28,212 associations |
| Expression Data | TSV | GEO-derived signatures |

## File Inventory

| File | Records | Description |
|------|---------|-------------|
| herbs.tsv | 1,037 | Herb metadata and IDs |
| ingredients.tsv | 12,933 | Chemical compounds |
| targets.tsv | 2,064 | Gene/protein targets |
| herb_ingredient.tsv | 49,258 | Herb-compound associations |
| ingredient_target.tsv | 28,212 | Compound-target links |
| diseases.tsv | 866 | Disease associations |
| expression/*.tsv | 2,000+ | Gene expression experiments |

## Download Steps

1. Navigate to http://herb.ac.cn/
2. Click on "Download" or "Data" section
3. Select desired data tables
4. Export as TSV or Excel format

## Verification

After download, verify record counts:

```python
import pandas as pd

# Check herb count
herbs = pd.read_csv('herbs.tsv', sep='\t')
assert len(herbs) >= 1037, "Expected at least 1037 herbs"

# Check ingredient count
ingredients = pd.read_csv('ingredients.tsv', sep='\t')
assert len(ingredients) >= 12933, "Expected at least 12933 ingredients"

print(f"Herbs: {len(herbs)}")
print(f"Ingredients: {len(ingredients)}")
```

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Periodic with HERB versions |
| Current Version | HERB 2.0 |
| Notification | Check website for updates |

---

## Dataset Versions

### Current Release: HERB 2.0

| Property | Value |
|----------|-------|
| Version | 2.0 |
| Release Date | 2021-01-01 |
| Total Size | ~500 MB |
| Expression Data | 2,000+ experiments |

### Version Contents

| Component | Records | Description |
|-----------|---------|-------------|
| herbs.tsv | 1,037 | Herb metadata |
| ingredients.tsv | 12,933 | Chemical compounds |
| targets.tsv | 2,064 | Gene/protein targets |
| herb_ingredient.tsv | 49,258 | Herb-compound associations |
| ingredient_target.tsv | 28,212 | Compound-target links |
| diseases.tsv | 866 | Disease associations |

### Previous Versions

| Version | Release | Status |
|---------|---------|--------|
| HERB 1.0 | 2020-01-01 | Archived |

---

## Notes

- Gene expression data is linked to GEO accessions
- Expression signatures enable connectivity analysis
- Large downloads may require stable connection
- Contact maintainers for commercial use inquiries
