---
id: download-etcm
title: "ETCM Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# ETCM Download Instructions

## Quick Start

```bash
# Visit the web interface - manual export required
# URL: http://www.tcmip.cn/ETCM/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ for data processing (optional)
- Spreadsheet software for Excel exports

## Download Methods

### Primary: Web Interface

Access the ETCM database through the main portal:
- **URL**: http://www.tcmip.cn/ETCM/
- **Data Access**: Export via search results

### Export Options

ETCM provides data export through the web interface:

| Export Type | Format | Description |
|-------------|--------|-------------|
| Search Results | Excel/CSV | Query-specific exports |
| Herb Data | Excel | Individual herb profiles |
| Formula Data | Excel | Formula compositions |
| Compound Data | Excel | Chemical information |

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| TCM Herbs | 403 | Herb metadata and TCM properties |
| TCM Formulas | 3,677 | Formula compositions |
| Compounds | 7,274 | Chemical compounds |
| Targets | 3,000+ | Predicted target proteins |
| Diseases | 500+ | Disease associations |
| Herb-Compound Links | 23,000+ | Associations |

## Expected Data Structure

### Herb Data Fields
- Herb ID (ETCM internal)
- Chinese name
- Pinyin name
- English name
- Latin botanical name
- TCM Nature (cold/cool/neutral/warm/hot)
- TCM Flavor (sweet/bitter/sour/salty/pungent)
- Meridian Tropism (organ systems)
- Therapeutic actions

### Compound Data Fields
- Compound ID (ETCM internal / PubChem CID)
- Compound name
- SMILES notation
- Molecular formula
- Molecular weight
- Source herbs

## Download Steps

1. Navigate to http://www.tcmip.cn/ETCM/
2. Use search function to query herbs, formulas, or compounds
3. Browse to individual entity pages
4. Export data via available download options
5. Combine exports for comprehensive dataset

## Verification

After download, verify against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Herbs | 403 |
| Formulas | 3,677 |
| Compounds | 7,274 |
| Targets | 3,000+ |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Periodic |
| Notification | Check website |

## Notes

- No public REST API available
- Manual export required for bulk data
- TCM properties (nature, flavor, meridian) are unique to this database
- Formula composition follows classical TCM theory
- Contact maintainers for academic collaborations
