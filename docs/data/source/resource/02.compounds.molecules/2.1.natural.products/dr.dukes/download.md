---
id: download-dr.dukes
title: "Dr. Duke's Database Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Dr. Duke's Phytochemical Database Download Instructions

## Quick Start

```bash
# Download main dataset
wget https://ndownloader.figshare.com/files/43363335 -O dukes_phytochem.csv

# Download data dictionary
wget https://ndownloader.figshare.com/files/43363338 -O dukes_dictionary.csv
```

## Prerequisites

- **wget** or **curl** for downloads
- Spreadsheet software or database tools for analysis
- ~50 MB disk space

## No Registration Required

Data is CC0 (Public Domain) and freely downloadable without registration.

## Download Methods

### Method 1: Figshare Direct Download (Recommended)

```bash
# Main phytochemical dataset
wget https://ndownloader.figshare.com/files/43363335 -O dukes_phytochem.csv

# Data dictionary
wget https://ndownloader.figshare.com/files/43363338 -O dukes_dictionary.csv

# Alternative with curl
curl -L https://ndownloader.figshare.com/files/43363335 -o dukes_phytochem.csv
```

### Method 2: USDA Ag Data Commons

```bash
# Navigate to the Ag Data Commons archive
# https://agdatacommons.nal.usda.gov/

# Search for "Duke phytochemical" and download available files
# Files include biological activities, ethnobotanical uses, and chemical data
```

### Method 3: Interactive Web Interface

```bash
# Access the interactive search interface
# https://phytochem.nal.usda.gov/

# Use for targeted queries rather than bulk download
# Supports searches by plant, chemical, and activity
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| dukes_phytochem.csv | ~30 MB | Main chemical-plant associations |
| dukes_dictionary.csv | ~1 MB | Field definitions and codes |
| biological_activities.csv | ~5 MB | Activity data for compounds |
| ethnobotany_uses.csv | ~10 MB | Traditional use records |

## Data Categories Available

| Category | Content |
|----------|---------|
| Plants | Taxonomic information, plant parts |
| Chemicals | Compound names, CAS numbers |
| Activities | Biological effects, dosages, LD50 |
| Ethnobotany | Traditional uses by culture |
| Concentrations | PPM values in plant parts |

## Post-Download Processing

```bash
# Preview the data structure
head -20 dukes_phytochem.csv

# Count total records
wc -l dukes_phytochem.csv

# Extract unique plant species
cut -d',' -f2 dukes_phytochem.csv | sort | uniq | wc -l

# Load into SQLite for analysis
sqlite3 dukes.db << 'EOF'
.mode csv
.import dukes_phytochem.csv phytochem
.schema
SELECT COUNT(*) FROM phytochem;
EOF

# Filter for specific plant
grep -i "curcuma" dukes_phytochem.csv > turmeric_compounds.csv
```

## Verification

```bash
# Check file integrity
file dukes_phytochem.csv
# Should report: ASCII text or UTF-8 text

# Validate CSV structure
head -1 dukes_phytochem.csv | tr ',' '\n' | wc -l

# Count records by category
awk -F',' '{print $1}' dukes_phytochem.csv | sort | uniq -c
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| v1.9.12.6-Beta | 2017-10-18 | ~50 MB | Current |
| Legacy | 1992-2016 | Varies | Archived |

### Version Notes

Dr. Duke's Database is an archived static resource covering data from 1992-2016:
- Based on Dr. James Duke's compilations from his USDA work
- Released under CC0 public domain license
- No new updates expected (historical archive)

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://phytochem.nal.usda.gov` |
| Rate Limit | N/A (web interface) |
| Auth Required | No |
| Documentation | https://phytochem.nal.usda.gov |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database Status | Archived (2016) |
| New Updates | None (static archive) |
| Data Quality | Stable, curated |

## Common Issues

- **Encoding issues**: Files are UTF-8; ensure your tools support it
- **Field delimiters**: Standard comma-separated; some fields may contain commas in quotes
- **Missing values**: Empty fields are common; handle nulls appropriately
- **Historical data**: Some references may be to older literature

## Alternative Access via Web Search

For targeted queries rather than bulk download:

```bash
# Search by plant name
# https://phytochem.nal.usda.gov/phytochem/plants/list?q=Curcuma

# Search by chemical name
# https://phytochem.nal.usda.gov/phytochem/chemicals/list?q=curcumin

# Search by biological activity
# https://phytochem.nal.usda.gov/phytochem/activities/list?q=anti-inflammatory
```

## License

CC0 (Public Domain) - Free for any use without restrictions
