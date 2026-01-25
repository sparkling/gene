---
id: download-dailymed
title: "DailyMed Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# DailyMed Download Instructions

## Quick Start

```bash
# Download all human prescription drug labels (SPL XML)
wget -r -np -nH --cut-dirs=4 \
  https://dailymed.nlm.nih.gov/dailymed/downloads/dm_spl_release_human_rx.zip
unzip dm_spl_release_human_rx.zip
```

## Prerequisites

- **wget** for downloads
- **unzip** for extraction
- ~50-100 GB disk space for full dataset
- XML parser for processing SPL files

## No Registration Required

Data is public domain (US Government) and freely downloadable.

## Download Methods

### Method 1: Bulk SPL Downloads (Recommended)

```bash
# Human prescription drugs
wget https://dailymed.nlm.nih.gov/dailymed/downloads/dm_spl_release_human_rx.zip

# Over-the-counter drugs
wget https://dailymed.nlm.nih.gov/dailymed/downloads/dm_spl_release_human_otc.zip

# Animal drugs
wget https://dailymed.nlm.nih.gov/dailymed/downloads/dm_spl_release_animal.zip

# Homeopathic products
wget https://dailymed.nlm.nih.gov/dailymed/downloads/dm_spl_release_homeopathic.zip

# Extract all
for f in dm_spl_release_*.zip; do
  unzip "$f" -d "${f%.zip}"
done
```

### Method 2: REST API (Individual Labels)

```bash
# Get label by drug name
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name=aspirin" | jq

# Get label by NDC
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?ndc=0002-3228-30" | jq

# Get label by Set ID
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/12345678-1234-1234-1234-123456789012.json" | jq

# Get label PDF
curl -O "https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid={set_id}&type=display"
```

### Method 3: FTP Access

```bash
# Connect via FTP
ftp ftp.nlm.nih.gov
# Navigate to: /nlmdata/dailymed/

# Or use wget for specific files
wget -r -np ftp://ftp.nlm.nih.gov/nlmdata/dailymed/
```

### Method 4: RSS Feeds for Updates

```bash
# Subscribe to RSS for new labels
# https://dailymed.nlm.nih.gov/dailymed/rss.cfm

# Get recent updates via API
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?updated_since=2024-01-01" | jq
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| dm_spl_release_human_rx.zip | ~40 GB | Human Rx labels |
| dm_spl_release_human_otc.zip | ~15 GB | OTC labels |
| dm_spl_release_animal.zip | ~2 GB | Animal drug labels |
| dm_spl_release_homeopathic.zip | ~5 GB | Homeopathic labels |
| dm_spl_daily_update_*.zip | ~100 MB | Daily incremental updates |

## SPL Package Contents

Each label package contains:

| File | Description |
|------|-------------|
| *.xml | SPL XML document |
| *.jpg/*.png | Product images |
| *.pdf | Rendered label (some packages) |

## Post-Download Processing

```bash
# List all SPL XML files
find dm_spl_release_human_rx -name "*.xml" > spl_files.txt
wc -l spl_files.txt

# Extract drug names from SPL files
for xml in $(head -100 spl_files.txt); do
  xmllint --xpath "//title/text()" "$xml" 2>/dev/null
done

# Extract specific section (indications)
xmllint --xpath "//section[code/@code='34067-9']/text" label.xml

# Parse to JSON (using xsltproc or custom tool)
xsltproc spl2json.xsl label.xml > label.json

# Build searchable index with SQLite
sqlite3 dailymed.db << 'EOF'
CREATE TABLE labels (
  set_id TEXT PRIMARY KEY,
  title TEXT,
  effective_date TEXT,
  ndc_codes TEXT
);
CREATE INDEX idx_title ON labels(title);
EOF

# Extract NDC codes from all labels
grep -ohP 'NDC [0-9]{4,5}-[0-9]{3,4}-[0-9]{1,2}' *.xml | sort | uniq
```

## Verification

```bash
# Check archive integrity
unzip -t dm_spl_release_human_rx.zip

# Count labels
find dm_spl_release_human_rx -name "*.xml" | wc -l

# Validate XML structure
xmllint --noout label.xml

# Check for required elements
xmllint --xpath "//id/@root" label.xml
xmllint --xpath "//title/text()" label.xml
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| Dec 2025 | 2025-12-19 | ~60 GB total | Current |
| Monthly releases | First Monday | Varies | Archived |

### Version Notes

DailyMed contains 154,512+ labeling records submitted to FDA:
- Prescription and nonprescription drugs
- Human and animal drugs
- Medical gases, devices, cosmetics, dietary supplements

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://dailymed.nlm.nih.gov/dailymed/services/v2` |
| Rate Limit | No hard limit |
| Auth Required | No |
| Documentation | https://dailymed.nlm.nih.gov/dailymed/app-support-web-services.cfm |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Full releases | Monthly |
| Daily updates | Daily |
| Real-time via API | Continuous |

## API Endpoints Reference

| Endpoint | Description |
|----------|-------------|
| /v2/spls.json | Search labels |
| /v2/drugnames.json | List drug names |
| /v2/drugclasses.json | Drug classes |
| /v2/ndcs.json | NDC lookup |
| /v2/rxcuis.json | RxNorm mapping |

## Common Issues

- **Large archives**: Use streaming extraction for limited disk space
- **XML parsing**: Use robust XML parser (some labels have complex structure)
- **Character encoding**: Files are UTF-8; handle special characters properly
- **Historical labels**: Superseded versions available via archive
- **Missing images**: Some packages may lack product images

## Integration with Other Resources

```bash
# Link to RxNorm via RXCUI
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?rxcui=161" | jq

# Cross-reference to openFDA
curl "https://api.fda.gov/drug/label.json?search=openfda.product_ndc:0002-3228" | jq
```

## License

Public Domain (US Government) - Free for any use
