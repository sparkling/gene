---
id: download-rxnorm
title: "RxNorm Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# RxNorm Download Instructions

## Quick Start

```bash
# Via RxNav API (no registration needed)
curl "https://rxnav.nlm.nih.gov/REST/rxcui/2670/allrelated.json" | jq

# Full download requires UMLS license (free registration)
```

## Prerequisites

- **UMLS License** (free) for full data downloads
- **curl** for API access
- ~5 GB disk space for full dataset
- RRF parser or database software

## Registration Requirements

- **API Access**: No registration required
- **Full Download**: Free UMLS license required
- Register at: https://uts.nlm.nih.gov/uts/signup-login

## Download Methods

### Method 1: RxNav REST API (Recommended for Targeted Access)

```bash
# Get drug information by RXCUI
curl "https://rxnav.nlm.nih.gov/REST/rxcui/2670.json" | jq

# Search by drug name
curl "https://rxnav.nlm.nih.gov/REST/drugs.json?name=acetaminophen" | jq

# Get all related concepts
curl "https://rxnav.nlm.nih.gov/REST/rxcui/2670/allrelated.json" | jq

# Get NDC codes for RXCUI
curl "https://rxnav.nlm.nih.gov/REST/rxcui/2670/ndcs.json" | jq

# Approximate term matching
curl "https://rxnav.nlm.nih.gov/REST/approximateTerm.json?term=tylenol%20500mg" | jq
```

### Method 2: UMLS Download (Full Dataset)

```bash
# 1. Login to UMLS at: https://www.nlm.nih.gov/research/umls/rxnorm/
# 2. Navigate to RxNorm files
# 3. Download current monthly release

# After download:
unzip rxnorm_full_*.zip -d rxnorm

# Key files:
# - RXNCONSO.RRF (concepts)
# - RXNREL.RRF (relationships)
# - RXNSAT.RRF (attributes)
# - RXNSTY.RRF (semantic types)
```

### Method 3: RxMix Batch Processing

```bash
# RxMix for batch operations
# https://mor.nlm.nih.gov/RxMix/

# Create batch file with drug names
echo "acetaminophen" > drugs.txt
echo "ibuprofen" >> drugs.txt
echo "aspirin" >> drugs.txt

# Submit via RxMix web interface for bulk processing
```

### Method 4: RxClass (Drug Classes)

```bash
# Get drugs by class
curl "https://rxnav.nlm.nih.gov/REST/rxclass/classMembers.json?classId=N02BE01&relaSource=ATC" | jq

# Get all drug classes
curl "https://rxnav.nlm.nih.gov/REST/rxclass/allClasses.json" | jq
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| RXNCONSO.RRF | ~2 GB | Concept names and sources |
| RXNREL.RRF | ~1 GB | Relationships |
| RXNSAT.RRF | ~800 MB | Attributes |
| RXNSTY.RRF | ~20 MB | Semantic types |
| RXNDOC.RRF | ~1 MB | Documentation |

## RRF Format

Files are pipe (|) delimited:

```
RXCUI|LAT|TS|LUI|STT|SUI|ISPREF|RXAUI|SAUI|SCUI|SDUI|SAB|TTY|CODE|STR|SRL|SUPPRESS|CVF
2670|ENG|S|L0014479|PF|S0042379|Y|A0043984||2670||RXNORM|SCD|2670|Acetaminophen 500 MG Oral Tablet|||
```

## Post-Download Processing

```bash
# Load into SQLite
sqlite3 rxnorm.db << 'EOF'
CREATE TABLE RXNCONSO (
  RXCUI TEXT, LAT TEXT, TS TEXT, LUI TEXT, STT TEXT,
  SUI TEXT, ISPREF TEXT, RXAUI TEXT PRIMARY KEY, SAUI TEXT,
  SCUI TEXT, SDUI TEXT, SAB TEXT, TTY TEXT, CODE TEXT,
  STR TEXT, SRL TEXT, SUPPRESS TEXT, CVF TEXT
);
.mode csv
.separator |
.import RXNCONSO.RRF RXNCONSO
CREATE INDEX idx_rxcui ON RXNCONSO(RXCUI);
CREATE INDEX idx_str ON RXNCONSO(STR);
CREATE INDEX idx_tty ON RXNCONSO(TTY);
EOF

# Query for specific term types
sqlite3 rxnorm.db "SELECT RXCUI, STR FROM RXNCONSO WHERE TTY='SCD' LIMIT 10;"

# Find all ingredients
sqlite3 rxnorm.db "SELECT DISTINCT STR FROM RXNCONSO WHERE TTY='IN';" > ingredients.txt

# Count concepts by type
sqlite3 rxnorm.db "SELECT TTY, COUNT(*) FROM RXNCONSO GROUP BY TTY ORDER BY COUNT(*) DESC;"
```

## Build Drug Dictionary

```bash
# Extract normalized drug names
sqlite3 rxnorm.db << 'EOF'
SELECT DISTINCT
  c1.RXCUI as rxcui,
  c1.STR as normalized_name,
  GROUP_CONCAT(DISTINCT c2.STR) as synonyms
FROM RXNCONSO c1
LEFT JOIN RXNCONSO c2 ON c1.RXCUI = c2.RXCUI
WHERE c1.SAB = 'RXNORM'
  AND c1.TTY IN ('SCD', 'SBD')
  AND c1.SUPPRESS = 'N'
GROUP BY c1.RXCUI
EOF
```

## Verification

```bash
# Check file integrity
wc -l RXNCONSO.RRF

# Verify field count
head -1 RXNCONSO.RRF | tr '|' '\n' | wc -l

# Check RXCUI uniqueness per atom
cut -d'|' -f8 RXNCONSO.RRF | sort | uniq -d | wc -l
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| RxNorm Dec 2025 | 2025-12-01 | ~5 GB | Current |
| RxNorm Nov 2025 | 2025-11-03 | ~5 GB | Archived |
| Previous monthly | First Monday each month | ~5 GB | Archived |

### Version Notes

RxNorm December 2025 edition includes:
- 14,592 base ingredients
- 5,151 brand names
- 17,544 clinical drugs
- 9,721 branded drugs
- 645 generic packs and 747 branded packs
- RxNorm celebrates its 20th anniversary in 2025

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://rxnav.nlm.nih.gov/REST` |
| Rate Limit | 20 req/sec (anonymous) |
| Auth Required | No (for API), Yes (for bulk download) |
| Documentation | https://lhncbc.nlm.nih.gov/RxNav |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Full release | Monthly (first Monday) |
| API data | Weekly (Wednesdays) |
| Weekly updates | Available |

## API Rate Limits

| Tier | Limit |
|------|-------|
| Anonymous | 20 requests/second |
| Registered | Higher limits available |

## Common Issues

- **UMLS registration**: Allow 1-2 days for license approval
- **Large files**: Use streaming parsers for memory efficiency
- **Pipe delimiter**: Handle empty fields carefully
- **Suppressible concepts**: Filter SUPPRESS='N' for active concepts
- **Multiple sources**: Use SAB='RXNORM' for RxNorm-native concepts

## Integration Mappings

```bash
# Map RxNorm to NDC
curl "https://rxnav.nlm.nih.gov/REST/rxcui/2670/ndcs.json" | jq

# Map RxNorm to ATC
curl "https://rxnav.nlm.nih.gov/REST/rxcui/2670/property.json?propName=ATC" | jq

# Map drug name to RXCUI
curl "https://rxnav.nlm.nih.gov/REST/rxcui.json?name=acetaminophen" | jq
```

## License

UMLS License (free) - Attribution required, research/clinical use
