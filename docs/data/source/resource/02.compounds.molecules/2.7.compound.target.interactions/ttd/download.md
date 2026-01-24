---
id: download-ttd
title: "TTD Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# TTD Download Instructions

## Quick Start

```bash
# Download full data from TTD website
# Navigate to: https://idrblab.net/ttd/full-data-download

# Individual files available for targets, drugs, and mappings
```

## Prerequisites

- Web browser for download page
- ~200 MB disk space
- Tab-separated file parser

## Registration

Free for academic use. Limited commercial access.

## Download Methods

### Method 1: Full Data Download (Recommended)

```bash
# Navigate to TTD download page
# https://idrblab.net/ttd/full-data-download

# Available downloads:
# 1. Target data (all targets)
# 2. Drug data (all drugs)
# 3. Target-drug mapping
# 4. Target-disease mapping
# 5. Pathway data
# 6. Cross-references

# Download each file as needed
```

### Method 2: Individual Data Categories

```bash
# Target Information
# - P1-01-TTD_target_download.txt (basic target info)
# - P1-02-TTD_drug_disease.txt (drug-disease links)
# - P1-03-TTD_target_pathway.txt (pathway associations)

# Drug Information
# - P2-01-TTD_drug_download.txt (basic drug info)
# - P2-02-TTD_drug_synonyms.txt (drug name synonyms)

# Cross-references
# - P3-01-TTD_crossmatching.txt (external IDs)
# - P3-02-TTD_UniProt_mapping.txt (UniProt links)
```

### Method 3: Web Interface Export

```bash
# 1. Go to: https://idrblab.net/ttd/
# 2. Search by target name, gene, or drug
# 3. Export individual record data
# 4. Use for targeted queries
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| TTD_target_download.txt | ~5 MB | All targets |
| TTD_drug_download.txt | ~20 MB | All drugs |
| TTD_target_drug.txt | ~15 MB | Target-drug mapping |
| TTD_drug_disease.txt | ~10 MB | Drug-disease links |
| TTD_crossmatching.txt | ~5 MB | Cross-references |
| TTD_UniProt_mapping.txt | ~1 MB | UniProt IDs |

## File Format

TTD uses tab-separated format with header rows:

```
# Target download format:
TARGETID	FORESSION	TARGNAME	TARGTYPE	BIOCLASS	...
TTDT00001	P00533	Epidermal growth factor receptor	Successful target	Kinase	...

# Drug download format:
DRUGID	DRUGNAME	DRUGTPYE	HIGHEST_STATUS	COMPANY	...
D0A9YA	Erlotinib	Small molecule drug	Approved	...
```

## Post-Download Processing

```bash
# Preview file structure
head -10 TTD_target_download.txt

# Count targets by status
awk -F'\t' 'NR>1 {print $4}' TTD_target_download.txt | sort | uniq -c

# Extract successful targets
head -1 TTD_target_download.txt > successful_targets.txt
grep "Successful" TTD_target_download.txt >> successful_targets.txt

# Count drugs by status
awk -F'\t' 'NR>1 {print $4}' TTD_drug_download.txt | sort | uniq -c

# Extract approved drugs
head -1 TTD_drug_download.txt > approved_drugs.txt
grep "Approved" TTD_drug_download.txt >> approved_drugs.txt

# Load into SQLite
sqlite3 ttd.db << 'EOF'
CREATE TABLE targets (
  ttd_id TEXT PRIMARY KEY,
  uniprot_id TEXT,
  name TEXT,
  target_type TEXT,
  biochem_class TEXT
);

CREATE TABLE drugs (
  drug_id TEXT PRIMARY KEY,
  name TEXT,
  drug_type TEXT,
  highest_status TEXT,
  company TEXT
);

CREATE TABLE target_drug (
  target_id TEXT,
  drug_id TEXT,
  activity TEXT,
  moa TEXT,
  PRIMARY KEY (target_id, drug_id)
);

.mode tabs
.import TTD_target_download.txt targets
.import TTD_drug_download.txt drugs
EOF

# Query for kinase targets
sqlite3 ttd.db "SELECT ttd_id, name FROM targets WHERE biochem_class LIKE '%Kinase%';"

# Find drugs for specific target
sqlite3 ttd.db "SELECT d.name, d.highest_status FROM drugs d JOIN target_drug td ON d.drug_id = td.drug_id WHERE td.target_id = 'TTDT00001';"
```

## Data Analysis Examples

```bash
# Count drugs by development status
awk -F'\t' 'NR>1 {status[$4]++} END {for (s in status) print s, status[s]}' TTD_drug_download.txt

# Find targets with most approved drugs
awk -F'\t' 'NR>1 {count[$1]++} END {for (t in count) print t, count[t]}' TTD_target_drug.txt | sort -k2 -rn | head -20

# Cross-reference with UniProt
join -t$'\t' -1 1 -2 1 \
  <(sort TTD_target_download.txt) \
  <(sort TTD_UniProt_mapping.txt) > targets_with_uniprot.txt
```

## Cross-Reference Mapping

TTD provides mappings to:

| Database | Field |
|----------|-------|
| UniProt | Target protein |
| ChEMBL | Drug compounds |
| DrugBank | Drug information |
| PubChem | Chemical structures |
| KEGG | Pathways |
| ICD-11 | Diseases |

## Verification

```bash
# Check file integrity
file TTD_target_download.txt
# Should be: ASCII text

# Count columns
head -1 TTD_target_download.txt | tr '\t' '\n' | wc -l

# Verify target count
wc -l TTD_target_download.txt

# Verify drug count
wc -l TTD_drug_download.txt

# Check for expected TTD IDs
grep -c "TTDT" TTD_target_download.txt
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | Semi-annual |
| Web interface | Continuous |
| Last major update | 2024 |

## Common Issues

- **Text encoding**: Files are UTF-8; handle special characters
- **Tab delimiters**: Ensure proper tab handling in parsers
- **Header rows**: First row contains column names
- **Missing values**: Some fields may be empty
- **ID formats**: TTDT+5 digits for targets, 6-char for drugs

## Integration Notes

Connect TTD to other resources:

```bash
# UniProt link
# Use TTD UniProt mapping to connect to protein data

# ChEMBL link
# Drug ChEMBL IDs provided for compound data

# DrugBank link
# DrugBank IDs for comprehensive drug info

# ICD-11 link
# Disease codes for medical ontology
```

## License

Free for academic use. Contact for commercial licensing.
