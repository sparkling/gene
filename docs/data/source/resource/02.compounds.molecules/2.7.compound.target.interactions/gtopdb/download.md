---
id: download-gtopdb
title: "GtoPdb Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# GtoPdb Download Instructions

## Quick Start

```bash
# Download complete database as PostgreSQL dump
wget https://www.guidetopharmacology.org/DATA/public_iuphardb.zip
unzip public_iuphardb.zip
```

## Prerequisites

- **wget** or **curl** for downloads
- **PostgreSQL** for database restore
- ~2 GB disk space
- Alternatively: CSV tools for flat file access

## No Registration Required

Data is CC BY-SA 4.0 licensed and freely downloadable.

## Download Methods

### Method 1: PostgreSQL Database (Recommended)

```bash
# Download PostgreSQL dump
wget https://www.guidetopharmacology.org/DATA/public_iuphardb.zip
unzip public_iuphardb.zip

# Create database and restore
createdb gtopdb
pg_restore -d gtopdb public_iuphardb_*.backup

# Or with psql for SQL format
psql -d gtopdb -f public_iuphardb_*.sql
```

### Method 2: CSV Downloads

```bash
# Download page: https://www.guidetopharmacology.org/download.jsp

# Targets
wget https://www.guidetopharmacology.org/DATA/targets_and_families.csv

# Ligands
wget https://www.guidetopharmacology.org/DATA/ligands.csv

# Interactions
wget https://www.guidetopharmacology.org/DATA/interactions.csv

# Approved drugs
wget https://www.guidetopharmacology.org/DATA/approved_drugs.csv
```

### Method 3: REST API

```bash
# List all targets
curl "https://www.guidetopharmacology.org/services/targets" | jq

# Get specific target
curl "https://www.guidetopharmacology.org/services/targets/290" | jq

# Get target ligands
curl "https://www.guidetopharmacology.org/services/targets/290/ligands" | jq

# List all ligands
curl "https://www.guidetopharmacology.org/services/ligands" | jq

# Get specific ligand
curl "https://www.guidetopharmacology.org/services/ligands/5085" | jq

# Search ligands
curl "https://www.guidetopharmacology.org/services/ligands?name=vemurafenib" | jq
```

### Method 4: Web Interface Export

```bash
# 1. Navigate to: https://www.guidetopharmacology.org
# 2. Search for target or ligand
# 3. Use export options on results page
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| public_iuphardb.zip | ~200 MB | PostgreSQL database |
| targets_and_families.csv | ~2 MB | Target information |
| ligands.csv | ~10 MB | Ligand data |
| interactions.csv | ~30 MB | Target-ligand interactions |
| approved_drugs.csv | ~1 MB | Approved drug subset |

## CSV Column Reference

### targets_and_families.csv

| Column | Description |
|--------|-------------|
| Target id | GtoPdb identifier |
| Target name | Full name |
| Target abbreviated name | Short name |
| Target type | receptor, enzyme, etc. |
| Family name | Target family |
| HGNC symbol | Gene symbol |
| UniProt ID | UniProt accession |
| Human SwissProt | Swiss-Prot entry |

### interactions.csv

| Column | Description |
|--------|-------------|
| Target id | Target identifier |
| Target name | Target name |
| Ligand id | Ligand identifier |
| Ligand name | Ligand name |
| Type | Interaction type |
| Action | Activation/Inhibition |
| Affinity parameter | pKi, pIC50, etc. |
| Affinity median | Median pX value |
| Affinity low | Low pX |
| Affinity high | High pX |
| Species | Assay species |
| Primary target | Primary/Secondary |

## Post-Download Processing

```bash
# Preview CSV data
head -5 targets_and_families.csv
head -5 interactions.csv

# Count targets by type
cut -d',' -f4 targets_and_families.csv | sort | uniq -c

# Count ligands
wc -l ligands.csv

# Filter for approved drugs with interactions
head -1 interactions.csv > approved_interactions.csv
grep -i "approved" ligands.csv | cut -d',' -f1 > approved_ids.txt
grep -f approved_ids.txt interactions.csv >> approved_interactions.csv

# Convert pKi to Ki (nM)
awk -F',' 'NR>1 && $8!="" {
  pki = $8;
  ki_nm = 10^(9-pki);
  print $0","ki_nm
}' interactions.csv > interactions_with_ki.csv

# Load into SQLite (alternative to PostgreSQL)
sqlite3 gtopdb.db << 'EOF'
.mode csv
.import targets_and_families.csv targets
.import ligands.csv ligands
.import interactions.csv interactions
CREATE INDEX idx_target ON interactions(target_id);
CREATE INDEX idx_ligand ON interactions(ligand_id);
EOF

# Query for kinase inhibitors
sqlite3 gtopdb.db "SELECT l.name, i.affinity_median FROM interactions i JOIN ligands l ON i.ligand_id = l.ligand_id WHERE i.type='inhibitor' AND i.affinity_median > 7 ORDER BY i.affinity_median DESC LIMIT 20;"
```

## PostgreSQL Queries

```sql
-- After restoring database
-- Get all GPCR targets
SELECT target_id, name, abbreviation
FROM targets
WHERE type = 'GPCR';

-- Get top affinity interactions
SELECT t.name as target, l.name as ligand,
       i.affinity_type, i.affinity_median
FROM interactions i
JOIN targets t ON i.target_id = t.target_id
JOIN ligands l ON i.ligand_id = l.ligand_id
WHERE i.affinity_median > 9
ORDER BY i.affinity_median DESC
LIMIT 50;

-- Find selective ligands
SELECT l.name, COUNT(DISTINCT i.target_id) as target_count
FROM ligands l
JOIN interactions i ON l.ligand_id = i.ligand_id
WHERE i.affinity_median > 7
GROUP BY l.ligand_id, l.name
HAVING COUNT(DISTINCT i.target_id) = 1;
```

## Verification

```bash
# Check CSV integrity
file targets_and_families.csv
# Should be: ASCII text

# Verify column count
head -1 targets_and_families.csv | tr ',' '\n' | wc -l

# Check PostgreSQL restore
psql -d gtopdb -c "SELECT COUNT(*) FROM targets;"
psql -d gtopdb -c "SELECT COUNT(*) FROM ligands;"
psql -d gtopdb -c "SELECT COUNT(*) FROM interactions;"
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| 2024.3 | 2024-12 | ~300 MB | Current |
| 2024.2 | 2024-09 | ~290 MB | Archived |
| 2024.1 | 2024-03 | ~280 MB | Archived |

### Version Notes

GtoPdb 2024.3 (Guide to Pharmacology) contains:
- 3,000+ human drug targets
- 11,000+ ligands
- 20,000+ quantitative interactions
- Expert-curated pharmacology data
- IUPHAR/BPS nomenclature standards

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.guidetopharmacology.org/services` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://www.guidetopharmacology.org/webServices.jsp |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database releases | Quarterly |
| Web updates | Continuous |
| API data | Real-time |

## API Rate Limits

| Access | Limit |
|--------|-------|
| REST API | Reasonable use |
| Bulk queries | Use downloads |

## Common Issues

- **pX notation**: Remember to convert pKi to Ki if needed
- **Multiple species**: Filter for human data if required
- **Primary/secondary**: Primary targets indicate main pharmacological target
- **Selectivity**: Non-selective ligands have multiple target interactions

## License

CC BY-SA 4.0 - Free for any use with attribution and share-alike
