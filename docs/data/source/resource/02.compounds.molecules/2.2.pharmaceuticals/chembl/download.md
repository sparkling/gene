---
id: download-chembl
title: "ChEMBL Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# ChEMBL Download Instructions

## Quick Start

```bash
# Download latest ChEMBL SQLite database
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34_sqlite.tar.gz
tar -xzf chembl_34_sqlite.tar.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **SQLite** or **PostgreSQL/MySQL** for database access
- **tar** for extraction
- Approximately 20-50GB disk space for full database

## No Registration Required

ChEMBL data is freely available under CC BY-SA 3.0 license without registration.

## Download Methods

### Method 1: SQLite Database (Recommended for Local Use)

```bash
# Download latest SQLite version
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34_sqlite.tar.gz

# Extract
tar -xzf chembl_34_sqlite.tar.gz

# Query the database
sqlite3 chembl_34/chembl_34_sqlite/chembl_34.db \
  "SELECT COUNT(*) FROM compound_records;"
```

### Method 2: PostgreSQL Database Dump

```bash
# Download PostgreSQL dump
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34_postgresql.tar.gz

# Extract
tar -xzf chembl_34_postgresql.tar.gz

# Create database and restore
createdb chembl_34
pg_restore -d chembl_34 chembl_34_postgresql/chembl_34_postgresql.dmp
```

### Method 3: MySQL Database Dump

```bash
# Download MySQL dump
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34_mysql.tar.gz

# Extract
tar -xzf chembl_34_mysql.tar.gz

# Create database and import
mysql -e "CREATE DATABASE chembl_34;"
mysql chembl_34 < chembl_34_mysql/chembl_34_mysql.dmp
```

### Method 4: Flat Files (SDF, TSV)

```bash
# Download compound structures (SDF)
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34.sdf.gz
gunzip chembl_34.sdf.gz

# Download UniChem mapping
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34_unichem_mapping.txt.gz
gunzip chembl_34_unichem_mapping.txt.gz

# Download target dictionary
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34_chemreps.txt.gz
```

### Method 5: FTP Bulk Download

```bash
# List all available files
curl -s https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/

# Download via FTP
ftp ftp.ebi.ac.uk
# cd /pub/databases/chembl/ChEMBLdb/latest/
# mget *
```

## File Inventory

### Database Files

| File | Size | Description |
|------|------|-------------|
| chembl_34_sqlite.tar.gz | ~1.5 GB | SQLite database |
| chembl_34_postgresql.tar.gz | ~3 GB | PostgreSQL dump |
| chembl_34_mysql.tar.gz | ~2.5 GB | MySQL dump |

### Structure Files

| File | Size | Description |
|------|------|-------------|
| chembl_34.sdf.gz | ~2 GB | Compound structures (SDF) |
| chembl_34.fa.gz | ~50 MB | Target sequences (FASTA) |
| chembl_34_chemreps.txt.gz | ~500 MB | Chemical representations (SMILES, InChI) |

### Mapping Files

| File | Size | Description |
|------|------|-------------|
| chembl_34_unichem_mapping.txt.gz | ~100 MB | UniChem cross-references |
| chembl_34_uniprot_mapping.txt.gz | ~10 MB | UniProt target mapping |

### Subset Files

| File | Size | Description |
|------|------|-------------|
| chembl_34_monomer_library.xml | ~5 MB | Monomer definitions |
| chembl_34_bio_component.xml | ~50 MB | Biological components |

## Post-Download Processing

```bash
# Index SQLite for better performance
sqlite3 chembl_34.db << 'EOF'
CREATE INDEX IF NOT EXISTS idx_activities_molregno ON activities(molregno);
CREATE INDEX IF NOT EXISTS idx_activities_assay_id ON activities(assay_id);
CREATE INDEX IF NOT EXISTS idx_compound_structures_molregno ON compound_structures(molregno);
ANALYZE;
EOF

# Extract bioactivity data
sqlite3 -header -separator $'\t' chembl_34.db << 'EOF' > bioactivities.tsv
SELECT
    m.chembl_id AS compound_chembl_id,
    t.chembl_id AS target_chembl_id,
    a.standard_type,
    a.standard_value,
    a.standard_units
FROM activities a
JOIN molecule_dictionary m ON a.molregno = m.molregno
JOIN target_dictionary t ON a.target_id = t.tid
WHERE a.standard_type IN ('IC50', 'EC50', 'Ki', 'Kd')
AND a.standard_value IS NOT NULL;
EOF

# Convert SDF to SMILES
obabel chembl_34.sdf -O chembl_34.smi
```

## Verification

```bash
# Verify archive integrity
tar -tzf chembl_34_sqlite.tar.gz | head

# Check database statistics
sqlite3 chembl_34.db << 'EOF'
SELECT 'Compounds:', COUNT(*) FROM molecule_dictionary;
SELECT 'Assays:', COUNT(*) FROM assays;
SELECT 'Activities:', COUNT(*) FROM activities;
SELECT 'Targets:', COUNT(*) FROM target_dictionary;
EOF

# Verify SDF file
grep -c '\$\$\$\$' chembl_34.sdf
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| ChEMBL 34 | 2024-06 | ~7 GB total | Current |
| ChEMBL 33 | 2023-12 | ~6.5 GB | Archived |
| ChEMBL 32 | 2023-06 | ~6 GB | Archived |

### Version Notes

ChEMBL 34 (latest) contains:
- 2.4M+ compounds with drug-like properties
- 20M+ bioactivity measurements
- 15,000+ targets
- Natural Product likeness scores
- Updated Chemical Probes flag
- Patent bioactivity data

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.ebi.ac.uk/chembl/api/data` |
| Rate Limit | No hard limit (be reasonable) |
| Auth Required | No |
| Documentation | https://chembl.gitbook.io/chembl-interface-documentation |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major releases | Quarterly |
| Data updates | Monthly |
| Bug fixes | As needed |

## Common Issues

- **Large file downloads**: Use wget with `-c` flag to resume interrupted downloads
- **SQLite locks**: Use WAL mode for concurrent access: `PRAGMA journal_mode=WAL;`
- **Memory with large queries**: Use LIMIT and pagination for large result sets
- **SDF parsing errors**: Some records may have non-standard features; use tolerant parsers
- **Version changes**: Schema may change between versions; check release notes

## API Access (Alternative)

```bash
# ChEMBL REST API (no download required)
curl "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25.json"

# Search for compounds
curl "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=aspirin"

# Get bioactivities
curl "https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=CHEMBL25"
```

## Data Schema Reference

Key tables in ChEMBL database:

| Table | Description |
|-------|-------------|
| molecule_dictionary | Compound registry |
| compound_structures | Chemical structures (SMILES, InChI) |
| activities | Bioactivity measurements |
| assays | Assay definitions |
| target_dictionary | Biological targets |
| drug_mechanism | Drug mechanisms of action |
| metabolism | Metabolic transformations |

## License

CC BY-SA 3.0 - Free for any use with attribution
