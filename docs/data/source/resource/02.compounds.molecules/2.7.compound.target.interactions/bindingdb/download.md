---
id: download-bindingdb
title: "BindingDB Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# BindingDB Download Instructions

## Quick Start

```bash
# Download complete BindingDB dataset
wget https://www.bindingdb.org/bind/downloads/BindingDB_All_202401.tsv.zip
unzip BindingDB_All_202401.tsv.zip
```

## Prerequisites

- **wget** or **curl** for downloads
- ~10-20 GB disk space for full dataset
- TSV/CSV parser or database software
- SDF viewer for structure files

## No Registration Required

Data is CC BY 3.0 licensed and freely downloadable.

## Download Methods

### Method 1: Complete Database (Recommended)

```bash
# Download full TSV (all data)
wget https://www.bindingdb.org/bind/downloads/BindingDB_All_202401.tsv.zip
unzip BindingDB_All_202401.tsv.zip

# Download SDF structures
wget https://www.bindingdb.org/bind/downloads/BindingDB_All_202401.sdf.zip
unzip BindingDB_All_202401.sdf.zip
```

### Method 2: Target-Specific Subsets

```bash
# Navigate to: https://www.bindingdb.org/download.jsp
# Available subsets:
# - By organism (Human, Mouse, etc.)
# - By target family (Kinases, GPCRs, etc.)
# - By data source (Patents, Publications)

# Example: Download human protein targets only
wget https://www.bindingdb.org/bind/downloads/BindingDB_Human_202401.tsv.zip
```

### Method 3: REST API

```bash
# Get ligands binding to UniProt target
curl "https://www.bindingdb.org/axis2/services/BDBService/getLigandsByUniprots?uniprot=P35354" | \
  xmllint --format -

# Get targets for a compound (by InChI Key)
curl "https://www.bindingdb.org/axis2/services/BDBService/getTargetByCompound?inchi_key=BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

# Get compound by BindingDB MonomerID
curl "https://www.bindingdb.org/axis2/services/BDBService/getLigandByMonomerID?monomerid=12345"
```

### Method 4: Web Interface Export

```bash
# 1. Go to: https://www.bindingdb.org
# 2. Search by target, compound, or structure
# 3. Select results and export to TSV
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| BindingDB_All_*.tsv.zip | ~3 GB | Complete data (TSV) |
| BindingDB_All_*.sdf.zip | ~5 GB | Complete structures (SDF) |
| BindingDB_Human_*.tsv.zip | ~2 GB | Human targets only |
| BindingDB_Patent_*.tsv.zip | ~1.5 GB | Patent-derived data |

## TSV Column Reference

| Column | Description |
|--------|-------------|
| BindingDB Reactant_set_id | Measurement ID |
| Ligand SMILES | Structure |
| Ligand InChI | InChI identifier |
| Ligand InChI Key | Structure hash |
| BindingDB MonomerID | Compound ID |
| Target Name | Protein name |
| Target Source Organism | Species |
| Ki (nM) | Binding constant |
| IC50 (nM) | Inhibition concentration |
| Kd (nM) | Dissociation constant |
| EC50 (nM) | Effective concentration |
| pH | Assay pH |
| Temp (C) | Temperature |
| Curation/DataSource | Source type |
| Article DOI | Reference |
| PubMed ID | PubMed reference |
| Patent Number | Patent reference |
| UniProt ID | Target UniProt |
| PDB IDs | Structure IDs |

## Post-Download Processing

```bash
# Preview TSV structure
head -5 BindingDB_All.tsv

# Count total records
wc -l BindingDB_All.tsv

# Extract unique targets
cut -f13 BindingDB_All.tsv | sort | uniq | wc -l

# Filter for specific target (e.g., kinases)
head -1 BindingDB_All.tsv > kinase_data.tsv
grep -i "kinase" BindingDB_All.tsv >> kinase_data.tsv

# Filter high-affinity binders (Ki < 100 nM)
awk -F'\t' 'NR==1 || ($7!="" && $7<100)' BindingDB_All.tsv > potent_binders.tsv

# Load into SQLite
sqlite3 bindingdb.db << 'EOF'
CREATE TABLE binding_data (
  reactant_id TEXT,
  smiles TEXT,
  inchi TEXT,
  inchi_key TEXT,
  monomer_id TEXT,
  ligand_name TEXT,
  ki_nm REAL,
  ic50_nm REAL,
  kd_nm REAL,
  ec50_nm REAL,
  target_name TEXT,
  organism TEXT,
  uniprot_id TEXT,
  pmid TEXT,
  patent TEXT
);
.mode tabs
.import BindingDB_All.tsv binding_data
CREATE INDEX idx_target ON binding_data(target_name);
CREATE INDEX idx_uniprot ON binding_data(uniprot_id);
CREATE INDEX idx_inchi ON binding_data(inchi_key);
EOF

# Query for specific target
sqlite3 bindingdb.db "SELECT DISTINCT ligand_name, ki_nm FROM binding_data WHERE uniprot_id='P35354' AND ki_nm < 100 ORDER BY ki_nm;"
```

## SDF Processing

```bash
# Count compounds in SDF
grep -c '^\$\$\$\$' BindingDB_All.sdf

# Convert to SMILES (requires Open Babel)
obabel BindingDB_All.sdf -O bindingdb.smi

# Extract subset by property
obabel BindingDB_All.sdf -O subset.sdf --filter "MW<500"
```

## Verification

```bash
# Check TSV integrity
file BindingDB_All.tsv

# Verify column count
head -1 BindingDB_All.tsv | tr '\t' '\n' | wc -l

# Check for expected columns
head -1 BindingDB_All.tsv | grep -o "Ki\|IC50\|Kd\|EC50"
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| Jan 2025 | 2025-01 | ~10 GB | Current |
| Monthly releases | Weekly | Varies | Rolling |

### Version Notes

BindingDB current release contains:
- 2.9M+ binding data entries
- 1.3M+ compounds
- 9,500+ protein targets
- Links to PDB, UniProt, ChEMBL, PubChem
- Patent and publication-derived data

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.bindingdb.org/axis2/services/BDBService` |
| Rate Limit | 60 req/min |
| Auth Required | No |
| Documentation | https://www.bindingdb.org/bind/BindingDBRESTfulAPI.jsp |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Full database | Weekly |
| Download files | Monthly |
| API data | Real-time |

## API Rate Limits

| Access | Limit |
|--------|-------|
| REST API | 60 requests/minute |
| Bulk queries | Contact for access |

## Common Issues

- **Large files**: Use streaming parser for multi-GB files
- **Missing values**: Many affinity columns may be empty
- **Multiple entries**: Same compound-target pair may have multiple measurements
- **Unit variations**: All affinities normalized to nM
- **TSV parsing**: Handle tabs in text fields properly

## Data Filtering Tips

```bash
# Keep only measurements with both IC50 and Ki
awk -F'\t' 'NR==1 || ($7!="" && $8!="")' BindingDB_All.tsv > dual_affinity.tsv

# Extract patent-only data
head -1 BindingDB_All.tsv > patent_data.tsv
awk -F'\t' '$NF!=""' BindingDB_All.tsv >> patent_data.tsv

# Get unique compound-target pairs
cut -f5,15 BindingDB_All.tsv | sort | uniq > compound_target_pairs.tsv
```

## License

CC BY 3.0 - Free for any use with attribution
