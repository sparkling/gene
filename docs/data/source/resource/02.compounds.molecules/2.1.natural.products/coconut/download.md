---
id: download-coconut
title: "COCONUT Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# COCONUT Download Instructions

## Quick Start

```bash
# Download complete database in SDF format
wget https://coconut.naturalproducts.net/download/coconut_db.sdf.gz
gunzip coconut_db.sdf.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- ~2 GB disk space for full dataset
- SDF viewer or chemistry toolkit (Open Babel, RDKit)
- JSON parser for metadata format

## No Registration Required

Data is CC BY 4.0 licensed and freely downloadable.

## Download Methods

### Method 1: Direct Download (Recommended)

```bash
# Complete database in SDF format (structures + metadata)
wget https://coconut.naturalproducts.net/download/coconut_db.sdf.gz
gunzip coconut_db.sdf.gz

# CSV format (tabular data)
wget https://coconut.naturalproducts.net/download/coconut_db.csv.gz
gunzip coconut_db.csv.gz

# JSON format (complete metadata)
wget https://coconut.naturalproducts.net/download/coconut_db.json.gz
gunzip coconut_db.json.gz

# SMILES format (structures only)
wget https://coconut.naturalproducts.net/download/coconut_structures.smiles.gz
gunzip coconut_structures.smiles.gz
```

### Method 2: REST API (Targeted Queries)

```bash
# Search by compound name
curl "https://coconut.naturalproducts.net/api/search?query=artemisinin" | jq

# Get compound by COCONUT ID
curl "https://coconut.naturalproducts.net/api/compound/CNP0123456" | jq

# Property-based search
curl "https://coconut.naturalproducts.net/api/search/properties?mw_min=200&mw_max=400&alogp_max=3" | jq

# Structure search (substructure)
curl -X POST "https://coconut.naturalproducts.net/api/structure/search" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "c1ccc2c(c1)cccc2", "searchType": "substructure"}' | jq

# Similarity search
curl -X POST "https://coconut.naturalproducts.net/api/structure/search" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(C)CCCC(C)CCCC", "searchType": "similarity", "threshold": 0.8}' | jq
```

### Method 3: GitHub Repository

```bash
# Clone the COCONUT repository for archived versions
git clone https://github.com/Steinbeck-Lab/coconut.git

# Check releases for specific versions
# https://github.com/Steinbeck-Lab/coconut/releases
```

### Method 4: Batch API Download

```bash
# Batch retrieval by IDs
curl -X POST "https://coconut.naturalproducts.net/api/batch/compounds" \
  -H "Content-Type: application/json" \
  -d '{"ids": ["CNP0123456", "CNP0234567", "CNP0345678"]}' | jq

# Batch download specific format
curl -X POST "https://coconut.naturalproducts.net/api/download" \
  -H "Content-Type: application/json" \
  -d '{"ids": ["CNP0123456"], "format": "sdf"}' > compounds.sdf
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| coconut_db.sdf.gz | ~500 MB | Structures + metadata |
| coconut_db.csv.gz | ~300 MB | Tabular properties |
| coconut_db.json.gz | ~600 MB | Complete JSON export |
| coconut_structures.smiles.gz | ~200 MB | SMILES only |

## SDF Properties

| Property | Description |
|----------|-------------|
| COCONUT_ID | CNP identifier |
| NAME | Compound name |
| MOLECULAR_FORMULA | Formula |
| MOLECULAR_WEIGHT | MW in Daltons |
| INCHI | InChI string |
| INCHI_KEY | InChI Key |
| SMILES | Canonical SMILES |
| QED | Drug-likeness score |
| ALOGP | Lipophilicity |
| SOURCES | Origin databases |

## Post-Download Processing

```bash
# Count compounds in SDF
grep -c '^\$\$\$\$' coconut_db.sdf

# Extract COCONUT IDs
grep "COCONUT_ID" coconut_db.sdf | head -20

# Convert to other formats (requires Open Babel)
obabel coconut_db.sdf -O coconut_db.mol2
obabel coconut_db.sdf -O coconut_db.inchi

# Extract structures as SMILES
obabel coconut_db.sdf -O coconut.smi

# Process JSON format
zcat coconut_db.json.gz | jq 'length'
zcat coconut_db.json.gz | jq '.[0]'

# Filter drug-like compounds (QED > 0.5)
zcat coconut_db.json.gz | jq '[.[] | select(.qed_drug_likeliness > 0.5)]' > druglike.json

# Extract compounds from specific source
zcat coconut_db.json.gz | jq '[.[] | select(.sources | contains(["ZINC"]))]' > zinc_compounds.json

# Load CSV into database
gunzip -c coconut_db.csv.gz > coconut.csv
sqlite3 coconut.db << 'EOF'
.mode csv
.import coconut.csv compounds
CREATE INDEX idx_coconut_id ON compounds(coconut_id);
CREATE INDEX idx_inchi_key ON compounds(inchi_key);
EOF

# Query by properties
sqlite3 coconut.db "SELECT coconut_id, name, molecular_weight, qed_drug_likeliness FROM compounds WHERE molecular_weight < 500 AND qed_drug_likeliness > 0.7 LIMIT 20;"
```

## RDKit Processing (Python)

```python
from rdkit import Chem
from rdkit.Chem import Descriptors
import gzip

# Read SDF file
with gzip.open('coconut_db.sdf.gz', 'rt') as f:
    suppl = Chem.ForwardSDMolSupplier(f)
    for mol in suppl:
        if mol is not None:
            coconut_id = mol.GetProp('COCONUT_ID')
            mw = Descriptors.MolWt(mol)
            print(f"{coconut_id}: MW = {mw:.2f}")

# Filter by substructure
pattern = Chem.MolFromSmarts('c1ccc2c(c1)cccc2')  # naphthalene
matches = [mol for mol in suppl if mol and mol.HasSubstructMatch(pattern)]
```

## Verification

```bash
# Check SDF integrity
grep -c '^\$\$\$\$' coconut_db.sdf

# Verify JSON format
zcat coconut_db.json.gz | jq 'length'

# Check CSV structure
head -5 coconut.csv

# Validate sample entries
zcat coconut_db.json.gz | jq '.[0:5] | .[].coconut_id'
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| COCONUT 2.0 | 2025-01-06 | ~1.5 GB | Current |
| COCONUT 1.x | 2021-2024 | ~1 GB | Archived |

### Version Notes

COCONUT 2.0 represents a comprehensive overhaul and curation of the database:
- Aggregates data from 63+ open natural product sources
- Uses ChEMBL curation pipeline with RDKit post-processing
- Standardizes molecular structures and metadata while preserving stereochemistry

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://coconut.naturalproducts.net/api` |
| Rate Limit | 60-120 req/min (see below) |
| Auth Required | No |
| Documentation | https://coconut.naturalproducts.net/api/docs |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Full database | Quarterly |
| Incremental updates | Monthly |
| API data | Real-time |

## API Rate Limits

| Access | Limit |
|--------|-------|
| Search queries | 60/minute |
| Compound retrieval | 120/minute |
| Batch operations | 10/minute |

## Common Issues

- **Large SDF files**: Use streaming readers (RDKit, Open Babel)
- **Memory usage**: Process JSON in chunks for large datasets
- **Structure validation**: Some compounds may have non-standard features
- **Missing properties**: Not all compounds have complete metadata

## Statistics Endpoint

```bash
# Get database statistics
curl "https://coconut.naturalproducts.net/api/stats" | jq

# Response includes:
# - total_compounds
# - total_organisms
# - total_sources
# - molecular_weight_distribution
```

## License

CC BY 4.0 - Free for any use with attribution
