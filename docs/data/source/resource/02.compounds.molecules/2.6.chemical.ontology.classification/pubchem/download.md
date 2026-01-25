---
id: download-pubchem
title: "PubChem Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# PubChem Download Instructions

## Quick Start

```bash
# Download compound SDF file (first chunk)
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_000000001_000025000.sdf.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **rsync** for bulk mirroring
- **RDKit** or **Open Babel** for structure processing
- 100GB-10TB storage depending on data scope

## No Registration Required

PubChem data is public domain (no restrictions).

## Download Methods

### Method 1: FTP Compound Downloads

```bash
# Download compound structures (SDF format, chunked)
wget -r -np -nH --cut-dirs=4 \
  https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/

# Download specific range
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_000000001_000025000.sdf.gz

# Download compound properties
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/Extras/CID-Property.gz

# Download SMILES
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/Extras/CID-SMILES.gz

# Download InChI
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/Extras/CID-InChI.gz
```

### Method 2: Substance Downloads

```bash
# Substance SDF files
wget -r -np -nH --cut-dirs=4 \
  https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/CURRENT-Full/SDF/

# Substance-to-compound mapping
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/Extras/CID-SID.gz
```

### Method 3: Bioassay Data

```bash
# Bioassay descriptions
wget -r -np -nH --cut-dirs=4 \
  https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/

# Bioassay data
wget -r -np -nH --cut-dirs=4 \
  https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/

# Specific assay
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/0000001_0001000.zip
```

### Method 4: PUG REST API

```bash
# Get compound by CID
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/JSON" \
  -o aspirin.json

# Get compound properties
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON" \
  -o aspirin_props.json

# Get by name
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/JSON" \
  -o aspirin_byname.json

# Get structure (SDF)
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF" \
  -o aspirin.sdf

# Get structure (SMILES)
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/CanonicalSMILES/TXT" \
  -o aspirin.smi
```

### Method 5: Batch API Downloads

```bash
# Multiple CIDs
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244,5988,3672/property/MolecularFormula,MolecularWeight/JSON" \
  -o batch_props.json

# POST for large lists
curl -X POST "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/MolecularFormula,MolecularWeight/JSON" \
  --data "cid=2244,5988,3672,2519,5090" \
  -o large_batch.json

# Substructure search (async)
curl -X POST "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/cids/JSON?SMILES=c1ccccc1" \
  -o benzene_substructure.json
```

### Method 6: Classification and Annotations

```bash
# Compound classification
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Annotation.gz

# MeSH annotations
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-MeSH

# Patent links
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Patent-Count.gz

# PubMed links
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-PMID.gz
```

## File Inventory

### Compound Data

| File Pattern | Size | Description |
|--------------|------|-------------|
| Compound_*_*.sdf.gz | ~10-50 MB each | Structure files |
| CID-SMILES.gz | ~5 GB | All SMILES |
| CID-InChI.gz | ~8 GB | All InChI |
| CID-Property.gz | ~3 GB | Properties |

### Substance Data

| File Pattern | Size | Description |
|--------------|------|-------------|
| Substance_*.sdf.gz | ~5-20 MB each | Substance structures |
| CID-SID.gz | ~5 GB | Compound-substance mapping |

### Bioassay Data

| File Type | Size | Description |
|-----------|------|-------------|
| Description/*.csv.gz | ~1 GB total | Assay descriptions |
| Data/*.zip | ~500 GB total | Assay results |

### Full Database

| Dataset | Total Size |
|---------|------------|
| Compounds (SDF) | ~1 TB |
| Substances (SDF) | ~2 TB |
| Bioassays | ~500 GB |

## Post-Download Processing

```bash
# Decompress
gunzip Compound_000000001_000025000.sdf.gz

# Parse SDF with RDKit
python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# Read SDF file
suppl = Chem.SDMolSupplier('Compound_000000001_000025000.sdf')

data = []
for mol in suppl:
    if mol is None:
        continue
    data.append({
        'cid': mol.GetProp('PUBCHEM_COMPOUND_CID'),
        'mw': Descriptors.MolWt(mol),
        'smiles': Chem.MolToSmiles(mol)
    })

df = pd.DataFrame(data)
df.to_csv('compounds.tsv', sep='\t', index=False)
print(f"Processed {len(df)} compounds")
EOF

# Parse CID-SMILES
zcat CID-SMILES.gz | head -1000 > sample_smiles.txt

# Convert to other formats
obabel Compound_000000001_000025000.sdf -O compounds.mol2 -m

# Create similarity database
python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import pickle

suppl = Chem.SDMolSupplier('Compound_000000001_000025000.sdf')
fps = {}

for mol in suppl:
    if mol is None:
        continue
    cid = mol.GetProp('PUBCHEM_COMPOUND_CID')
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    fps[cid] = fp

with open('fingerprints.pkl', 'wb') as f:
    pickle.dump(fps, f)
EOF
```

## Verification

```bash
# Check SDF format
head -100 Compound_000000001_000025000.sdf

# Count compounds in SDF
grep -c '^\$\$\$\$' Compound_000000001_000025000.sdf

# Verify SMILES file
zcat CID-SMILES.gz | wc -l
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| PubChem 2025 | 2025-01 (NAR paper) | Multi-TB | Current |
| Continuous updates | Daily | N/A | Ongoing |

### Version Notes

PubChem 2025 update statistics (as of Sept 2024):
- 119 million compounds
- 322 million substances
- 295 million bioactivities
- 1000+ data sources
- 130+ new sources added in past two years

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` |
| Rate Limit | 5 req/sec |
| Auth Required | No |
| Documentation | https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest |

## Update Schedule

| Data Type | Frequency |
|-----------|-----------|
| Compounds | Daily |
| Bioassays | Daily |
| Full refresh | Weekly |

## Common Issues

- **Large files**: Download in chunks; use rsync for resume
- **Rate limits**: API has limits; batch requests and add delays
- **Structure parsing**: Some compounds may fail; handle errors gracefully
- **Stereochemistry**: Not all compounds have defined stereochemistry
- **Duplicates**: Multiple SIDs may map to same CID

## API Limits

| Limit Type | Value |
|------------|-------|
| Requests/second | 5 |
| CIDs per request | 100 (POST: 10000) |
| Async list size | 100000 |
| Result timeout | 30 seconds |

## PUG REST Endpoints

| Endpoint | Description |
|----------|-------------|
| /compound/cid/{cid} | Get by CID |
| /compound/name/{name} | Get by name |
| /compound/smiles/{smiles} | Get by SMILES |
| /compound/inchi/{inchi} | Get by InChI |
| /compound/substructure | Substructure search |
| /compound/similarity | Similarity search |

## Related Resources

- [ChEBI](../chebi/) - Chemical ontology
- [ChEMBL](../../2.2.pharmaceuticals/chembl/) - Bioactivity data
- [DrugBank](../../2.2.pharmaceuticals/drugbank/) - Drug data
