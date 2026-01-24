---
id: download-swiss-model
title: "SWISS-MODEL Repository Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# SWISS-MODEL Repository Download Instructions

## Quick Start

```bash
# Download model for a specific UniProt accession
curl -o p53_model.pdb "https://swissmodel.expasy.org/repository/uniprot/P04637.pdb"
```

## Prerequisites

- **curl** or **wget** for downloads
- **Python** with requests library for batch downloads
- Web browser for interactive modeling
- 1-100GB storage depending on scope

## No Registration Required

Pre-computed models are freely available. Server modeling requires optional account.

## Download Methods

### Method 1: Single Protein Model Download

```bash
# Download best model for UniProt accession (PDB format)
curl -o model.pdb "https://swissmodel.expasy.org/repository/uniprot/P04637.pdb"

# Download specific model by model ID
curl -o model.pdb "https://swissmodel.expasy.org/repository/uniprot/P04637/1.pdb"

# Download in mmCIF format
curl -o model.cif "https://swissmodel.expasy.org/repository/uniprot/P04637/1.cif"
```

### Method 2: Repository API

```bash
# Get all models for a UniProt entry (JSON)
curl "https://swissmodel.expasy.org/repository/uniprot/P04637.json" > p53_models.json

# Parse and download all models
python3 << 'EOF'
import requests
import json

# Get model list
response = requests.get("https://swissmodel.expasy.org/repository/uniprot/P04637.json")
data = response.json()

# Download each model
for i, model in enumerate(data['result']['structures'], 1):
    url = model['coordinates']
    pdb_response = requests.get(url)
    with open(f"p53_model_{i}.pdb", 'w') as f:
        f.write(pdb_response.text)
    print(f"Downloaded model {i}: template {model['template']}, coverage {model['coverage']:.2f}")
EOF
```

### Method 3: Batch Download by UniProt List

```bash
# Create list of UniProt accessions
echo -e "P04637\nP00533\nP42345" > uniprot_list.txt

# Download models for each
while read acc; do
  echo "Downloading $acc..."
  curl -s -o "${acc}_model.pdb" \
    "https://swissmodel.expasy.org/repository/uniprot/${acc}.pdb"
  sleep 1  # Be polite to the server
done < uniprot_list.txt
```

### Method 4: Python Batch Download

```python
import requests
import time
import os

def download_swiss_model(uniprot_id, output_dir="."):
    """Download all SWISS-MODEL models for a UniProt ID."""

    # Get model metadata
    url = f"https://swissmodel.expasy.org/repository/uniprot/{uniprot_id}.json"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"No models found for {uniprot_id}")
        return []

    data = response.json()
    models = data.get('result', {}).get('structures', [])

    downloaded = []
    for i, model in enumerate(models, 1):
        coord_url = model.get('coordinates')
        if coord_url:
            pdb_response = requests.get(coord_url)
            if pdb_response.status_code == 200:
                filename = f"{uniprot_id}_model_{i}.pdb"
                filepath = os.path.join(output_dir, filename)
                with open(filepath, 'w') as f:
                    f.write(pdb_response.text)
                downloaded.append({
                    'file': filepath,
                    'template': model.get('template'),
                    'qmean': model.get('qmean'),
                    'coverage': model.get('coverage')
                })

    return downloaded

# Example usage
uniprot_ids = ['P04637', 'P00533', 'P42345']
for uid in uniprot_ids:
    models = download_swiss_model(uid, output_dir='swiss_models')
    print(f"{uid}: Downloaded {len(models)} models")
    time.sleep(1)
```

### Method 5: Using the Modeling Server (Custom Models)

```bash
# Submit sequence for modeling
curl -X POST "https://swissmodel.expasy.org/automodel" \
  -H "Content-Type: application/json" \
  -d '{
    "target_sequences": ">my_protein\nMEEPQSDPSVEPPLSQETFSDLWKLL...",
    "project_name": "my_project"
  }' > job_response.json

# Get job ID from response
JOB_ID=$(cat job_response.json | jq -r '.project_id')

# Check job status (poll until complete)
curl "https://swissmodel.expasy.org/project/${JOB_ID}/models/summary/"

# Download completed models
curl -o model.pdb "https://swissmodel.expasy.org/project/${JOB_ID}/models/01/model.pdb"
```

### Method 6: Interactive Server Access

1. Visit https://swissmodel.expasy.org
2. Paste sequence or upload FASTA
3. Configure modeling parameters
4. Submit and monitor progress
5. Download results via web interface

## File Inventory

### Per-Model Files

| File | Size | Description |
|------|------|-------------|
| {acc}_model_{n}.pdb | 50-500 KB | Model coordinates |
| {acc}_model_{n}.cif | 100-800 KB | mmCIF format |

### Repository API Response

| Field | Description |
|-------|-------------|
| structures | Array of available models |
| coverage_image | Coverage diagram URL |
| sequence_length | Target sequence length |

### Server Project Files

| File | Description |
|------|-------------|
| model.pdb | Final model coordinates |
| alignment.txt | Target-template alignment |
| report.pdf | Modeling report |
| qmean_local.txt | Per-residue quality |

## Post-Download Processing

```bash
# Extract QMEAN from B-factor column
python3 << 'EOF'
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure("model", "p53_model.pdb")

qmean_scores = []
for model in structure:
    for chain in model:
        for residue in chain:
            ca = residue['CA'] if 'CA' in residue else None
            if ca:
                qmean_scores.append(ca.bfactor)

print(f"Mean local QMEAN: {np.mean(qmean_scores):.3f}")
print(f"Min: {np.min(qmean_scores):.3f}, Max: {np.max(qmean_scores):.3f}")
EOF

# Filter models by quality
python3 << 'EOF'
import json
import requests

response = requests.get("https://swissmodel.expasy.org/repository/uniprot/P04637.json")
data = response.json()

# Keep only high-quality models (QMEAN > -2)
good_models = [m for m in data['result']['structures'] if m.get('qmean', -10) > -2]
print(f"Good models: {len(good_models)} / {len(data['result']['structures'])}")

for m in good_models:
    print(f"  Template: {m['template']}, QMEAN: {m['qmean']:.2f}, Coverage: {m['coverage']:.2f}")
EOF

# Combine multiple models for full coverage
# (Requires manual alignment if regions don't overlap)
cat p53_model_1.pdb > p53_combined.pdb

# Renumber residues
python3 << 'EOF'
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure("model", "p53_model.pdb")

# Renumber residues starting from 1
new_residue_num = 1
for model in structure:
    for chain in model:
        for residue in chain:
            residue.id = (' ', new_residue_num, ' ')
            new_residue_num += 1

io = PDBIO()
io.set_structure(structure)
io.save("p53_renumbered.pdb")
EOF
```

## Verification

```bash
# Check PDB format
head -50 p53_model.pdb

# Verify QMEAN in header
grep "QMEAN" p53_model.pdb

# Count atoms
grep "^ATOM" p53_model.pdb | wc -l

# Check model coverage
grep "REMARK" p53_model.pdb | grep -i "coverage"

# Validate structure geometry
# Using MolProbity or similar tool
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New models | Continuous |
| Template updates | With PDB releases (weekly) |
| Repository refresh | Monthly |

## Common Issues

- **No models found**: Sequence may lack homologous templates
- **Low coverage**: Multiple models may cover different regions
- **Low QMEAN**: Low sequence identity to template
- **Server queue**: Complex modeling may take hours
- **Rate limiting**: Add delays between requests

## Quality Filtering

```python
# Filter by quality criteria
def filter_models(models, min_qmean=-2, min_coverage=0.3, min_identity=0.3):
    """Filter models by quality metrics."""
    return [
        m for m in models
        if m.get('qmean', -10) >= min_qmean
        and m.get('coverage', 0) >= min_coverage
        and m.get('identity', 0) >= min_identity
    ]
```

## Comparison with Other Sources

| Aspect | SWISS-MODEL | AlphaFold | PDB |
|--------|-------------|-----------|-----|
| Method | Homology | AI prediction | Experimental |
| Coverage | Template-dependent | ~100% | Variable |
| Quality | QMEAN | pLDDT | Resolution |
| Best for | Known folds | Novel folds | Ground truth |

## Related Resources

- [PDB](../pdb/) - Template structures
- [AlphaFold DB](../alphafold.db/) - AI predictions
- [UniProt](../../7.1.protein.sequences.annotations/uniprot/) - Source sequences
