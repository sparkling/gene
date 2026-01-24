---
id: download-swissadme
title: "SwissADME Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# SwissADME Download Instructions

## Quick Start

```bash
# SwissADME is a web-only tool - no bulk download available
# Access via web interface: http://www.swissadme.ch

# For batch processing, prepare SMILES file and use web interface
# or programmatic submission (see below)
```

## Prerequisites

- **Web browser** for interactive use
- **Python 3.x** with `requests` and `BeautifulSoup4` for programmatic access
- **SMILES strings** for input molecules
- No registration required

## No Registration Required

SwissADME is freely accessible without registration for all users including commercial use.

## Access Methods

### Method 1: Web Interface (Recommended)

```bash
# 1. Open browser to http://www.swissadme.ch
# 2. Enter SMILES (up to 100 molecules, one per line)
# 3. Submit and wait for processing
# 4. Download results as CSV
```

**Input Format:**
- One SMILES per line
- Optional: Include molecule name after SMILES (tab-separated)
- Maximum 100 molecules per submission

### Method 2: Programmatic Access (Python)

```python
#!/usr/bin/env python3
"""
SwissADME programmatic access
Note: Web scraping - respect usage limits
"""
import requests
from time import sleep

def query_swissadme(smiles_list, delay=2):
    """
    Query SwissADME for ADME predictions

    Args:
        smiles_list: List of SMILES strings (max 100)
        delay: Seconds between requests

    Returns:
        Response content
    """
    url = "http://www.swissadme.ch/index.php"

    # Format SMILES input
    smiles_input = "\n".join(smiles_list[:100])

    # Submit form
    data = {
        'smiles': smiles_input
    }

    response = requests.post(url, data=data)
    sleep(delay)  # Be respectful

    return response.text

# Example usage
smiles = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    "CC(=O)NC1=CC=C(C=C1)O",      # Paracetamol
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
]

result = query_swissadme(smiles)
print("Results received")
```

### Method 3: Batch Processing Script

```python
#!/usr/bin/env python3
"""
Batch SwissADME processor
Handles large compound sets in chunks
"""
import pandas as pd
import requests
from bs4 import BeautifulSoup
from time import sleep
import csv

def batch_swissadme(input_file, output_file, chunk_size=50, delay=5):
    """
    Process large SMILES file through SwissADME

    Args:
        input_file: CSV/TSV with SMILES column
        output_file: Output CSV path
        chunk_size: Molecules per request (max 100)
        delay: Seconds between requests
    """
    # Read input
    df = pd.read_csv(input_file, sep='\t')
    smiles_col = 'SMILES' if 'SMILES' in df.columns else df.columns[0]

    all_results = []

    # Process in chunks
    for i in range(0, len(df), chunk_size):
        chunk = df[smiles_col].iloc[i:i+chunk_size].tolist()
        print(f"Processing compounds {i+1}-{min(i+chunk_size, len(df))}")

        # Submit to SwissADME
        url = "http://www.swissadme.ch/index.php"
        smiles_input = "\n".join(chunk)

        response = requests.post(url, data={'smiles': smiles_input})

        # Parse results (simplified - actual parsing depends on output format)
        # Store raw HTML or parse CSV link
        all_results.append({
            'chunk': i // chunk_size,
            'count': len(chunk),
            'status': response.status_code
        })

        sleep(delay)

    # Save status
    pd.DataFrame(all_results).to_csv(output_file, index=False)
    print(f"Processed {len(df)} compounds in {len(all_results)} batches")

# Example
# batch_swissadme('compounds.tsv', 'swissadme_status.csv')
```

### Method 4: Alternative - Local ADME Tools

For bulk local processing, consider these alternatives:

```bash
# RDKit descriptors (open source)
pip install rdkit-pypi

# ADMETlab 2.0 (free web tool with API)
# https://admetlab3.scbdd.com/

# pkCSM (free web tool)
# https://biosig.lab.uq.edu.au/pkcsm/
```

```python
# Local descriptor calculation with RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def calculate_adme_descriptors(smiles):
    """Calculate basic ADME-relevant descriptors locally"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'RotatableBonds': Descriptors.NumRotatableBonds(mol),
        'Lipinski_Violations': sum([
            Descriptors.MolWt(mol) > 500,
            Descriptors.MolLogP(mol) > 5,
            Descriptors.NumHDonors(mol) > 5,
            Descriptors.NumHAcceptors(mol) > 10
        ])
    }

# Process file locally
import pandas as pd

df = pd.read_csv('compounds.tsv', sep='\t')
results = df['SMILES'].apply(calculate_adme_descriptors)
results_df = pd.DataFrame(results.tolist())
results_df.to_csv('adme_descriptors.csv', index=False)
```

## File Inventory

SwissADME is web-only; no downloadable database files exist.

| Resource | Type | Description |
|----------|------|-------------|
| Web Interface | Interactive | Main prediction tool |
| CSV Export | On-demand | Results for submitted molecules |
| BOILED-Egg Plot | PNG | Visual GI/BBB assessment |

## Input File Format

### SMILES Input (swissadme_input.txt)

```
CC(=O)OC1=CC=CC=C1C(=O)O	Aspirin
CC(=O)NC1=CC=C(C=C1)O	Paracetamol
CN1C=NC2=C1C(=O)N(C(=O)N2C)C	Caffeine
CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F	Celecoxib
```

### Batch Input CSV

```csv
name,smiles
Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Paracetamol,CC(=O)NC1=CC=C(C=C1)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

## Post-Processing

### Parse CSV Results

```python
import pandas as pd

# Load SwissADME CSV export
df = pd.read_csv('swissadme_results.csv')

# Filter drug-like compounds
drug_like = df[
    (df['Lipinski_Violations'] == 0) &
    (df['GI_Absorption'] == 'High') &
    (df['PAINS_Alerts'] == 0)
]

print(f"Drug-like compounds: {len(drug_like)}/{len(df)}")

# Export hits
drug_like.to_csv('drug_like_hits.csv', index=False)
```

### Analyze BOILED-Egg Results

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('swissadme_results.csv')

# Create BOILED-Egg plot
fig, ax = plt.subplots(figsize=(10, 8))

# Plot compounds
colors = []
for _, row in df.iterrows():
    if row['GI_Absorption'] == 'High' and row['BBB_Permeant'] == 'Yes':
        colors.append('yellow')  # Yolk - BBB permeable
    elif row['GI_Absorption'] == 'High':
        colors.append('white')   # White - GI absorbed
    else:
        colors.append('gray')    # Outside - poor absorption

ax.scatter(df['WLOGP'], df['TPSA'], c=colors, edgecolors='black', s=100)
ax.set_xlabel('WLOGP')
ax.set_ylabel('TPSA')
ax.set_title('BOILED-Egg Model')

# Add reference regions
from matplotlib.patches import Ellipse
white = Ellipse((2.0, 70), 9, 140, fill=False, color='gray', linestyle='--')
yolk = Ellipse((3.0, 40), 6, 80, fill=False, color='orange', linestyle='--')
ax.add_patch(white)
ax.add_patch(yolk)

plt.savefig('boiled_egg_results.png', dpi=150)
```

## Verification

```python
# Verify SwissADME results
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

df = pd.read_csv('swissadme_results.csv')

# Verify MW calculations
for _, row in df.head().iterrows():
    mol = Chem.MolFromSmiles(row['Canonical_SMILES'])
    if mol:
        calc_mw = Descriptors.MolWt(mol)
        reported_mw = row['MW']
        diff = abs(calc_mw - reported_mw)
        print(f"{row['Molecule']}: Reported={reported_mw:.2f}, Calculated={calc_mw:.2f}, Diff={diff:.2f}")
```

```bash
# Check CSV structure
head -5 swissadme_results.csv

# Count columns
head -1 swissadme_results.csv | tr ',' '\n' | wc -l
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Model Updates | Periodic |
| Web Interface | Continuous |
| Algorithm Improvements | As published |

## Rate Limits and Best Practices

| Consideration | Recommendation |
|---------------|----------------|
| Batch Size | Max 100 compounds per submission |
| Request Delay | 2-5 seconds between submissions |
| Large Sets | Split into multiple batches |
| Attribution | Cite SwissADME paper |

## Common Issues

- **Invalid SMILES**: Validate SMILES with RDKit before submission
- **Timeout**: Reduce batch size for complex molecules
- **Results parsing**: CSV export is most reliable format
- **No API**: Web scraping needed for automation; respect usage policies

## Integration with Other Tools

```bash
# Combine with ChEMBL data
python3 << 'EOF'
import pandas as pd

# Load ChEMBL compounds
chembl = pd.read_csv('chembl_compounds.tsv', sep='\t')

# Prepare for SwissADME
swissadme_input = chembl[['canonical_smiles', 'pref_name']].copy()
swissadme_input.columns = ['SMILES', 'Name']

# Save batches (100 compounds each)
for i in range(0, len(swissadme_input), 100):
    batch = swissadme_input.iloc[i:i+100]
    batch.to_csv(f'swissadme_batch_{i//100:03d}.txt', sep='\t', index=False, header=False)
    print(f"Created batch {i//100}: {len(batch)} compounds")
EOF
```

## Related Resources

- [SuperCYP](../supercyp/) - CYP450 interactions
- [ChEMBL](../../2.2.pharmaceuticals/chembl/) - Bioactivity data
- [DrugBank](../../2.2.pharmaceuticals/drugbank/) - Drug properties
- [ADMETlab](https://admetlab3.scbdd.com/) - Alternative ADME predictor
