---
id: download-npclassifier
title: "NPClassifier Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# NPClassifier Download Instructions

## Quick Start

```bash
# Classify single compound via API
curl "https://npclassifier.ucsd.edu/classify?smiles=CC1%3DCC2%3DC(C%3DC1)NC%3DC2CCN"

# Clone and run locally
git clone https://github.com/mwang87/NP-Classifier.git
cd NP-Classifier
pip install -r requirements.txt
```

## Prerequisites

- **curl** or **Python requests** for API access
- **Python 3.7+** for local installation
- **PyTorch** for local model inference
- **RDKit** for structure handling
- No registration required

## No Registration Required

NPClassifier API is freely accessible under MIT license for all users.

## Download Methods

### Method 1: REST API (Single Compound)

```bash
# Simple GET request
curl "https://npclassifier.ucsd.edu/classify?smiles=CC1%3DCC2%3DC(C%3DC1)NC%3DC2CCN"

# With proper URL encoding (Python)
python3 -c "
import urllib.parse
smiles = 'CC1=CC2=C(C=C1)NC=C2CCN'
encoded = urllib.parse.quote(smiles, safe='')
print(f'https://npclassifier.ucsd.edu/classify?smiles={encoded}')
"

# Save output
curl -s "https://npclassifier.ucsd.edu/classify?smiles=CC1%3DCC2%3DC(C%3DC1)NC%3DC2CCN" \
  -o tryptamine_class.json
```

### Method 2: Python API Client

```python
#!/usr/bin/env python3
"""
NPClassifier API client with batch processing
"""
import requests
import urllib.parse
import pandas as pd
from time import sleep

API_URL = "https://npclassifier.ucsd.edu/classify"

def classify_smiles(smiles, retries=3):
    """Classify single SMILES via API"""
    encoded = urllib.parse.quote(smiles, safe='')
    url = f"{API_URL}?smiles={encoded}"

    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 400:
                return {'error': 'Invalid SMILES', 'smiles': smiles}
        except Exception as e:
            if attempt == retries - 1:
                return {'error': str(e), 'smiles': smiles}
            sleep(1)

    return {'error': 'Failed after retries', 'smiles': smiles}

def extract_top_predictions(result):
    """Extract top predictions from result"""
    if 'error' in result:
        return {
            'smiles': result.get('smiles'),
            'error': result.get('error')
        }

    # Get top pathway
    pathways = result.get('pathway_results', [])
    top_pathway = pathways[0] if pathways else {}

    # Get top superclass
    superclasses = result.get('superclass_results', [])
    top_superclass = superclasses[0] if superclasses else {}

    # Get top class
    classes = result.get('class_results', [])
    top_class = classes[0] if classes else {}

    return {
        'smiles': result.get('smiles'),
        'pathway': top_pathway.get('pathway'),
        'pathway_score': top_pathway.get('pathway_score'),
        'superclass': top_superclass.get('superclass'),
        'superclass_score': top_superclass.get('superclass_score'),
        'class': top_class.get('class'),
        'class_score': top_class.get('class_score'),
        'isglycoside': result.get('isglycoside')
    }

def batch_classify(smiles_list, delay=0.5):
    """Classify list of SMILES with rate limiting"""
    results = []
    for i, smiles in enumerate(smiles_list):
        result = classify_smiles(smiles)
        results.append(extract_top_predictions(result))

        if (i + 1) % 10 == 0:
            print(f"Processed {i + 1}/{len(smiles_list)}")

        sleep(delay)

    return results

# Example usage
if __name__ == "__main__":
    test_compounds = [
        "CC1=CC2=C(C=C1)NC=C2CCN",  # Tryptamine
        "COC1=CC(=CC(=C1O)OC)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC",  # Curcumin
        "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",  # Cholesterol
    ]

    results = batch_classify(test_compounds)
    df = pd.DataFrame(results)
    df.to_csv('npclassifier_results.csv', index=False)
    print(df)
```

### Method 3: Local Installation

```bash
# Clone repository
git clone https://github.com/mwang87/NP-Classifier.git
cd NP-Classifier

# Create virtual environment
python3 -m venv npclass_env
source npclass_env/bin/activate

# Install dependencies
pip install -r requirements.txt

# Download model weights (if not included)
# Models are typically in the models/ directory

# Test installation
python classify.py --smiles "CC1=CC2=C(C=C1)NC=C2CCN"
```

### Method 4: Local Batch Processing

```python
#!/usr/bin/env python3
"""
Local NPClassifier batch processing (requires local installation)
"""
import sys
sys.path.append('/path/to/NP-Classifier')

from NPClassifier import NPClassifier
import pandas as pd

# Initialize classifier
classifier = NPClassifier()

def classify_file(input_file, output_file, smiles_col='SMILES'):
    """Process file of compounds"""
    df = pd.read_csv(input_file, sep='\t')

    results = []
    for idx, row in df.iterrows():
        smiles = row[smiles_col]
        try:
            result = classifier.classify(smiles)
            results.append({
                'input_smiles': smiles,
                'pathway': result['pathway_results'][0]['pathway'] if result['pathway_results'] else None,
                'pathway_score': result['pathway_results'][0]['pathway_score'] if result['pathway_results'] else None,
                'superclass': result['superclass_results'][0]['superclass'] if result['superclass_results'] else None,
                'class': result['class_results'][0]['class'] if result['class_results'] else None,
                'isglycoside': result.get('isglycoside')
            })
        except Exception as e:
            results.append({
                'input_smiles': smiles,
                'error': str(e)
            })

        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1}/{len(df)}")

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    return results_df

# Example
# classify_file('natural_products.tsv', 'np_classified.csv')
```

### Method 5: GNPS Integration

For mass spectrometry data:

```bash
# NPClassifier is integrated into GNPS workflows
# Access via: https://gnps.ucsd.edu

# After running GNPS molecular networking:
# 1. Go to "View All Spectra with Annotation"
# 2. NPClassifier predictions are included for annotated compounds

# Direct GNPS library spectrum classification
curl "https://gnps.ucsd.edu/ProteoSAFe/gnpslibrary.jsp?library=GNPS-LIBRARY&compound=CCMSLIB00000001"
```

### Method 6: Parallel API Processing

```python
#!/usr/bin/env python3
"""
Parallel NPClassifier API calls with async
"""
import asyncio
import aiohttp
import urllib.parse
import pandas as pd

API_URL = "https://npclassifier.ucsd.edu/classify"

async def classify_async(session, smiles, semaphore):
    """Async classification with semaphore for rate limiting"""
    async with semaphore:
        encoded = urllib.parse.quote(smiles, safe='')
        url = f"{API_URL}?smiles={encoded}"

        try:
            async with session.get(url, timeout=30) as response:
                if response.status == 200:
                    data = await response.json()
                    return {'smiles': smiles, 'result': data}
                else:
                    return {'smiles': smiles, 'error': f'HTTP {response.status}'}
        except Exception as e:
            return {'smiles': smiles, 'error': str(e)}

async def batch_classify_async(smiles_list, max_concurrent=5):
    """Process batch with limited concurrency"""
    semaphore = asyncio.Semaphore(max_concurrent)

    async with aiohttp.ClientSession() as session:
        tasks = [
            classify_async(session, smiles, semaphore)
            for smiles in smiles_list
        ]
        results = await asyncio.gather(*tasks)

    return results

def process_results(results):
    """Convert async results to DataFrame"""
    rows = []
    for r in results:
        if 'error' in r:
            rows.append({'smiles': r['smiles'], 'error': r['error']})
        else:
            data = r['result']
            pathways = data.get('pathway_results', [])
            superclasses = data.get('superclass_results', [])
            classes = data.get('class_results', [])

            rows.append({
                'smiles': r['smiles'],
                'pathway': pathways[0]['pathway'] if pathways else None,
                'pathway_score': pathways[0]['pathway_score'] if pathways else None,
                'superclass': superclasses[0]['superclass'] if superclasses else None,
                'class': classes[0]['class'] if classes else None,
                'isglycoside': data.get('isglycoside')
            })

    return pd.DataFrame(rows)

# Example usage
if __name__ == "__main__":
    # Load SMILES
    smiles_list = [
        "CC1=CC2=C(C=C1)NC=C2CCN",
        "COC1=C(C=CC(=C1)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC)O",
        "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
    ]

    # Run async batch
    results = asyncio.run(batch_classify_async(smiles_list))
    df = process_results(results)
    print(df)
```

## File Inventory

| Resource | Type | Description |
|----------|------|-------------|
| REST API | Service | Classification endpoint |
| GitHub Repo | Source | Model code and weights |
| GNPS | Integration | Mass spec workflows |

### GitHub Repository Contents

| File/Directory | Description |
|----------------|-------------|
| models/ | Trained model weights |
| NPClassifier.py | Main classifier class |
| classify.py | CLI interface |
| requirements.txt | Python dependencies |
| notebooks/ | Example Jupyter notebooks |

## Post-Download Processing

### Aggregate Classification Statistics

```python
import pandas as pd

# Load classifications
df = pd.read_csv('npclassifier_results.csv')

# Pathway distribution
print("Pathway Distribution:")
print(df['pathway'].value_counts())

# Superclass distribution
print("\nSuperclass Distribution:")
print(df['superclass'].value_counts().head(20))

# Glycoside percentage
glycoside_pct = df['isglycoside'].mean() * 100
print(f"\nGlycoside percentage: {glycoside_pct:.1f}%")

# High confidence predictions
high_conf = df[df['pathway_score'] > 0.9]
print(f"\nHigh confidence (>0.9): {len(high_conf)}/{len(df)} ({100*len(high_conf)/len(df):.1f}%)")
```

### Combine with Other Databases

```python
import pandas as pd

# Load NPClassifier results
npc = pd.read_csv('npclassifier_results.csv')

# Load COCONUT natural products
coconut = pd.read_csv('coconut_compounds.tsv', sep='\t')

# Merge on SMILES (after canonicalization)
from rdkit import Chem

def canonicalize(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

npc['canonical_smiles'] = npc['smiles'].apply(canonicalize)
coconut['canonical_smiles'] = coconut['SMILES'].apply(canonicalize)

merged = coconut.merge(npc, on='canonical_smiles', how='left')
print(f"Merged: {len(merged)} records")

# Analyze pathway by source organism
if 'source_organism' in merged.columns:
    pathway_by_org = merged.groupby(['source_organism', 'pathway']).size().unstack(fill_value=0)
    print(pathway_by_org.head(10))
```

### Build Pathway Network

```python
import networkx as nx
import pandas as pd

df = pd.read_csv('npclassifier_results.csv')

# Create hierarchy network
G = nx.DiGraph()

for _, row in df.iterrows():
    if pd.notna(row['pathway']):
        G.add_node(row['pathway'], level='pathway')
    if pd.notna(row['superclass']):
        G.add_node(row['superclass'], level='superclass')
        if pd.notna(row['pathway']):
            G.add_edge(row['pathway'], row['superclass'])
    if pd.notna(row['class']):
        G.add_node(row['class'], level='class')
        if pd.notna(row['superclass']):
            G.add_edge(row['superclass'], row['class'])

# Export
nx.write_graphml(G, 'npclassifier_hierarchy.graphml')
print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
```

## Verification

```bash
# Test API availability
curl -I "https://npclassifier.ucsd.edu/classify?smiles=CCO"

# Validate response structure
curl -s "https://npclassifier.ucsd.edu/classify?smiles=CC1%3DCC2%3DC(C%3DC1)NC%3DC2CCN" | \
  python3 -c "import sys, json; d=json.load(sys.stdin); print('pathway_results:', len(d.get('pathway_results',[]))); print('superclass_results:', len(d.get('superclass_results',[]))); print('class_results:', len(d.get('class_results',[])))"
```

```python
# Verify known classifications
import requests
import urllib.parse

test_cases = [
    ("CC1=CC2=C(C=C1)NC=C2CCN", "Alkaloids", "Tryptamine"),
    ("CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C", "Terpenoids", "Cholesterol"),
    ("OC1=CC=C(C=C1)C2=CC(=O)C3=C(O)C=C(O)C=C3O2", "Shikimates and Phenylpropanoids", "Apigenin"),
]

for smiles, expected_pathway, name in test_cases:
    encoded = urllib.parse.quote(smiles, safe='')
    url = f"https://npclassifier.ucsd.edu/classify?smiles={encoded}"
    response = requests.get(url)
    data = response.json()

    actual = data['pathway_results'][0]['pathway'] if data.get('pathway_results') else 'None'
    status = "OK" if expected_pathway == actual else "MISMATCH"
    print(f"{name}: Expected '{expected_pathway}', Got '{actual}' [{status}]")
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| NPClassifier v1.0+ | 2021 | ~200 MB model | Current |
| Updates | Ongoing | N/A | Active |

### Version Notes

NPClassifier is a deep neural network classifier for natural products:
- 7 Pathways: Fatty acids, Polyketides, Shikimates, Terpenoids, Alkaloids, Amino acids, Carbohydrates
- 70 Superclasses
- 672 Classes
- Trained on 73,607 NPs from PubChem, ChEBI, UNPD, etc.
- Glycoside detection included

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://npclassifier.ucsd.edu/classify` |
| Rate Limit | 1-2 req/sec recommended |
| Auth Required | No |
| Documentation | https://github.com/mwang87/NP-Classifier |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Model Updates | As published |
| API | Continuous |
| Training Data | Annual (COCONUT/LOTUS) |

## Rate Limits

| Limit Type | Recommendation |
|------------|----------------|
| Sequential requests | 1-2 per second |
| Parallel requests | Max 5 concurrent |
| Batch size | 100-1000 per session |

## Common Issues

- **Invalid SMILES**: Validate with RDKit before submission
- **URL encoding**: Use proper URL encoding for special characters
- **Mixed pathways**: Compounds can have multiple pathway origins
- **Low scores**: May indicate non-natural product or novel scaffold
- **Timeout**: Use local installation for large batches

## Integration with GNPS

```python
# After GNPS molecular networking job completes
# Download feature quantification table with annotations

import requests
import pandas as pd

# Get GNPS job results (replace with actual task ID)
GNPS_TASK = "your_gnps_task_id"

# Download annotation results
url = f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={GNPS_TASK}&view=view_all_annotations_DB"
annotations = requests.get(url).json()

# NPClassifier predictions are included for library matches
for ann in annotations['blockData']:
    if 'NPClassifier_pathway' in ann:
        print(f"Compound: {ann.get('Compound_Name')}")
        print(f"  Pathway: {ann.get('NPClassifier_pathway')}")
        print(f"  Superclass: {ann.get('NPClassifier_superclass')}")
```

## Related Resources

- [ClassyFire](../classyfire/) - General chemical taxonomy
- [COCONUT](../../2.1.natural.products/coconut/) - Natural products database
- [GNPS](https://gnps.ucsd.edu) - Mass spectrometry platform
- [ChEBI](../chebi/) - Chemical ontology
