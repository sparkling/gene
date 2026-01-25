---
id: download-classyfire
title: "ClassyFire Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# ClassyFire Download Instructions

## Quick Start

```bash
# Classify single compound via InChIKey lookup (pre-computed)
curl "http://classyfire.wishartlab.com/entities/BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json" \
  -o aspirin_classification.json

# Classify from SMILES (async, requires polling)
curl -X POST "http://classyfire.wishartlab.com/entities.json" \
  -H "Content-Type: application/json" \
  -d '{"label":"aspirin","smiles":"CC(=O)OC1=CC=CC=C1C(=O)O"}'
```

## Prerequisites

- **curl** or **Python requests** for API calls
- **jq** for JSON parsing (optional)
- InChIKeys for instant lookup (recommended)
- No registration required

## No Registration Required

ClassyFire API is freely accessible. Academic use is free; commercial use requires contact.

## Download Methods

### Method 1: InChIKey Lookup (Recommended)

Pre-computed classifications for 80M+ compounds:

```bash
# Single lookup
curl "http://classyfire.wishartlab.com/entities/BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json" \
  -o classification.json

# Parse with jq
curl -s "http://classyfire.wishartlab.com/entities/BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json" | \
  jq '{kingdom: .kingdom.name, superclass: .superclass.name, class: .class.name}'
```

### Method 2: SMILES Classification (Async)

```bash
# Step 1: Submit compound
RESPONSE=$(curl -s -X POST "http://classyfire.wishartlab.com/entities.json" \
  -H "Content-Type: application/json" \
  -d '{"label":"my_compound","smiles":"CC(=O)OC1=CC=CC=C1C(=O)O"}')

QUERY_ID=$(echo $RESPONSE | jq -r '.id')
echo "Query ID: $QUERY_ID"

# Step 2: Poll for completion (wait for "Done" status)
while true; do
  STATUS=$(curl -s "http://classyfire.wishartlab.com/entities/$QUERY_ID/status.json" | jq -r '.status')
  echo "Status: $STATUS"
  if [ "$STATUS" = "Done" ]; then break; fi
  sleep 5
done

# Step 3: Get results
curl -s "http://classyfire.wishartlab.com/entities/$QUERY_ID.json" -o result.json
```

### Method 3: Batch Classification

```bash
# Submit batch (newline-separated SMILES)
curl -X POST "http://classyfire.wishartlab.com/queries.json" \
  -H "Content-Type: application/json" \
  -d '{
    "label": "my_batch",
    "query_input": "CC(=O)OC1=CC=CC=C1C(=O)O\nCC(=O)NC1=CC=C(C=C1)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C"
  }'

# Check batch status
curl "http://classyfire.wishartlab.com/queries/{QUERY_ID}/status.json"

# Get all entities
curl "http://classyfire.wishartlab.com/queries/{QUERY_ID}/entities.json" -o batch_results.json
```

### Method 4: Python Batch Processor

```python
#!/usr/bin/env python3
"""
ClassyFire batch classification with retry logic
"""
import requests
import json
import time
import pandas as pd

API_BASE = "http://classyfire.wishartlab.com"

def classify_inchikey(inchikey):
    """Instant lookup for pre-computed classification"""
    url = f"{API_BASE}/entities/{inchikey}.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return None

def classify_smiles_batch(smiles_list, label="batch"):
    """Submit batch of SMILES for classification"""
    url = f"{API_BASE}/queries.json"

    # Join SMILES with newlines
    query_input = "\n".join(smiles_list)

    payload = {
        "label": label,
        "query_input": query_input
    }

    response = requests.post(url, json=payload)
    if response.status_code in [200, 201]:
        return response.json()['id']
    raise Exception(f"Submission failed: {response.status_code}")

def wait_for_completion(query_id, max_wait=600, poll_interval=10):
    """Poll until classification complete"""
    url = f"{API_BASE}/queries/{query_id}/status.json"
    elapsed = 0

    while elapsed < max_wait:
        response = requests.get(url)
        data = response.json()

        if data['status'] == 'Done':
            return True
        elif data['status'] == 'Error':
            raise Exception("Classification failed")

        time.sleep(poll_interval)
        elapsed += poll_interval
        print(f"  Waiting... {elapsed}s")

    raise TimeoutError("Classification timed out")

def get_batch_results(query_id):
    """Retrieve batch classification results"""
    url = f"{API_BASE}/queries/{query_id}/entities.json"
    response = requests.get(url)
    return response.json()

def extract_classification(entity):
    """Extract key classification fields"""
    return {
        'smiles': entity.get('smiles'),
        'inchikey': entity.get('inchikey'),
        'kingdom': entity.get('kingdom', {}).get('name'),
        'superclass': entity.get('superclass', {}).get('name'),
        'class': entity.get('class', {}).get('name'),
        'subclass': entity.get('subclass', {}).get('name') if entity.get('subclass') else None,
        'direct_parent': entity.get('direct_parent', {}).get('name'),
        'kingdom_id': entity.get('kingdom', {}).get('chemont_id'),
        'superclass_id': entity.get('superclass', {}).get('chemont_id'),
        'class_id': entity.get('class', {}).get('chemont_id'),
    }

# Example usage
if __name__ == "__main__":
    # Test compounds
    smiles_list = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CC(=O)NC1=CC=C(C=C1)O",      # Paracetamol
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    ]

    print("Submitting batch...")
    query_id = classify_smiles_batch(smiles_list, "test_batch")
    print(f"Query ID: {query_id}")

    print("Waiting for completion...")
    wait_for_completion(query_id)

    print("Fetching results...")
    results = get_batch_results(query_id)

    # Convert to DataFrame
    classifications = [extract_classification(e) for e in results]
    df = pd.DataFrame(classifications)
    df.to_csv('classyfire_results.csv', index=False)
    print(f"Saved {len(df)} classifications")
```

### Method 5: Large-Scale InChIKey Lookup

```python
#!/usr/bin/env python3
"""
Efficient bulk InChIKey lookup for ClassyFire
"""
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from time import sleep

API_BASE = "http://classyfire.wishartlab.com"

def lookup_inchikey(inchikey, retries=3):
    """Lookup with retry logic"""
    for attempt in range(retries):
        try:
            url = f"{API_BASE}/entities/{inchikey}.json"
            response = requests.get(url, timeout=30)

            if response.status_code == 200:
                data = response.json()
                return {
                    'inchikey': inchikey,
                    'kingdom': data.get('kingdom', {}).get('name'),
                    'superclass': data.get('superclass', {}).get('name'),
                    'class': data.get('class', {}).get('name'),
                    'subclass': data.get('subclass', {}).get('name') if data.get('subclass') else None,
                    'direct_parent': data.get('direct_parent', {}).get('name'),
                    'status': 'found'
                }
            elif response.status_code == 404:
                return {'inchikey': inchikey, 'status': 'not_found'}

        except Exception as e:
            if attempt == retries - 1:
                return {'inchikey': inchikey, 'status': 'error', 'error': str(e)}
            sleep(1)

    return {'inchikey': inchikey, 'status': 'failed'}

def bulk_lookup(inchikeys, workers=5, delay=0.2):
    """
    Parallel InChIKey lookup

    Args:
        inchikeys: List of InChIKeys
        workers: Parallel workers (keep low to respect rate limits)
        delay: Delay between requests
    """
    results = []

    with ThreadPoolExecutor(max_workers=workers) as executor:
        for i, result in enumerate(executor.map(lookup_inchikey, inchikeys)):
            results.append(result)
            if (i + 1) % 100 == 0:
                print(f"Processed {i + 1}/{len(inchikeys)}")
            sleep(delay)

    return results

# Example
if __name__ == "__main__":
    # Load InChIKeys from file
    df = pd.read_csv('compounds.tsv', sep='\t')
    inchikeys = df['inchikey'].dropna().tolist()[:1000]  # Limit for testing

    print(f"Looking up {len(inchikeys)} InChIKeys...")
    results = bulk_lookup(inchikeys)

    results_df = pd.DataFrame(results)
    results_df.to_csv('classyfire_lookups.csv', index=False)

    # Statistics
    found = results_df[results_df['status'] == 'found']
    print(f"Found: {len(found)}/{len(results)} ({100*len(found)/len(results):.1f}%)")
```

### Method 6: Taxonomy Browser Download

```bash
# Download taxonomy structure
python3 << 'EOF'
import requests
import json
from collections import deque

API_BASE = "http://classyfire.wishartlab.com"

def get_taxonomy_node(chemont_id):
    """Get single taxonomy node"""
    url = f"{API_BASE}/tax_nodes/{chemont_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return None

def get_children(chemont_id):
    """Get child nodes"""
    url = f"{API_BASE}/tax_nodes/{chemont_id}/children.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return []

def download_taxonomy(root_id="C0000000", max_depth=3):
    """Download taxonomy tree starting from root"""
    taxonomy = {}
    queue = deque([(root_id, 0)])

    while queue:
        node_id, depth = queue.popleft()
        if depth > max_depth:
            continue

        node = get_taxonomy_node(node_id)
        if node:
            taxonomy[node_id] = {
                'name': node.get('name'),
                'description': node.get('description'),
                'chemont_id': node.get('chemont_id'),
                'depth': depth
            }

            if depth < max_depth:
                children = get_children(node_id)
                for child in children:
                    child_id = child.get('chemont_id', '').replace('CHEMONTID:', 'C')
                    if child_id:
                        queue.append((child_id, depth + 1))

        print(f"Downloaded: {node_id} (depth {depth})")

    return taxonomy

# Download top 3 levels
taxonomy = download_taxonomy(max_depth=2)

with open('chemont_taxonomy.json', 'w') as f:
    json.dump(taxonomy, f, indent=2)

print(f"Downloaded {len(taxonomy)} taxonomy nodes")
EOF
```

## File Inventory

ClassyFire is primarily API-based; no bulk database downloads available.

| Resource | Type | Description |
|----------|------|-------------|
| REST API | Service | Classification endpoint |
| InChIKey Index | Pre-computed | 80M+ compounds |
| ChemOnt Ontology | Taxonomy | 4,800+ categories |

## Post-Download Processing

### Parse Classification Results

```python
import pandas as pd
import json

# Load batch results
with open('classyfire_results.json') as f:
    data = json.load(f)

# Extract hierarchy
rows = []
for entity in data:
    rows.append({
        'smiles': entity.get('smiles'),
        'inchikey': entity.get('inchikey'),
        'kingdom': entity.get('kingdom', {}).get('name'),
        'superclass': entity.get('superclass', {}).get('name'),
        'class': entity.get('class', {}).get('name'),
        'subclass': entity.get('subclass', {}).get('name') if entity.get('subclass') else None,
        'direct_parent': entity.get('direct_parent', {}).get('name'),
        'substituents': '|'.join(entity.get('substituents', [])),
        'molecular_framework': entity.get('molecular_framework')
    })

df = pd.DataFrame(rows)
df.to_csv('classifications_flat.csv', index=False)

# Summary by superclass
print(df['superclass'].value_counts())
```

### Build Classification Graph

```python
import networkx as nx
import json

with open('classyfire_results.json') as f:
    data = json.load(f)

G = nx.DiGraph()

for entity in data:
    # Add nodes for each level
    if entity.get('kingdom'):
        G.add_node(entity['kingdom']['chemont_id'],
                   name=entity['kingdom']['name'],
                   level='kingdom')

    if entity.get('superclass'):
        G.add_node(entity['superclass']['chemont_id'],
                   name=entity['superclass']['name'],
                   level='superclass')
        G.add_edge(entity['kingdom']['chemont_id'],
                   entity['superclass']['chemont_id'])

    # Add more levels...

# Export
nx.write_graphml(G, 'classyfire_hierarchy.graphml')
```

## Verification

```bash
# Check API availability
curl -I "http://classyfire.wishartlab.com/entities/BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json"

# Validate JSON response
curl -s "http://classyfire.wishartlab.com/entities/BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json" | jq .

# Check classification has expected fields
curl -s "http://classyfire.wishartlab.com/entities/BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json" | \
  jq 'keys'
```

```python
# Verify classification accuracy
import requests

# Known classifications
test_cases = [
    ("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "Aspirin", "Phenyl acetates"),  # Aspirin
    ("QGZKDVFQNNGYKY-UHFFFAOYSA-N", "Glutamine", "Alpha-amino acids"),
    ("RYYVLZVUVIJVGH-UHFFFAOYSA-N", "Caffeine", "Xanthines"),
]

for inchikey, name, expected_parent in test_cases:
    url = f"http://classyfire.wishartlab.com/entities/{inchikey}.json"
    response = requests.get(url)
    data = response.json()
    actual = data.get('direct_parent', {}).get('name', 'Unknown')
    status = "OK" if expected_parent in actual else "MISMATCH"
    print(f"{name}: Expected '{expected_parent}', Got '{actual}' [{status}]")
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| ChemOnt 2.1 | 2016-11 | ~10 MB | Current |
| Pre-computed Index | Continuous | 80M+ compounds | Current |

### Version Notes

ClassyFire provides automated chemical classification:
- ChemOnt taxonomy with 4,825+ chemical classes
- 11 hierarchical levels (Kingdom to direct parent)
- 80M+ pre-computed compound classifications via InChIKey
- Ruby API available for programmatic access

## API Access

| Property | Value |
|----------|-------|
| Base URL | `http://classyfire.wishartlab.com` |
| Rate Limit | ~1000 InChIKey lookups/hour |
| Auth Required | No |
| Documentation | http://classyfire.wishartlab.com |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Pre-computed Index | Monthly |
| ChemOnt Updates | As needed |
| API | Continuous |

## Rate Limits

| Limit Type | Value |
|------------|-------|
| InChIKey lookups | ~1000/hour |
| Batch submissions | ~100 compounds/batch |
| Parallel requests | Keep <5 concurrent |

## Common Issues

- **429 Rate Limited**: Reduce request frequency
- **Classification timeout**: Use InChIKey lookup instead
- **Missing subclass**: Not all compounds have all levels
- **Alternative parents**: Check alternative_parents for complete picture
- **Special characters**: URL-encode SMILES in GET requests

## Integration Examples

```python
# Map ChEMBL compounds to ClassyFire
import pandas as pd

# Load ChEMBL
chembl = pd.read_csv('chembl_compounds.tsv', sep='\t')

# Classify via InChIKey
from classyfire_utils import bulk_lookup

inchikeys = chembl['standard_inchi_key'].dropna().tolist()
classifications = bulk_lookup(inchikeys[:1000])

# Merge
chembl_classified = chembl.merge(
    pd.DataFrame(classifications),
    left_on='standard_inchi_key',
    right_on='inchikey'
)

chembl_classified.to_csv('chembl_classyfire.csv', index=False)
```

## Related Resources

- [ChEBI](../chebi/) - Chemical ontology
- [NPClassifier](../npclassifier/) - Natural product classification
- [PubChem](../pubchem/) - Chemical repository
