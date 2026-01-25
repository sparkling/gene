---
id: download-kampodb
title: "KampoDB Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# KampoDB Download Instructions

## Quick Start

```bash
# REST API for programmatic access
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/"

# List all formulas
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/"

# Get compounds in Kakkonto
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/KT/compound"
```

## Prerequisites

- Python 3.8+ with requests library
- cURL or similar HTTP client
- JSON processing tools

## Download Methods

### Primary: REST API (Recommended)

KampoDB provides comprehensive REST API with JSON responses:

**Base URL**: https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/

### API Endpoint Categories

| Category | Endpoint | Description |
|----------|----------|-------------|
| List All | `/api/{entity}/` | All IDs and names |
| Entity Info | `/api/{entity}/{id}/info` | Detailed information |
| Relations | `/api/{entity}/{id}/{target}` | Hierarchy navigation |
| Enrichment | `/api/{entity}/{id}/pathway` | KEGG annotations |
| Docking | `/api/docking/compound/{id}/protein/{id}` | Binding data |

## File Inventory

| Entity | Count | Description |
|--------|-------|-------------|
| Kampo Formulas | 298 | Traditional prescriptions |
| Crude Drugs | 180 | Medicinal materials |
| Compounds | 3,002 | PubChem CIDs |
| Proteins | 62,906 | NCBI Gene IDs |
| Docking Simulations | 3,063,505 | Binding predictions |

## API Examples

```bash
# List all formulas
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/"

# Get formula info
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/KT/info"

# Get compounds in formula
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/KT/compound"

# Get docking score
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/docking/compound/969516/protein/2475"

# Get protein pathways
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/protein/2475/pathway"
```

## Python Download Script

```python
import requests
import json
import time

BASE_URL = "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api"

def get_all_formulas():
    """Get list of all Kampo formulas."""
    response = requests.get(f"{BASE_URL}/formula/")
    return response.json()

def get_formula_compounds(formula_id):
    """Get all compounds in a formula."""
    response = requests.get(f"{BASE_URL}/formula/{formula_id}/compound")
    return response.json()

def get_compound_targets(compound_id):
    """Get target proteins for a compound."""
    response = requests.get(f"{BASE_URL}/compound/{compound_id}/protein")
    return response.json()

def get_docking_score(compound_id, protein_id):
    """Get molecular docking results."""
    response = requests.get(
        f"{BASE_URL}/docking/compound/{compound_id}/protein/{protein_id}"
    )
    return response.json()

# Download all formulas
formulas = get_all_formulas()
print(f"Total formulas: {len(formulas)}")

# Download compounds for each formula
all_compounds = []
for formula in formulas[:5]:  # Sample first 5
    compounds = get_formula_compounds(formula['id'])
    all_compounds.extend(compounds)
    time.sleep(1)  # Be respectful to server

print(f"Sample compounds: {len(all_compounds)}")
```

## Bulk Download Strategy

Since no bulk download endpoint exists:

1. Query `/api/formula/` for all formula IDs
2. Iterate through formulas to get compounds
3. Query compounds for targets
4. Build local database from API responses
5. Cache responses to minimize requests

```python
import json

# Save API responses locally
def cache_response(endpoint, data):
    filename = endpoint.replace('/', '_') + '.json'
    with open(f'kampodb_cache/{filename}', 'w') as f:
        json.dump(data, f)

# Build comprehensive dataset
formulas = get_all_formulas()
cache_response('formulas', formulas)

for formula in formulas:
    compounds = get_formula_compounds(formula['id'])
    cache_response(f"formula_{formula['id']}_compounds", compounds)
    time.sleep(0.5)  # Rate limiting
```

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Formulas | 298 |
| Crude Drugs | 180 |
| Compounds | 3,002 |
| Proteins | 62,906 |
| Docking Simulations | 3,063,505 |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Initial Release | 2018 |
| Status | Active |
| License | CC BY-SA 4.0 |

---

## Dataset Versions

### Current Release: KampoDB

| Property | Value |
|----------|-------|
| Version | 1.0 |
| Release Date | 2018-01-01 |
| Total Size | ~500 MB (via API) |
| Docking Simulations | 3,063,505 |

### Version Contents

| Component | Records | Description |
|-----------|---------|-------------|
| Kampo Formulas | 298 | Traditional prescriptions |
| Crude Drugs | 180 | Medicinal materials |
| Compounds | 3,002 | PubChem CIDs |
| Proteins | 62,906 | NCBI Gene IDs |
| Docking Simulations | 3,063,505 | Binding predictions |

### Data Characteristics

| Property | Value |
|----------|-------|
| Compound IDs | PubChem CIDs |
| Protein IDs | NCBI Gene IDs |
| Hierarchy | 4-layer (Formula > Drug > Compound > Protein) |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| Base URL | `https://wakanmoview.inm.u-toyama.ac.jp/kampo/api` |
| Authentication | None required |
| Rate Limit | Use reasonable delays |
| Response Format | JSON |

### API Endpoints

| Operation | Endpoint | Example |
|-----------|----------|---------|
| List Formulas | `/formula/` | List all formula IDs |
| Formula Info | `/formula/{id}/info` | `/formula/KT/info` |
| Formula Compounds | `/formula/{id}/compound` | `/formula/KT/compound` |
| Compound Targets | `/compound/{id}/protein` | Get protein targets |
| Protein Pathways | `/protein/{id}/pathway` | KEGG annotations |
| Docking | `/docking/compound/{cid}/protein/{pid}` | Binding data |

---

## Notes

- No bulk download endpoint (query individual entities)
- Rate limiting not documented; use reasonable delays
- Compound IDs are PubChem CIDs
- Gene IDs require manual UniProt mapping
- 4-layer hierarchy fully traversable via API
- Extensive docking data (3M+ simulations)
