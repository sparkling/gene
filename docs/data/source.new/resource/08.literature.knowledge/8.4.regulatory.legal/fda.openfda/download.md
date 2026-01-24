---
id: download-fda-openfda
title: "FDA openFDA Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# FDA openFDA Download Instructions

## Quick Start

```bash
# Query drug adverse events
curl "https://api.fda.gov/drug/event.json?limit=100" -o adverse_events.json
```

## Prerequisites

- **curl** for API access
- **jq** for JSON processing
- API key for higher rate limits (free registration)
- Approximately 1-50GB for bulk downloads

## Registration (Optional, Recommended)

1. Register at https://open.fda.gov/apis/authentication/
2. Receive API key via email
3. Include key in requests for higher limits

## Download Methods

### Method 1: Drug Adverse Event Reports (FAERS)

```bash
# Basic query
curl "https://api.fda.gov/drug/event.json?limit=100" \
  -o adverse_events.json

# Search by drug name
curl "https://api.fda.gov/drug/event.json?search=patient.drug.medicinalproduct:aspirin&limit=100" \
  -o aspirin_events.json

# Search by reaction
curl "https://api.fda.gov/drug/event.json?search=patient.reaction.reactionmeddrapt:nausea&limit=100" \
  -o nausea_events.json

# Date range
curl "https://api.fda.gov/drug/event.json?search=receivedate:[20200101+TO+20201231]&limit=100" \
  -o events_2020.json

# With API key (higher limits)
curl "https://api.fda.gov/drug/event.json?api_key=YOUR_KEY&limit=1000" \
  -o events_batch.json
```

### Method 2: Drug Labeling (SPL)

```bash
# Search drug labels
curl "https://api.fda.gov/drug/label.json?search=openfda.brand_name:advil&limit=10" \
  -o advil_labels.json

# Search by active ingredient
curl "https://api.fda.gov/drug/label.json?search=openfda.substance_name:acetaminophen&limit=100" \
  -o acetaminophen_labels.json

# Get warnings and precautions
curl "https://api.fda.gov/drug/label.json?search=warnings:liver&limit=100" \
  -o liver_warnings.json
```

### Method 3: Drug NDC Directory

```bash
# Query NDC directory
curl "https://api.fda.gov/drug/ndc.json?limit=100" \
  -o ndc_directory.json

# Search by product type
curl "https://api.fda.gov/drug/ndc.json?search=product_type:HUMAN+PRESCRIPTION+DRUG&limit=100" \
  -o prescription_drugs.json

# Search by route
curl "https://api.fda.gov/drug/ndc.json?search=route:ORAL&limit=100" \
  -o oral_drugs.json
```

### Method 4: Device Adverse Events (MAUDE)

```bash
# Device adverse events
curl "https://api.fda.gov/device/event.json?limit=100" \
  -o device_events.json

# Search by device
curl "https://api.fda.gov/device/event.json?search=device.generic_name:pacemaker&limit=100" \
  -o pacemaker_events.json
```

### Method 5: Bulk Downloads

```bash
# Download bulk data files (FAERS)
wget https://download.open.fda.gov/drug/event/drug-event-0001-of-0001.json.zip

# List all available downloads
curl "https://api.fda.gov/download.json" | jq '.results.drug'

# Download specific endpoint bulk file
DOWNLOAD_URL=$(curl -s "https://api.fda.gov/download.json" | \
  jq -r '.results.drug.label.partitions[0].file')
wget "$DOWNLOAD_URL" -O drug_labels_bulk.zip
```

### Method 6: Pagination for Large Queries

```bash
# Paginate through results
python3 << 'EOF'
import requests
import json
import time

base_url = "https://api.fda.gov/drug/event.json"
api_key = "YOUR_API_KEY"  # Optional
limit = 1000
skip = 0
all_results = []

while True:
    url = f"{base_url}?limit={limit}&skip={skip}"
    if api_key:
        url += f"&api_key={api_key}"

    response = requests.get(url)
    data = response.json()

    if 'results' not in data or len(data['results']) == 0:
        break

    all_results.extend(data['results'])
    skip += limit

    if skip >= data.get('meta', {}).get('results', {}).get('total', 0):
        break

    time.sleep(0.5)  # Rate limiting

with open('all_events.json', 'w') as f:
    json.dump(all_results, f)

print(f"Downloaded {len(all_results)} records")
EOF
```

## File Inventory

### Drug Endpoints

| Endpoint | Description | Total Records |
|----------|-------------|---------------|
| drug/event | Adverse events (FAERS) | ~20 million |
| drug/label | Drug labeling (SPL) | ~150,000 |
| drug/ndc | NDC directory | ~200,000 |
| drug/drugsfda | Approved drugs | ~40,000 |
| drug/enforcement | Recalls | ~20,000 |

### Device Endpoints

| Endpoint | Description | Total Records |
|----------|-------------|---------------|
| device/event | MAUDE reports | ~15 million |
| device/classification | Device classes | ~6,000 |
| device/510k | 510(k) clearances | ~200,000 |
| device/pma | PMA approvals | ~50,000 |
| device/recall | Device recalls | ~100,000 |

### Food Endpoints

| Endpoint | Description |
|----------|-------------|
| food/event | Adverse events |
| food/enforcement | Recalls |

## Post-Download Processing

```bash
# Parse FAERS data
python3 << 'EOF'
import json
import pandas as pd

with open('adverse_events.json') as f:
    data = json.load(f)

events = []
for result in data.get('results', []):
    patient = result.get('patient', {})
    for drug in patient.get('drug', []):
        for reaction in patient.get('reaction', []):
            events.append({
                'safetyreportid': result.get('safetyreportid'),
                'drug_name': drug.get('medicinalproduct'),
                'drug_indication': drug.get('drugindication'),
                'reaction': reaction.get('reactionmeddrapt'),
                'outcome': reaction.get('reactionoutcome'),
                'receivedate': result.get('receivedate')
            })

df = pd.DataFrame(events)
df.to_csv('adverse_events_parsed.tsv', sep='\t', index=False)
print(f"Parsed {len(df)} drug-reaction pairs")
EOF

# Extract drug-reaction counts
jq '[.results[].patient.reaction[].reactionmeddrapt] | group_by(.) | map({reaction: .[0], count: length}) | sort_by(-.count)[:20]' \
  adverse_events.json > top_reactions.json

# Link to DrugBank IDs
python3 << 'EOF'
import json

with open('drug_labels.json') as f:
    data = json.load(f)

mappings = []
for result in data.get('results', []):
    openfda = result.get('openfda', {})
    mappings.append({
        'brand_name': openfda.get('brand_name', [''])[0],
        'generic_name': openfda.get('generic_name', [''])[0],
        'rxcui': openfda.get('rxcui', []),
        'unii': openfda.get('unii', []),
        'nui': openfda.get('nui', [])
    })

import pandas as pd
df = pd.DataFrame(mappings)
df.to_csv('drug_mappings.tsv', sep='\t', index=False)
EOF
```

## Verification

```bash
# Check API response
curl -s "https://api.fda.gov/drug/event.json?limit=1" | jq '.meta'

# Count total records for endpoint
curl -s "https://api.fda.gov/drug/event.json?limit=1" | jq '.meta.results.total'

# Verify bulk download
unzip -l drug_labels_bulk.zip
```

## Update Schedule

| Endpoint | Frequency |
|----------|-----------|
| drug/event | Quarterly |
| drug/label | Weekly |
| drug/ndc | Daily |
| drug/enforcement | Weekly |

## Common Issues

- **Rate limits**: 240 requests/minute without key; 120,000/day with key
- **Skip limit**: Maximum skip is 26,000; use date filters for more
- **Missing data**: Not all fields populated; handle nulls
- **Encoding**: Some text fields have HTML entities
- **Bulk file size**: Large JSON files; stream process

## Rate Limits

| With API Key | Without Key |
|--------------|-------------|
| 120,000/day | 1,000/day |
| 240/minute | 40/minute |
| 1,000 per request | 100 per request |

## Query Syntax

| Operator | Example | Description |
|----------|---------|-------------|
| : | field:value | Exact match |
| + | term+term | AND |
| - | -term | NOT |
| "" | "exact phrase" | Phrase search |
| [] | [20200101 TO 20201231] | Range |
| * | aspir* | Wildcard |
| .exact | field.exact | Exact field match |

## Related Resources

- [DrugBank](../../02.compounds.molecules/2.2.pharmaceuticals/drugbank/) - Drug database
- [ClinicalTrials.gov](../clinicaltrials.gov/) - Clinical trials
- [DailyMed](../../02.compounds.molecules/2.2.pharmaceuticals/dailymed/) - Drug labeling
