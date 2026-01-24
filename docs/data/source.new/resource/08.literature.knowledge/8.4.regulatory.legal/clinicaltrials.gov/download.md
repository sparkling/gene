---
id: download-clinicaltrials
title: "ClinicalTrials.gov Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# ClinicalTrials.gov Download Instructions

## Quick Start

```bash
# Download all studies via API
curl "https://clinicaltrials.gov/api/v2/studies?format=json&pageSize=1000" -o studies_page1.json
```

## Prerequisites

- **curl** or **wget** for downloads
- **jq** for JSON processing
- Approximately 5-20GB disk space for full database

## No Registration Required

ClinicalTrials.gov data is public domain.

## Download Methods

### Method 1: REST API (v2, Recommended)

```bash
# Get studies in JSON format
curl "https://clinicaltrials.gov/api/v2/studies?format=json&pageSize=100" \
  -o studies.json

# Search by condition
curl "https://clinicaltrials.gov/api/v2/studies?query.cond=breast+cancer&format=json&pageSize=100" \
  -o breast_cancer_trials.json

# Search by intervention
curl "https://clinicaltrials.gov/api/v2/studies?query.intr=pembrolizumab&format=json&pageSize=100" \
  -o pembrolizumab_trials.json

# Filter by status
curl "https://clinicaltrials.gov/api/v2/studies?filter.overallStatus=RECRUITING&format=json&pageSize=100" \
  -o recruiting_trials.json

# Get specific trial
curl "https://clinicaltrials.gov/api/v2/studies/NCT04368728" \
  -o NCT04368728.json
```

### Method 2: Bulk Download (All Trials)

```bash
# Paginate through all studies
PAGE_TOKEN=""
PAGE=1

while true; do
  if [ -z "$PAGE_TOKEN" ]; then
    RESPONSE=$(curl -s "https://clinicaltrials.gov/api/v2/studies?format=json&pageSize=1000")
  else
    RESPONSE=$(curl -s "https://clinicaltrials.gov/api/v2/studies?format=json&pageSize=1000&pageToken=${PAGE_TOKEN}")
  fi

  echo "$RESPONSE" > "studies_page${PAGE}.json"

  PAGE_TOKEN=$(echo "$RESPONSE" | jq -r '.nextPageToken // empty')
  if [ -z "$PAGE_TOKEN" ]; then
    break
  fi

  PAGE=$((PAGE + 1))
  sleep 1
done
```

### Method 3: Download by Date Range

```bash
# Studies updated in date range
curl "https://clinicaltrials.gov/api/v2/studies?filter.lastUpdatePostDate=2024-01-01,2024-12-31&format=json&pageSize=1000" \
  -o studies_2024.json

# Studies started in date range
curl "https://clinicaltrials.gov/api/v2/studies?filter.startDate=2020-01-01,2024-12-31&format=json&pageSize=1000" \
  -o studies_started_2020_2024.json
```

### Method 4: XML Download (Legacy)

```bash
# Download studies in XML format
curl "https://clinicaltrials.gov/ct2/results/download_studies?down_count=10000&down_flds=all&down_fmt=xml" \
  -o studies.zip

unzip studies.zip -d studies_xml/
```

### Method 5: CSV Export (via Web Interface)

1. Navigate to https://clinicaltrials.gov/search
2. Apply desired filters
3. Click "Download" button
4. Select CSV format

### Method 6: Study Records Data API

```bash
# Get study protocol
curl "https://clinicaltrials.gov/api/v2/studies/NCT04368728/protocol" \
  -o protocol.json

# Get study results
curl "https://clinicaltrials.gov/api/v2/studies/NCT04368728/results" \
  -o results.json
```

## File Inventory

### API Downloads

| Query Type | Estimated Records |
|------------|-------------------|
| All studies | ~500,000 |
| Recruiting | ~50,000 |
| Completed | ~150,000 |
| With results | ~60,000 |

### Data Fields

| Field Category | Description |
|----------------|-------------|
| Protocol | Study design, arms, outcomes |
| Status | Current status, dates |
| Sponsors | Sponsors, collaborators |
| Eligibility | Criteria, demographics |
| Results | Outcome measures, adverse events |

## Post-Download Processing

```bash
# Parse JSON responses
python3 << 'EOF'
import json
import pandas as pd
import glob

all_studies = []

for json_file in glob.glob('studies_page*.json'):
    with open(json_file) as f:
        data = json.load(f)
        for study in data.get('studies', []):
            protocol = study.get('protocolSection', {})
            identification = protocol.get('identificationModule', {})
            status = protocol.get('statusModule', {})
            conditions = protocol.get('conditionsModule', {})

            all_studies.append({
                'nct_id': identification.get('nctId'),
                'title': identification.get('briefTitle'),
                'status': status.get('overallStatus'),
                'phase': status.get('phases', ['NA'])[0] if status.get('phases') else 'NA',
                'conditions': '; '.join(conditions.get('conditions', [])),
                'start_date': status.get('startDateStruct', {}).get('date'),
                'completion_date': status.get('completionDateStruct', {}).get('date')
            })

df = pd.DataFrame(all_studies)
df.to_csv('all_trials.tsv', sep='\t', index=False)
print(f"Total trials: {len(df)}")
EOF

# Extract drug trials
python3 << 'EOF'
import json
import pandas as pd

with open('studies.json') as f:
    data = json.load(f)

drug_trials = []
for study in data.get('studies', []):
    protocol = study.get('protocolSection', {})
    arms = protocol.get('armsInterventionsModule', {})
    interventions = arms.get('interventions', [])

    for intervention in interventions:
        if intervention.get('type') == 'DRUG':
            drug_trials.append({
                'nct_id': protocol.get('identificationModule', {}).get('nctId'),
                'drug_name': intervention.get('name'),
                'description': intervention.get('description')
            })

df = pd.DataFrame(drug_trials)
df.to_csv('drug_trials.tsv', sep='\t', index=False)
EOF

# Filter by sponsor type
jq '.studies[] | select(.protocolSection.sponsorCollaboratorsModule.leadSponsor.class == "NIH")' \
  studies.json > nih_sponsored.json
```

## Verification

```bash
# Check JSON structure
cat studies_page1.json | jq 'keys'

# Count studies
cat studies_page1.json | jq '.studies | length'

# Check specific study
curl -s "https://clinicaltrials.gov/api/v2/studies/NCT04368728" | jq '.protocolSection.identificationModule'
```

## Update Schedule

| Data Type | Frequency |
|-----------|-----------|
| Study registrations | Real-time |
| Status updates | Real-time |
| Results posting | Within 1 year of completion |

## Common Issues

- **Pagination**: API returns max 1000 per page; use pageToken
- **Rate limits**: Be respectful; add delays between requests
- **Missing results**: Not all completed trials have posted results
- **Data quality**: Sponsor-reported data; verify critical information
- **Nested JSON**: Complex structure; parse carefully

## API Parameters Reference

| Parameter | Description |
|-----------|-------------|
| query.cond | Search conditions |
| query.intr | Search interventions |
| query.term | General search |
| filter.overallStatus | Trial status |
| filter.phase | Trial phase |
| filter.geo | Geographic location |
| sort | Sort field |
| pageSize | Results per page (max 1000) |
| pageToken | Pagination token |
| format | json or csv |

## Status Values

| Status | Description |
|--------|-------------|
| NOT_YET_RECRUITING | Approved but not started |
| RECRUITING | Currently enrolling |
| ACTIVE_NOT_RECRUITING | Ongoing, closed to enrollment |
| COMPLETED | Finished |
| TERMINATED | Stopped early |
| WITHDRAWN | Never enrolled |
| SUSPENDED | Temporarily halted |

## Related Resources

- [FDA OpenFDA](../fda.openfda/) - Drug approvals
- [DrugBank](../../02.compounds.molecules/2.2.pharmaceuticals/drugbank/) - Drug information
- [PubMed](../8.1.scientific.literature/pubmed/) - Related publications
