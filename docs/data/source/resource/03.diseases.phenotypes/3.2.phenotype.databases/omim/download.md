---
id: download-omim
title: "OMIM Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# OMIM Download Instructions

## Quick Start

```bash
# After registration, download via FTP
wget --user=your_email --password=your_key \
  ftp://ftp.omim.org/OMIM/mim2gene.txt
```

## Prerequisites

- **OMIM API key** (required - free registration for academic use)
- **wget** or **curl** for downloads
- Approximately 1GB disk space

## Registration Process

### Step 1: Apply for Access

1. Navigate to https://www.omim.org/downloads
2. Click "Request Downloads"
3. Complete application form (academic affiliation required)
4. Agree to Terms of Use

### Step 2: Receive API Key

1. Applications reviewed within 1-2 business days
2. Receive API key via email
3. API key works for both downloads and API

## Access Restrictions

| Use Type | Access |
|----------|--------|
| Academic research | Free with registration |
| Non-commercial | Free with registration |
| Commercial | Requires license |
| API access | Rate limited |

## Download Methods

### Method 1: FTP Downloads

```bash
# Set credentials
OMIM_EMAIL="your@email.com"
OMIM_KEY="your_api_key"

# Gene-MIM mapping
wget --user=${OMIM_EMAIL} --password=${OMIM_KEY} \
  ftp://ftp.omim.org/OMIM/mim2gene.txt

# MIM titles
wget --user=${OMIM_EMAIL} --password=${OMIM_KEY} \
  ftp://ftp.omim.org/OMIM/mimTitles.txt

# Gene map
wget --user=${OMIM_EMAIL} --password=${OMIM_KEY} \
  ftp://ftp.omim.org/OMIM/genemap2.txt

# Phenotype series
wget --user=${OMIM_EMAIL} --password=${OMIM_KEY} \
  ftp://ftp.omim.org/OMIM/phenotypicSeries.txt

# OMIM text (full entries)
wget --user=${OMIM_EMAIL} --password=${OMIM_KEY} \
  ftp://ftp.omim.org/OMIM/omim.txt.gz
```

### Method 2: API Access

```bash
# Set API key
OMIM_KEY="your_api_key"

# Get entry by MIM number
curl "https://api.omim.org/api/entry?mimNumber=100100&include=text&format=json&apiKey=${OMIM_KEY}" \
  -o entry_100100.json

# Search for entries
curl "https://api.omim.org/api/entry/search?search=breast+cancer&include=all&limit=20&format=json&apiKey=${OMIM_KEY}" \
  -o breast_cancer_search.json

# Get gene map
curl "https://api.omim.org/api/geneMap/search?phenotypeMapExists=true&limit=100&format=json&apiKey=${OMIM_KEY}" \
  -o gene_map.json

# Get allelic variants
curl "https://api.omim.org/api/allelicVariant?mimNumber=100100&format=json&apiKey=${OMIM_KEY}" \
  -o variants_100100.json
```

### Method 3: Batch API Downloads

```bash
# Download multiple entries
MIM_NUMBERS="100100,100200,100300"
curl "https://api.omim.org/api/entry?mimNumber=${MIM_NUMBERS}&include=text&format=json&apiKey=${OMIM_KEY}" \
  -o batch_entries.json

# Iterate through all gene map entries
python3 << 'EOF'
import requests
import json
import time

api_key = "your_api_key"
base_url = "https://api.omim.org/api"

# Get gene map with pagination
offset = 0
limit = 100
all_entries = []

while True:
    url = f"{base_url}/geneMap/search?phenotypeMapExists=true&start={offset}&limit={limit}&format=json&apiKey={api_key}"
    response = requests.get(url)
    data = response.json()

    entries = data['omim']['searchResponse']['geneMapList']
    if not entries:
        break

    all_entries.extend(entries)
    offset += limit
    time.sleep(0.5)  # Rate limiting

    if offset > data['omim']['searchResponse']['totalResults']:
        break

with open('all_gene_map.json', 'w') as f:
    json.dump(all_entries, f)
print(f"Downloaded {len(all_entries)} entries")
EOF
```

### Method 4: Specific Data Types

```bash
# Clinical synopses
curl "https://api.omim.org/api/clinicalSynopsis?mimNumber=100100&format=json&apiKey=${OMIM_KEY}" \
  -o clinical_synopsis.json

# Phenotypic series
curl "https://api.omim.org/api/phenotypicSeries?phenotypicSeriesNumber=PS100100&format=json&apiKey=${OMIM_KEY}" \
  -o phenotypic_series.json

# External links
curl "https://api.omim.org/api/entry?mimNumber=100100&include=externalLinks&format=json&apiKey=${OMIM_KEY}" \
  -o external_links.json
```

## File Inventory

### FTP Files

| File | Size | Description |
|------|------|-------------|
| omim.txt.gz | ~200 MB | Complete OMIM text |
| mim2gene.txt | ~1 MB | MIM to gene mapping |
| mimTitles.txt | ~3 MB | Entry titles |
| genemap2.txt | ~5 MB | Gene-phenotype map |
| phenotypicSeries.txt | ~500 KB | Phenotype groupings |
| morbidmap.txt | ~2 MB | Disease-gene associations |

### Derived Data

| File | Size | Description |
|------|------|-------------|
| allelicVariants.txt | ~5 MB | All variants |
| clinicalSynopses.txt | ~10 MB | Clinical descriptions |

## Post-Download Processing

```bash
# Parse gene map
python3 << 'EOF'
import pandas as pd

# Read genemap2
genemap = pd.read_csv('genemap2.txt', sep='\t', skiprows=3,
                       names=['chromosome', 'genomic_start', 'genomic_end',
                              'cyto_location', 'computed_cyto', 'mim_number',
                              'gene_symbols', 'gene_name', 'approved_symbol',
                              'entrez_id', 'ensembl_id', 'comments', 'phenotypes',
                              'mouse_info'])

# Extract disease genes
disease_genes = genemap[genemap['phenotypes'].notna()]
print(f"Disease-associated genes: {len(disease_genes)}")

# Save cleaned version
disease_genes[['mim_number', 'approved_symbol', 'phenotypes']].to_csv(
    'disease_genes.tsv', sep='\t', index=False)
EOF

# Parse mim2gene
awk -F'\t' 'NR>1 && $2=="gene" {print $1"\t"$4"\t"$5}' mim2gene.txt > gene_mim.tsv

# Extract phenotypes from genemap2
python3 << 'EOF'
import re
import pandas as pd

genemap = pd.read_csv('genemap2.txt', sep='\t', skiprows=3)

phenotype_data = []
for _, row in genemap.iterrows():
    if pd.isna(row.iloc[12]):
        continue
    gene_mim = row.iloc[5]
    gene_symbol = row.iloc[8]

    # Parse phenotypes (format: "Name, MIM, inheritance")
    phenos = str(row.iloc[12]).split(';')
    for pheno in phenos:
        match = re.search(r'(.+?),\s*(\d{6})', pheno)
        if match:
            phenotype_data.append({
                'gene_mim': gene_mim,
                'gene_symbol': gene_symbol,
                'phenotype_name': match.group(1).strip(),
                'phenotype_mim': match.group(2)
            })

df = pd.DataFrame(phenotype_data)
df.to_csv('gene_phenotype_map.tsv', sep='\t', index=False)
EOF
```

## Verification

```bash
# Check file format
head -20 genemap2.txt

# Count entries
wc -l mim2gene.txt

# Verify API response
cat entry_100100.json | jq '.omim.entryList[0].entry.titles'
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| OMIM Jan 2026 | 2026-01-25 | ~300 MB | Current |
| Daily updates | Continuous | Varies | Active |

### Version Notes

OMIM current database statistics:
- 27,500+ entries (genes and phenotypes)
- 17,300+ gene entries
- 8,800+ phenotype entries with known molecular basis
- 4,700+ phenotypic series
- Daily curation by medical geneticists

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://api.omim.org/api` |
| Rate Limit | 4 req/sec |
| Auth Required | Yes (API key required) |
| Documentation | https://omim.org/help/api |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | Daily |
| FTP files | Weekly |
| New entries | Continuous |

## Common Issues

- **API rate limits**: 4 requests/second; use delays in batch scripts
- **Missing entries**: Not all MIM numbers exist; handle 404s gracefully
- **Text formatting**: Full text has special characters; use proper encoding
- **Phenotype parsing**: Complex format; use regex carefully
- **Gene symbols**: May be outdated; cross-reference with HGNC

## MIM Number Types

| Prefix | Type |
|--------|------|
| * | Gene |
| + | Gene with phenotype |
| # | Phenotype (molecular basis known) |
| % | Phenotype (molecular basis unknown) |
| None | Phenotype or other |
| ^ | Moved/removed entry |

## API Parameters

| Parameter | Description |
|-----------|-------------|
| mimNumber | Entry ID(s), comma-separated |
| include | text, clinicalSynopsis, externalLinks, geneMap, etc. |
| format | json, xml |
| search | Search query |
| limit | Results per page |
| start | Pagination offset |

## Related Resources

- [ClinVar](../../01.genetics.genomics/1.1.variant.repositories/clinvar/) - Clinical variants
- [HPO](../hpo/) - Phenotype ontology
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/) - Rare diseases
