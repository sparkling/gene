---
id: download-panelapp
title: "PanelApp Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# PanelApp Download Instructions

## Quick Start

```bash
# Download all panels list
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/" -o panelapp_panels.json

# Download specific panel (e.g., Intellectual Disability)
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/245/" -o panel_245.json
```

## Prerequisites

- **curl** or **wget** for downloads
- **JSON parser** (jq recommended)
- Approximately 500MB for all panels

## No Registration Required

PanelApp is openly accessible with no registration needed.

## Download Methods

### Method 1: REST API - Panel List

```bash
# Get all panels
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/" \
  -H "Accept: application/json" \
  -o panelapp_all_panels.json

# Get panels with pagination
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/?page=1" -o panels_page1.json

# Count total panels
cat panelapp_all_panels.json | jq '.count'

# List panel names and IDs
cat panelapp_all_panels.json | jq '.results[] | {id, name}'
```

### Method 2: Download Individual Panels

```bash
# Download panel by ID
PANEL_ID=245  # Intellectual disability
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/${PANEL_ID}/" \
  -H "Accept: application/json" \
  -o "panel_${PANEL_ID}.json"

# Download panel genes only
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/${PANEL_ID}/genes/" \
  -o "panel_${PANEL_ID}_genes.json"

# Download panel in TSV format
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/${PANEL_ID}/?format=tsv" \
  -o "panel_${PANEL_ID}.tsv"
```

### Method 3: Download All Panels Programmatically

```bash
# Fetch all panel IDs and download each
python3 << 'EOF'
import requests
import json
import time

base_url = "https://panelapp.genomicsengland.co.uk/api/v1"

# Get all panels
response = requests.get(f"{base_url}/panels/")
panels = response.json()

print(f"Total panels: {panels['count']}")

# Download each panel
for panel in panels['results'][:10]:  # Limit to 10 for demo
    panel_id = panel['id']
    panel_name = panel['name'].replace('/', '_')

    # Get full panel details
    detail_response = requests.get(f"{base_url}/panels/{panel_id}/")

    with open(f"panel_{panel_id}_{panel_name[:30]}.json", 'w') as f:
        json.dump(detail_response.json(), f, indent=2)

    print(f"Downloaded: {panel_id} - {panel['name']}")
    time.sleep(0.5)  # Rate limiting
EOF
```

### Method 4: PanelApp Australia

```bash
# Australian mirror
BASE_URL="https://panelapp.agha.umccr.org/api/v1"

# Get panels
curl "${BASE_URL}/panels/" -o panelapp_au_panels.json

# Get specific panel
curl "${BASE_URL}/panels/123/" -o panelapp_au_panel_123.json
```

### Method 5: Search and Filter

```bash
# Search panels by name
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/?name=intellectual" \
  -o panels_intellectual.json

# Search by gene
curl "https://panelapp.genomicsengland.co.uk/api/v1/genes/?gene_symbol=MECP2" \
  -o gene_mecp2_panels.json

# Filter by disease group
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/?disease_group=Neurology%20and%20neurodevelopmental%20disorders" \
  -o panels_neurology.json

# Get green genes only (high evidence)
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/245/genes/?confidence_level=3" \
  -o panel_245_green_genes.json
```

### Method 6: Export Formats

```bash
# TSV export
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/245/?format=tsv" \
  -o panel_245.tsv

# BED format (gene coordinates)
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/245/genes/?format=bed" \
  -o panel_245.bed

# Download signed-off versions
curl "https://panelapp.genomicsengland.co.uk/api/v1/panels/signedoff/" \
  -o signed_off_panels.json
```

## File Inventory

### API Endpoints

| Endpoint | Description |
|----------|-------------|
| /panels/ | List all panels |
| /panels/{id}/ | Panel details |
| /panels/{id}/genes/ | Panel genes |
| /genes/ | Search genes |
| /panels/signedoff/ | Signed-off versions |

### Panel JSON Structure

| Field | Description |
|-------|-------------|
| id | Panel identifier |
| name | Panel name |
| disease_group | Disease category |
| disease_sub_group | Subcategory |
| version | Current version |
| genes | List of gene entries |
| strs | Short tandem repeats |
| regions | Genomic regions |

### Gene Entry Fields

| Field | Description |
|-------|-------------|
| gene_symbol | HGNC symbol |
| gene_name | Full gene name |
| hgnc_id | HGNC identifier |
| confidence_level | 1-3 (Red-Amber-Green) |
| mode_of_inheritance | AR, AD, XL, etc. |
| phenotypes | Associated phenotypes |
| publications | Supporting PMIDs |
| evidence | Evidence sources |

## Post-Download Processing

```bash
# Extract green genes from panel
python3 << 'EOF'
import json
import pandas as pd

with open('panel_245.json') as f:
    panel = json.load(f)

genes = []
for gene in panel['genes']:
    if gene['confidence_level'] == 3:  # Green
        genes.append({
            'symbol': gene['gene_data']['gene_symbol'],
            'hgnc_id': gene['gene_data']['hgnc_id'],
            'moi': gene.get('mode_of_inheritance', ''),
            'phenotypes': ';'.join(gene.get('phenotypes', [])),
        })

df = pd.DataFrame(genes)
df.to_csv('panel_245_green_genes.tsv', sep='\t', index=False)
print(f"Green genes: {len(df)}")
EOF

# Create BED file for panel
python3 << 'EOF'
import json
import requests

# Get gene coordinates from Ensembl
def get_gene_coords(symbol):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
    response = requests.get(url, headers={"Content-Type": "application/json"})
    if response.status_code == 200:
        data = response.json()
        return data.get('seq_region_name'), data.get('start'), data.get('end')
    return None, None, None

with open('panel_245.json') as f:
    panel = json.load(f)

with open('panel_245_coords.bed', 'w') as bed:
    for gene in panel['genes'][:10]:  # Limit for demo
        symbol = gene['gene_data']['gene_symbol']
        chrom, start, end = get_gene_coords(symbol)
        if chrom and start and end:
            confidence = gene['confidence_level']
            bed.write(f"chr{chrom}\t{start}\t{end}\t{symbol}\t{confidence}\n")
EOF

# Aggregate all panels for a gene
python3 << 'EOF'
import json
import requests

gene = "MECP2"
response = requests.get(
    f"https://panelapp.genomicsengland.co.uk/api/v1/genes/?gene_symbol={gene}"
)
data = response.json()

print(f"Gene {gene} appears in {data['count']} panels:")
for entry in data['results']:
    print(f"  - {entry['panel']['name']} (confidence: {entry['confidence_level']})")
EOF

# Compare panels
python3 << 'EOF'
import json

# Load two panels
with open('panel_245.json') as f:
    panel1 = json.load(f)
with open('panel_285.json') as f:
    panel2 = json.load(f)

genes1 = {g['gene_data']['gene_symbol'] for g in panel1['genes']}
genes2 = {g['gene_data']['gene_symbol'] for g in panel2['genes']}

overlap = genes1.intersection(genes2)
unique1 = genes1 - genes2
unique2 = genes2 - genes1

print(f"Panel 1 genes: {len(genes1)}")
print(f"Panel 2 genes: {len(genes2)}")
print(f"Overlap: {len(overlap)}")
print(f"Unique to Panel 1: {len(unique1)}")
print(f"Unique to Panel 2: {len(unique2)}")
EOF
```

## Verification

```bash
# Check panel download
cat panel_245.json | jq '.name, .version, .genes | length'

# List gene confidence levels
cat panel_245.json | jq '.genes[].confidence_level' | sort | uniq -c

# Verify TSV format
head panel_245.tsv

# Check signed-off version
cat panel_245.json | jq '.version_created'
```

## Update Schedule

| Update Type | Frequency |
|-------------|-----------|
| Expert reviews | Continuous |
| New panels | As needed |
| Version sign-off | Periodic |

## Common Issues

- **Panel versions**: Panels have multiple versions; use signed-off for clinical
- **Confidence changes**: Gene evidence level can change over time
- **Gene symbols**: Use HGNC IDs for stable identifiers
- **API pagination**: Large result sets are paginated
- **Rate limiting**: Add delays between bulk requests

## Confidence Levels

| Level | Color | Clinical Use |
|-------|-------|--------------|
| 3 | Green | Diagnostic grade |
| 2 | Amber | Moderate evidence |
| 1 | Red | Low evidence |

## Mode of Inheritance Codes

| Code | Meaning |
|------|---------|
| MONOALLELIC | Autosomal dominant |
| BIALLELIC | Autosomal recessive |
| X-LINKED | X-linked |
| BOTH | Multiple inheritance modes |
| MITOCHONDRIAL | Mitochondrial |

## Integration Examples

```bash
# Map PanelApp to HPO
python3 << 'EOF'
import json
import pandas as pd

with open('panel_245.json') as f:
    panel = json.load(f)

# Extract HPO terms
gene_hpo = []
for gene in panel['genes']:
    symbol = gene['gene_data']['gene_symbol']
    for pheno in gene.get('phenotypes', []):
        gene_hpo.append({
            'gene': symbol,
            'phenotype': pheno
        })

df = pd.DataFrame(gene_hpo)
df.to_csv('panel_gene_phenotypes.tsv', sep='\t', index=False)
EOF

# Cross-reference with OMIM
python3 << 'EOF'
import json

with open('panel_245.json') as f:
    panel = json.load(f)

omim_genes = []
for gene in panel['genes']:
    omim = gene['gene_data'].get('omim_gene', [])
    if omim:
        omim_genes.append({
            'symbol': gene['gene_data']['gene_symbol'],
            'omim': omim
        })

print(f"Genes with OMIM: {len(omim_genes)}")
EOF
```

## Related Resources

- [Orphanet](../orphanet/) - Rare disease database
- [DECIPHER](../decipher/) - CNV database
- [OMIM](../../3.2.phenotype.databases/omim/) - Mendelian genetics
- [ClinGen](../../../../01.genetics.genomics/1.1.variant.repositories/clinvar/) - Gene curation
