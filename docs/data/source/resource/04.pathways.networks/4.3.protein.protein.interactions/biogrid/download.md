---
id: download-biogrid
title: "BioGRID Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# BioGRID - Download Documentation

## Overview

BioGRID provides free access via REST API (requires API key) and bulk downloads. All data is available under the MIT license.

## REST API

### Registration

Get free API key at: https://webservice.thebiogrid.org/

### Base URL

```
https://webservice.thebiogrid.org
```

### Search by Gene

```bash
# Search by gene symbol
curl "https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=TP53&taxId=9606&accesskey=YOUR_KEY&format=json"

# Multiple genes
curl "https://webservice.thebiogrid.org/interactions/?geneList=TP53|MDM2|CDKN1A&taxId=9606&accesskey=YOUR_KEY&format=tab2"

# Search any field
curl "https://webservice.thebiogrid.org/interactions/?searchNames=true&searchSynonyms=true&geneList=p53&accesskey=YOUR_KEY"
```

### Search by BioGRID ID

```bash
curl "https://webservice.thebiogrid.org/interactions/?interactorList=112315&accesskey=YOUR_KEY&format=json"
```

### Search by Publication

```bash
curl "https://webservice.thebiogrid.org/interactions/?pubmedList=1535557&accesskey=YOUR_KEY&format=json"
```

### Filter Parameters

```bash
# Physical interactions only
curl "https://webservice.thebiogrid.org/interactions/?geneList=TP53&evidenceList=physical&taxId=9606&accesskey=YOUR_KEY"

# Specific experimental system
curl "https://webservice.thebiogrid.org/interactions/?geneList=TP53&experimentalSystemList=Two-hybrid|Affinity Capture-MS&accesskey=YOUR_KEY"

# Exclude self-interactions
curl "https://webservice.thebiogrid.org/interactions/?geneList=TP53&selfInteractionsExcluded=true&accesskey=YOUR_KEY"

# High throughput only
curl "https://webservice.thebiogrid.org/interactions/?geneList=TP53&throughputTag=high&accesskey=YOUR_KEY"
```

### Pagination

```bash
# Start and max results
curl "https://webservice.thebiogrid.org/interactions/?geneList=TP53&start=0&max=1000&accesskey=YOUR_KEY"
```

## API Parameters

| Parameter | Description | Values |
|-----------|-------------|--------|
| geneList | Gene symbols (pipe-separated) | TP53\|MDM2 |
| interactorList | BioGRID IDs | 112315 |
| pubmedList | PubMed IDs | 1535557 |
| taxId | Organism taxonomy ID | 9606 |
| searchNames | Search official names | true/false |
| searchSynonyms | Search synonyms | true/false |
| evidenceList | Interaction type | physical\|genetic |
| experimentalSystemList | Detection methods | Two-hybrid |
| selfInteractionsExcluded | Exclude self | true/false |
| throughputTag | Throughput level | high/low |
| start | Pagination start | 0 |
| max | Max results | 10000 |
| format | Output format | json/tab2/tab3 |

## Bulk Downloads

### Download Portal

```
https://downloads.thebiogrid.org/
```

### Current Release

```bash
# All organisms (TAB 2.0)
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.xxx/BIOGRID-ALL-4.4.xxx.tab2.zip

# All organisms (TAB 3.0)
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.xxx/BIOGRID-ALL-4.4.xxx.tab3.zip

# All organisms (PSI-MI TAB 2.5)
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.xxx/BIOGRID-ALL-4.4.xxx.mitab.zip
```

### Single Organism

```bash
# Human only
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.xxx/BIOGRID-ORGANISM-Homo_sapiens-4.4.xxx.tab2.zip

# Yeast only
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.xxx/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.xxx.tab2.zip
```

### Identifiers Only

```bash
# ID mapping file
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.xxx/BIOGRID-IDENTIFIERS-4.4.xxx.tab.zip
```

## Output Formats

| Format | Parameter | Description |
|--------|-----------|-------------|
| TAB 2.0 | tab2 | 26 columns |
| TAB 3.0 | tab3 | 35 columns with ontology |
| PSI-MITAB | mitab | 15 columns standard |
| JSON | json | Structured format |

## Python Examples

### API Access

```python
import requests
import time

class BioGRIDClient:
    def __init__(self, api_key):
        self.base_url = "https://webservice.thebiogrid.org/interactions/"
        self.api_key = api_key

    def get_interactions(self, gene, organism=9606, evidence_type=None):
        """Get interactions for a gene."""
        params = {
            "geneList": gene,
            "taxId": organism,
            "searchNames": "true",
            "accesskey": self.api_key,
            "format": "json",
            "max": 10000
        }

        if evidence_type:
            params["evidenceList"] = evidence_type

        response = requests.get(self.base_url, params=params)
        return response.json()

    def get_interactors(self, gene, organism=9606):
        """Get list of interacting proteins."""
        data = self.get_interactions(gene, organism, "physical")

        interactors = set()
        for interaction in data.values():
            if interaction["OFFICIAL_SYMBOL_A"] != gene:
                interactors.add(interaction["OFFICIAL_SYMBOL_A"])
            if interaction["OFFICIAL_SYMBOL_B"] != gene:
                interactors.add(interaction["OFFICIAL_SYMBOL_B"])

        return list(interactors)
```

### Bulk Data Processing

```python
import pandas as pd
import zipfile
import io

def load_biogrid_tab2(filepath):
    """Load BioGRID TAB 2.0 file."""
    columns = [
        'INTERACTION_ID', 'ENTREZ_GENE_A', 'ENTREZ_GENE_B',
        'BIOGRID_A', 'BIOGRID_B', 'SYSTEMATIC_NAME_A', 'SYSTEMATIC_NAME_B',
        'SYMBOL_A', 'SYMBOL_B', 'SYNONYMS_A', 'SYNONYMS_B',
        'EXPERIMENTAL_SYSTEM', 'EXPERIMENTAL_SYSTEM_TYPE',
        'AUTHOR', 'PUBMED_ID', 'ORGANISM_A', 'ORGANISM_B',
        'THROUGHPUT', 'SCORE', 'MODIFICATION', 'QUALIFICATIONS',
        'TAGS', 'SOURCE_DATABASE', 'SWISSPROT_A', 'TREMBL_A', 'SWISSPROT_B'
    ]

    with zipfile.ZipFile(filepath) as z:
        for name in z.namelist():
            if name.endswith('.txt'):
                with z.open(name) as f:
                    df = pd.read_csv(f, sep='\t', names=columns, skiprows=1)
                    return df

def filter_human_physical(df):
    """Filter for human physical interactions."""
    return df[
        (df['ORGANISM_A'] == 9606) &
        (df['ORGANISM_B'] == 9606) &
        (df['EXPERIMENTAL_SYSTEM_TYPE'] == 'physical')
    ]
```

### Build Network

```python
import networkx as nx

def build_ppi_network(interactions_df):
    """Build NetworkX graph from interactions."""
    G = nx.Graph()

    for _, row in interactions_df.iterrows():
        gene_a = row['SYMBOL_A']
        gene_b = row['SYMBOL_B']

        if gene_a and gene_b and gene_a != '-' and gene_b != '-':
            if G.has_edge(gene_a, gene_b):
                G[gene_a][gene_b]['weight'] += 1
                G[gene_a][gene_b]['pmids'].add(row['PUBMED_ID'])
            else:
                G.add_edge(gene_a, gene_b,
                          weight=1,
                          pmids={row['PUBMED_ID']},
                          methods={row['EXPERIMENTAL_SYSTEM']})

    return G
```

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| API with key | No strict limit |
| Recommended | 1 request/second |
| Bulk download | Unlimited |

## Data Size Estimates

| Organism | Interactions |
|----------|-------------|
| Human | ~700,000 |
| Yeast | ~500,000 |
| All organisms | ~2,200,000 |
| File size (all) | ~500 MB zipped |

---

## Dataset Versions

### Current Release: BioGRID 4.4.xxx

| Property | Value |
|----------|-------|
| Version | 4.4.230 (example) |
| Release Date | 2024-01-15 |
| Total Size | ~500 MB (zipped) |
| Total Interactions | ~2.2 million |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| BIOGRID-ALL-*.tab2.zip | ~200 MB | ~2.2M | All organisms TAB 2.0 |
| BIOGRID-ALL-*.tab3.zip | ~250 MB | ~2.2M | All organisms TAB 3.0 |
| BIOGRID-ALL-*.mitab.zip | ~150 MB | ~2.2M | PSI-MI TAB 2.5 |
| BIOGRID-ORGANISM-Homo_sapiens-*.tab2.zip | ~50 MB | ~700K | Human only |

### Previous Versions

| Version | Release | Interactions | Status |
|---------|---------|--------------|--------|
| 4.4.229 | 2024-01-01 | ~2.2M | Archived |
| 4.4.228 | 2023-12-01 | ~2.2M | Archived |
| 4.4.227 | 2023-11-01 | ~2.1M | Archived |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| Base URL | `https://webservice.thebiogrid.org` |
| Authentication | API key required (free) |
| Rate Limit | 1 request/second recommended |
| Response Format | JSON, TAB2, TAB3 |

### API Endpoints

| Operation | Endpoint | Parameters |
|-----------|----------|------------|
| Interactions | `/interactions/` | `geneList`, `taxId`, `evidenceList` |
| By Gene | `/interactions/` | `geneList=TP53&taxId=9606` |
| By PubMed | `/interactions/` | `pubmedList=1535557` |
| By BioGRID ID | `/interactions/` | `interactorList=112315` |

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| accesskey | Your API key (required) |
| format | Output format (json, tab2, tab3) |

---

## Update Frequency

| Release Type | Frequency |
|--------------|-----------|
| Major release | Monthly |
| Updates | Continuous |

## See Also

- [Schema Documentation](./schema.md)
- [BioGRID REST API](https://wiki.thebiogrid.org/doku.php/biogridrest)
