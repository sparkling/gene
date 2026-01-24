---
id: download-intact
title: "IntAct Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# IntAct - Download Documentation

## Overview

IntAct provides molecular interaction data via PSICQUIC REST API and FTP bulk downloads. As an IMEx consortium member, data follows PSI-MI standards.

## PSICQUIC API

### Base URL

```
https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search
```

### Search by Identifier

```bash
# Search by UniProt AC
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/identifier:P04637?format=tab25"

# Search by gene name (alias)
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/alias:TP53?format=tab27"

# Search by IntAct ID
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/interaction_id:EBI-77734?format=tab25"
```

### Filter Queries

```bash
# Human-human interactions
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/taxidA:9606%20AND%20taxidB:9606?format=tab25"

# By detection method (two-hybrid = MI:0018)
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/detmethod:%22MI:0018%22?format=tab25"

# By interaction type (direct = MI:0407)
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/type:%22MI:0407%22?format=tab25"

# By publication
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/pubid:1535557?format=tab25"

# Combined query
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/identifier:P04637%20AND%20taxidB:9606?format=tab25"
```

### Pagination

```bash
# First 100 results
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/taxidA:9606?firstResult=0&maxResults=100&format=tab25"

# Next 100 results
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/taxidA:9606?firstResult=100&maxResults=100&format=tab25"
```

### Output Formats

```bash
# PSI-MITAB 2.5 (15 columns)
curl "...?format=tab25"

# PSI-MITAB 2.7 (42 columns)
curl "...?format=tab27"

# XML PSI-MI 2.5
curl "...?format=xml25"

# JSON
curl "...?format=json"

# Count only
curl "...?format=count"
```

## FTP Bulk Downloads

### FTP Server

```
ftp://ftp.ebi.ac.uk/pub/databases/intact/current/
```

### Full Database

```bash
# PSI-MITAB 2.7 (all interactions)
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt

# PSI-MITAB with negative interactions
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact_negative.txt

# PSI-MI XML
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/pmidMIF25.zip

# XGMML (for Cytoscape)
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/various/intact.xgmml
```

### Species-Specific

```bash
# Human interactions
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.txt

# Mouse interactions
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/mouse.txt

# Yeast interactions
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/yeast.txt
```

### Complex Portal

```bash
# Curated complexes
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/homo_sapiens.tsv
```

## PSICQUIC Query Syntax

### Field Names

| Field | Description | Example |
|-------|-------------|---------|
| identifier | Any ID | identifier:P04637 |
| id | Same as identifier | id:P04637 |
| alias | Gene names/synonyms | alias:TP53 |
| pubid | PubMed ID | pubid:1535557 |
| pubauth | Publication author | pubauth:Smith |
| taxidA/taxidB | Taxonomy ID | taxidA:9606 |
| species | Organism name | species:human |
| type | Interaction type MI | type:"MI:0407" |
| detmethod | Detection method MI | detmethod:"MI:0018" |
| interaction_id | IntAct ID | interaction_id:EBI-77734 |

### Boolean Operators

```
identifier:P04637 AND identifier:P38398
identifier:P04637 OR identifier:P38398
identifier:P04637 NOT taxidB:10090
```

### Wildcards

```
alias:TP5*
identifier:P0463?
```

## Python Examples

### PSICQUIC API Access

```python
import requests
import pandas as pd
from io import StringIO

class IntActClient:
    def __init__(self):
        self.base_url = "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search"

    def query(self, query_str, format="tab27", max_results=None):
        """Execute PSICQUIC query."""
        params = {"format": format}
        if max_results:
            params["maxResults"] = max_results

        url = f"{self.base_url}/query/{query_str}"
        response = requests.get(url, params=params)

        return response.text

    def get_interactions(self, protein_id, species=9606):
        """Get interactions for a protein."""
        query = f"identifier:{protein_id} AND taxidA:{species} AND taxidB:{species}"
        data = self.query(query)

        # Parse MITAB
        columns = [
            'ID_A', 'ID_B', 'Alt_ID_A', 'Alt_ID_B', 'Alias_A', 'Alias_B',
            'Detection_Method', 'Author', 'Publication', 'TaxID_A', 'TaxID_B',
            'Type', 'Source', 'Interaction_ID', 'Confidence'
        ]

        df = pd.read_csv(StringIO(data), sep='\t', header=None,
                        names=columns[:15] if data.count('\t') < 20 else None)

        return df

    def count(self, query_str):
        """Count interactions matching query."""
        response = self.query(query_str, format="count")
        return int(response.strip())

# Example usage
client = IntActClient()
interactions = client.get_interactions("P04637")
count = client.count("taxidA:9606 AND taxidB:9606")
print(f"Human-human interactions: {count}")
```

### Parse MITAB File

```python
def parse_mitab(filepath):
    """Parse PSI-MITAB file."""
    columns_27 = [
        'ID_A', 'ID_B', 'Alt_ID_A', 'Alt_ID_B', 'Alias_A', 'Alias_B',
        'Detection_Method', 'Author', 'Publication', 'TaxID_A', 'TaxID_B',
        'Type', 'Source', 'Interaction_ID', 'Confidence',
        'Expansion', 'Bio_Role_A', 'Bio_Role_B', 'Exp_Role_A', 'Exp_Role_B',
        'Type_A', 'Type_B', 'Xref_A', 'Xref_B', 'Xref_Int',
        'Annotation_A', 'Annotation_B', 'Annotation_Int',
        'Taxid_Host', 'Params', 'Creation_Date', 'Update_Date',
        'Checksum_A', 'Checksum_B', 'Checksum_Int', 'Negative',
        'Feature_A', 'Feature_B', 'Stoichiometry_A', 'Stoichiometry_B',
        'Participant_Ident_A', 'Participant_Ident_B'
    ]

    df = pd.read_csv(filepath, sep='\t', header=None, names=columns_27,
                     comment='#', low_memory=False)

    return df

def extract_uniprot_ids(mitab_df):
    """Extract UniProt IDs from MITAB data."""
    def parse_id(id_str):
        if pd.isna(id_str):
            return None
        for part in str(id_str).split('|'):
            if 'uniprotkb:' in part:
                return part.split('uniprotkb:')[1].split('(')[0]
        return None

    mitab_df['UniProt_A'] = mitab_df['ID_A'].apply(parse_id)
    mitab_df['UniProt_B'] = mitab_df['ID_B'].apply(parse_id)

    return mitab_df
```

### Build Network

```python
import networkx as nx

def build_network(mitab_df, score_threshold=0.4):
    """Build NetworkX graph from IntAct data."""
    G = nx.Graph()

    for _, row in mitab_df.iterrows():
        id_a = row.get('UniProt_A') or row['ID_A']
        id_b = row.get('UniProt_B') or row['ID_B']

        if not id_a or not id_b:
            continue

        # Parse MI score if present
        score = 0.5  # default
        conf = str(row.get('Confidence', ''))
        if 'intact-miscore:' in conf:
            try:
                score = float(conf.split('intact-miscore:')[1].split('|')[0])
            except:
                pass

        if score >= score_threshold:
            if G.has_edge(id_a, id_b):
                G[id_a][id_b]['weight'] += 1
            else:
                G.add_edge(id_a, id_b, weight=1, score=score)

    return G
```

## Complex Portal

### Web Interface

```
https://www.ebi.ac.uk/complexportal/
```

### API

```bash
# Get complex by ID
curl "https://www.ebi.ac.uk/intact/complex-ws/export/CPX-1?format=json"

# Search complexes
curl "https://www.ebi.ac.uk/intact/complex-ws/search?query=p53&species=9606"
```

## IMEx Consortium Federated Query

Query all IMEx databases simultaneously:

```bash
# Query all IMEx members
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/imex/webservices/current/search/query/identifier:P04637?format=tab25"
```

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| PSICQUIC API | No strict limit |
| Recommended | 1 request/second |
| FTP | Unlimited |

## Update Frequency

| Release | Frequency |
|---------|-----------|
| Monthly release | Monthly |
| Incremental | Weekly |

## See Also

- [Schema Documentation](./schema.md)
- [IntAct Documentation](https://www.ebi.ac.uk/intact/documentation)
