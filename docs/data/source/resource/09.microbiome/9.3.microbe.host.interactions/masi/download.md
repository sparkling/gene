---
id: download-masi
title: "MASI Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# MASI Download Instructions

## Quick Start

```bash
# MASI may have limited public availability
# Contact maintainers or check publication supplementary materials
```

## Prerequisites

- Web browser for interface access (if available)
- **Python** for data processing
- Access credentials (if required)
- 50-200MB storage for dataset

## Access Information

MASI availability may be limited. Check:
1. Publication supplementary materials
2. Contact corresponding authors
3. Related databases (HMDB, KEGG)

## Download Methods

### Method 1: Publication Supplementary Data

```bash
# Check supplementary files from MASI publication
# Download supplementary tables (usually Excel/TSV)

# Example pattern (verify actual URLs):
# wget https://[journal]/supplementary/masi_interactions.xlsx
```

### Method 2: Web Interface (if available)

1. Check for web portal access
2. Register if required
3. Query by metabolite or target
4. Export search results

### Method 3: API Access (if available)

```python
import requests

# Example API pattern (verify with actual documentation)
base_url = "https://masi.database.org/api"

def get_interactions_by_metabolite(metabolite_name):
    """Query interactions by metabolite."""
    url = f"{base_url}/interactions"
    params = {"metabolite": metabolite_name}
    response = requests.get(url, params=params)
    return response.json()

def get_interactions_by_target(target_symbol):
    """Query interactions by host target."""
    url = f"{base_url}/interactions"
    params = {"target": target_symbol}
    response = requests.get(url, params=params)
    return response.json()

# Example usage
butyrate_interactions = get_interactions_by_metabolite("butyrate")
gpr43_interactions = get_interactions_by_target("GPR43")
```

### Method 4: Build from Related Sources

If MASI is not directly available, build similar dataset from:

```python
import pandas as pd

# Combine data from multiple sources
# 1. HMDB for metabolite information
# 2. KEGG for pathway data
# 3. Literature mining

# Example: Extract microbiome-signaling from HMDB
hmdb_metabolites = pd.read_csv("hmdb_metabolites.csv")

# Filter for microbially-derived
microbial = hmdb_metabolites[
    hmdb_metabolites['origin'].str.contains('microbial', case=False, na=False)
]

# Extract protein targets
targets = microbial[['name', 'protein_associations']].dropna()
print(f"Found {len(targets)} metabolites with protein associations")
```

## File Inventory

### Expected Data Files (if available)

| File | Size (est.) | Description |
|------|-------------|-------------|
| interactions.tsv | ~2 MB | Main interaction table |
| metabolites.tsv | ~500 KB | Metabolite annotations |
| targets.tsv | ~300 KB | Host target annotations |
| pathways.tsv | ~200 KB | Pathway mappings |

## Post-Download Processing

```python
import pandas as pd

# Load interactions (adjust for actual format)
df = pd.read_csv("interactions.tsv", sep="\t")

# Summary statistics
print(f"Total interactions: {len(df)}")
print(f"Unique metabolites: {df['metabolite_name'].nunique()}")
print(f"Unique targets: {df['target_symbol'].nunique()}")

# Group by metabolite class
if 'metabolite_class' in df.columns:
    print("\nInteractions by metabolite class:")
    print(df['metabolite_class'].value_counts())

# Group by effect type
if 'effect' in df.columns:
    print("\nEffect distribution:")
    print(df['effect'].value_counts())

# Group by pathway
if 'pathway' in df.columns:
    print("\nTop pathways:")
    print(df['pathway'].value_counts().head(10))
```

### Create Network Representation

```python
import pandas as pd
import networkx as nx

# Load data
df = pd.read_csv("interactions.tsv", sep="\t")

# Create bipartite network
G = nx.Graph()

# Add nodes
for metabolite in df['metabolite_name'].unique():
    G.add_node(metabolite, bipartite=0, type='metabolite')

for target in df['target_symbol'].unique():
    G.add_node(target, bipartite=1, type='target')

# Add edges
for _, row in df.iterrows():
    G.add_edge(
        row['metabolite_name'],
        row['target_symbol'],
        effect=row.get('effect', 'unknown'),
        pmid=row.get('pmid', None)
    )

print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Find hub metabolites
metabolite_degree = {n: d for n, d in G.degree() if G.nodes[n]['type'] == 'metabolite'}
top_metabolites = sorted(metabolite_degree.items(), key=lambda x: x[1], reverse=True)[:10]
print("\nTop metabolites by target count:")
for m, d in top_metabolites:
    print(f"  {m}: {d} targets")

# Find hub targets
target_degree = {n: d for n, d in G.degree() if G.nodes[n]['type'] == 'target'}
top_targets = sorted(target_degree.items(), key=lambda x: x[1], reverse=True)[:10]
print("\nTop targets by metabolite count:")
for t, d in top_targets:
    print(f"  {t}: {d} metabolites")

# Export for visualization
nx.write_gexf(G, "masi_network.gexf")
```

### Pathway Enrichment Analysis

```python
import pandas as pd
from scipy import stats

# Load data
df = pd.read_csv("interactions.tsv", sep="\t")

# Group by pathway
pathway_counts = df.groupby('pathway').agg({
    'metabolite_name': 'nunique',
    'target_symbol': 'nunique',
    'interaction_id': 'count'
}).rename(columns={
    'metabolite_name': 'metabolite_count',
    'target_symbol': 'target_count',
    'interaction_id': 'interaction_count'
})

pathway_counts = pathway_counts.sort_values('interaction_count', ascending=False)
print("Pathway summary:")
print(pathway_counts.head(15))
```

## Verification

```bash
# Check file format
head interactions.tsv

# Verify metabolite IDs link to PubChem
python3 << 'EOF'
import pandas as pd
import requests

df = pd.read_csv("interactions.tsv", sep="\t")
sample_cid = df['pubchem_cid'].dropna().iloc[0]

response = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(sample_cid)}/property/MolecularFormula/JSON")
if response.status_code == 200:
    print(f"Verified PubChem CID {sample_cid}")
EOF
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | Check with maintainers |
| Literature integration | Ongoing |

## Common Issues

- **Limited availability**: May require direct contact
- **Data format**: Check publication for format details
- **ID mapping**: Metabolite IDs may vary between sources
- **Incomplete annotations**: Some interactions lack mechanism details
- **Literature access**: PMIDs may require journal subscription

## Alternative Data Sources

```python
# Build similar dataset from public sources
from Bio import Entrez

Entrez.email = "your@email.com"

# Search PubMed for microbiome-signaling literature
handle = Entrez.esearch(
    db="pubmed",
    term="gut microbiome[MeSH] AND signal transduction[MeSH]",
    retmax=100
)
results = Entrez.read(handle)
pmids = results["IdList"]
print(f"Found {len(pmids)} relevant publications")

# Extract metabolite-target pairs from abstracts
# (Requires NLP/text mining approach)
```

## Related Resources

- [VMH](../vmh/_index.md) - Metabolic reconstructions
- [HMDB](../../../../06.nutrition.food/6.4.metabolomics/hmdb/_index.md) - Metabolite database
- [KEGG](../../../../04.pathways.networks/4.1.metabolic.pathways/kegg/_index.md) - Pathway data
- [gutMGene](../../9.1.gut.microbiome/gutmgene/_index.md) - Gene expression
