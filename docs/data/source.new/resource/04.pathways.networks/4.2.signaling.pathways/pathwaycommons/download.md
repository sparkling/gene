---
id: download-pathwaycommons
title: "Pathway Commons Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# Pathway Commons Download Instructions

## Quick Start

```bash
# Download all data as SIF (Simple Interaction Format)
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.sif.gz
gunzip PathwayCommons12.All.hgnc.sif.gz

# Download all data as BioPAX
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.BIOPAX.owl.gz
gunzip PathwayCommons12.All.BIOPAX.owl.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- ~5-15GB disk space for full downloads
- For BioPAX processing: Paxtools library (Java) or RDFlib (Python)
- For SIF analysis: Cytoscape, NetworkX, or igraph

## No Registration Required

All Pathway Commons data is freely available without registration.

## Download Methods

### Method 1: HTTP Bulk Downloads (Recommended)

```bash
# Create download directory
mkdir -p pathwaycommons && cd pathwaycommons

# Download SIF format (simplest - for network analysis)
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.sif.gz
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.uniprot.sif.gz

# Download Extended SIF (includes pubmed, mediators)
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.txt.gz

# Download Full BioPAX (complete pathway data)
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.BIOPAX.owl.gz

# Download GSEA Gene Sets
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.gmt.gz
```

### Method 2: Individual Data Sources

```bash
# Download specific source (Reactome only)
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reactome.BIOPAX.owl.gz
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.reactome.hgnc.sif.gz

# Download KEGG pathways
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.kegg.BIOPAX.owl.gz
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.kegg.hgnc.sif.gz

# Download WikiPathways
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.wp.BIOPAX.owl.gz

# Download BioGRID interactions
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.biogrid.BIOPAX.owl.gz

# Download IntAct interactions
wget https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.intact.BIOPAX.owl.gz
```

### Method 3: PC2 API Download

```bash
# Get specific pathway as BioPAX
curl -o apoptosis.owl "https://www.pathwaycommons.org/pc2/get?uri=http://identifiers.org/reactome/R-HSA-109581&format=BIOPAX"

# Get pathway as SIF
curl -o apoptosis.sif "https://www.pathwaycommons.org/pc2/get?uri=http://identifiers.org/reactome/R-HSA-109581&format=SIF"

# Get neighborhood of gene (SIF)
curl -o tp53_network.sif "https://www.pathwaycommons.org/pc2/graph?source=TP53&kind=neighborhood&format=SIF"

# Get paths between two genes
curl -o paths.sif "https://www.pathwaycommons.org/pc2/graph?source=TP53&target=MDM2&kind=pathsbetween&format=SIF"

# Get Extended SIF with mediators
curl -o paths_ext.txt "https://www.pathwaycommons.org/pc2/graph?source=TP53&target=MDM2&kind=pathsbetween&format=EXTENDED_BINARY_SIF"

# Get JSON-LD format
curl -o pathway.json "https://www.pathwaycommons.org/pc2/get?uri=http://identifiers.org/reactome/R-HSA-109581&format=JSONLD"
```

### Method 4: Paxtools (Java Library)

```bash
# Download Paxtools
wget https://github.com/BioPAX/Paxtools/releases/download/v5.1.0/paxtools-5.1.0.jar

# Convert BioPAX to SIF
java -jar paxtools-5.1.0.jar toSIF input.owl output.sif \
  seqDb=uniprot,hgnc \
  chemDb=chebi

# Merge multiple BioPAX files
java -jar paxtools-5.1.0.jar merge output.owl input1.owl input2.owl

# Validate BioPAX
java -jar paxtools-5.1.0.jar validate input.owl
```

## File Inventory

### SIF Format Files

| File | Size | Description |
|------|------|-------------|
| PathwayCommons12.All.hgnc.sif.gz | ~25 MB | All interactions (HGNC symbols) |
| PathwayCommons12.All.uniprot.sif.gz | ~25 MB | All interactions (UniProt IDs) |
| PathwayCommons12.All.hgnc.txt.gz | ~100 MB | Extended SIF with annotations |
| PathwayCommons12.reactome.hgnc.sif.gz | ~5 MB | Reactome only |
| PathwayCommons12.kegg.hgnc.sif.gz | ~2 MB | KEGG only |

### BioPAX Files

| File | Size | Description |
|------|------|-------------|
| PathwayCommons12.All.BIOPAX.owl.gz | ~1.5 GB | Complete BioPAX database |
| PathwayCommons12.reactome.BIOPAX.owl.gz | ~300 MB | Reactome pathways |
| PathwayCommons12.kegg.BIOPAX.owl.gz | ~50 MB | KEGG pathways |
| PathwayCommons12.wp.BIOPAX.owl.gz | ~100 MB | WikiPathways |
| PathwayCommons12.biogrid.BIOPAX.owl.gz | ~400 MB | BioGRID interactions |

### Gene Set Files

| File | Size | Description |
|------|------|-------------|
| PathwayCommons12.All.hgnc.gmt.gz | ~5 MB | GSEA gene set format |

## Post-Download Processing

### Parse SIF Network

```bash
# Count interactions by type
gunzip -c PathwayCommons12.All.hgnc.sif.gz | \
  cut -f2 | sort | uniq -c | sort -rn

# Extract specific interaction type
gunzip -c PathwayCommons12.All.hgnc.sif.gz | \
  grep "CONTROLS-PHOSPHORYLATION-OF" > phosphorylation_network.sif

# Get unique genes
gunzip -c PathwayCommons12.All.hgnc.sif.gz | \
  cut -f1,3 | tr '\t' '\n' | sort -u > all_genes.txt

# Extract TP53 interactions
gunzip -c PathwayCommons12.All.hgnc.sif.gz | \
  grep -w "TP53" > tp53_interactions.sif
```

### Parse Extended SIF

```bash
# View header
gunzip -c PathwayCommons12.All.hgnc.txt.gz | head -1

# Extract interactions with PubMed IDs
gunzip -c PathwayCommons12.All.hgnc.txt.gz | \
  awk -F'\t' '$5 != ""' > interactions_with_pubmed.txt

# Filter by source database
gunzip -c PathwayCommons12.All.hgnc.txt.gz | \
  grep -i "reactome" > reactome_interactions.txt
```

### Process BioPAX with Python

```python
# Install: pip install pybiopax rdflib
from rdflib import Graph

# Load BioPAX
g = Graph()
g.parse("PathwayCommons12.All.BIOPAX.owl", format="xml")

# Query pathways
query = """
PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
SELECT ?pathway ?name WHERE {
  ?pathway a bp:Pathway .
  ?pathway bp:displayName ?name .
} LIMIT 100
"""
results = g.query(query)
for row in results:
    print(row.name)

# Extract proteins from pathway
query = """
PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
SELECT DISTINCT ?protein ?name WHERE {
  ?pathway a bp:Pathway .
  ?pathway bp:displayName "Apoptosis" .
  ?pathway bp:pathwayComponent+ ?component .
  ?component bp:participant ?protein .
  ?protein a bp:Protein .
  ?protein bp:displayName ?name .
}
"""
```

### Convert to NetworkX

```python
import networkx as nx

# Load SIF as graph
G = nx.DiGraph()
with open("PathwayCommons12.All.hgnc.sif") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 3:
            source, interaction, target = parts
            G.add_edge(source, target, type=interaction)

# Analyze network
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")

# Find shortest path
path = nx.shortest_path(G, "TP53", "AKT1")
print(f"Path: {' -> '.join(path)}")

# Get neighbors
neighbors = list(G.neighbors("TP53"))
print(f"TP53 targets: {neighbors[:10]}")
```

### Load into Cytoscape

```bash
# Install Cytoscape, then:
# 1. File -> Import -> Network from File
# 2. Select .sif file
# 3. Use default SIF format settings

# Command line (with cyREST)
curl -X POST "http://localhost:1234/v1/networks?source=file" \
  -H "Content-Type: application/json" \
  -d '{"sourceLocation": "/path/to/network.sif"}'
```

## Verification

```bash
# Check file integrity
gunzip -t PathwayCommons12.All.hgnc.sif.gz && echo "SIF OK"
gunzip -t PathwayCommons12.All.BIOPAX.owl.gz && echo "BioPAX OK"

# Count SIF records
gunzip -c PathwayCommons12.All.hgnc.sif.gz | wc -l
# Expected: ~2,400,000

# Verify SIF format (3 columns)
gunzip -c PathwayCommons12.All.hgnc.sif.gz | \
  awk -F'\t' 'NF != 3 {print "Error at line " NR; exit 1}'

# Check BioPAX validity (basic)
gunzip -c PathwayCommons12.All.BIOPAX.owl.gz | \
  head -100 | grep -q "biopax-level3" && echo "BioPAX Level 3 confirmed"

# Count pathways in BioPAX
gunzip -c PathwayCommons12.All.BIOPAX.owl.gz | \
  grep -c "<bp:Pathway"
# Expected: ~5,700

# Verify data sources present
gunzip -c PathwayCommons12.All.hgnc.txt.gz | \
  cut -f4 | sort -u | head -20
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major version | Annual |
| Data refresh | Quarterly |
| Source updates | Continuous |

## Common Issues

- **Large BioPAX files**: Use streaming XML parser or Paxtools
- **Memory limits**: Process SIF format instead of BioPAX for network analysis
- **URI resolution**: Some URIs may be outdated; check Identifiers.org
- **Gene symbol changes**: Use UniProt version for stable identifiers
- **Missing edges**: Some sources have restrictive licenses (KEGG partial)

## API Query Examples

```bash
# Search for TP53 pathways
curl "https://www.pathwaycommons.org/pc2/search?q=TP53&type=Pathway&organism=9606"

# Get top pathways for a gene set
curl "https://www.pathwaycommons.org/pc2/top_pathways?q=TP53+BRCA1+ATM"

# Traverse relationships
curl "https://www.pathwaycommons.org/pc2/traverse?uri=http://identifiers.org/uniprot/P04637&path=ProteinReference/entityReferenceOf/participantOf"

# Download specific data source
curl -o reactome.owl "https://www.pathwaycommons.org/pc2/get?source=reactome"
```

## Data Formats Reference

### SIF (Simple Interaction Format)

```
SOURCE_GENE	INTERACTION_TYPE	TARGET_GENE
TP53	INTERACTS_WITH	MDM2
TP53	CONTROLS-STATE-CHANGE-OF	BCL2
```

### Extended SIF

```
PARTICIPANT_A	INTERACTION_TYPE	PARTICIPANT_B	INTERACTION_DATA_SOURCE	INTERACTION_PUBMED_ID	PATHWAY_NAMES	MEDIATOR_IDS
TP53	controls-phosphorylation-of	MDM2	reactome	11900253	TP53 Regulation	CHEK2
```

### GMT (Gene Set)

```
Apoptosis%Reactome%R-HSA-109581	http://identifiers.org/reactome/R-HSA-109581	TP53	BCL2	CASP3	...
```

## License

Free for all uses. Attribution to Pathway Commons encouraged. Original data sources retain their licenses.
