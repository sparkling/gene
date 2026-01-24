---
id: download-reactome
title: "Reactome Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Reactome Download Instructions

## Quick Start

```bash
# Download pathway hierarchy
wget https://reactome.org/download/current/ReactomePathways.txt

# Download all pathways (BioPAX format)
wget https://reactome.org/download/current/biopax.zip
```

## Prerequisites

- **wget** or **curl** for downloads
- **Java** for BioPAX tools
- **Python** for API access
- Approximately 5-20GB disk space

## No Registration Required

Reactome data is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Core Data Files

```bash
# Pathway hierarchy
wget https://reactome.org/download/current/ReactomePathways.txt

# Pathway relationships (parent-child)
wget https://reactome.org/download/current/ReactomePathwaysRelation.txt

# Gene-to-pathway mapping (human)
wget https://reactome.org/download/current/NCBI2Reactome.txt

# UniProt to pathway mapping
wget https://reactome.org/download/current/UniProt2Reactome.txt

# All reactions
wget https://reactome.org/download/current/ReactionLikeEvent.txt
```

### Method 2: BioPAX Format

```bash
# Complete pathway data (BioPAX Level 3)
wget https://reactome.org/download/current/biopax.zip
unzip biopax.zip

# Species-specific (human)
wget https://reactome.org/download/current/biopax/Homo_sapiens.owl

# Individual pathway
# Use API to get specific pathway BioPAX
curl "https://reactome.org/ContentService/exporter/event/R-HSA-109582.owl" -o p53_pathway.owl
```

### Method 3: SBML/SBGN Formats

```bash
# SBML files
wget https://reactome.org/download/current/homo_sapiens.sbml.tar.bz2
tar -xjf homo_sapiens.sbml.tar.bz2

# SBGN files
wget https://reactome.org/download/current/homo_sapiens.sbgn.tar.bz2
```

### Method 4: GMT Files (Gene Set Format)

```bash
# GMT format for GSEA
wget https://reactome.org/download/current/ReactomePathways.gmt.zip
unzip ReactomePathways.gmt.zip

# Human GMT
# Included in the zip as "ReactomePathways.gmt"
```

### Method 5: Graph Database Dump (Neo4j)

```bash
# Neo4j database dump
wget https://reactome.org/download/current/reactome.graphdb.dump

# Restore to Neo4j
neo4j-admin load --from=reactome.graphdb.dump --database=reactome
```

### Method 6: REST API Downloads

```bash
# Get pathway details
curl "https://reactome.org/ContentService/data/pathway/R-HSA-109582/containedEvents"

# Get pathway hierarchy for human
curl "https://reactome.org/ContentService/data/eventsHierarchy/9606" -o human_hierarchy.json

# Export pathway as image
curl "https://reactome.org/ContentService/exporter/diagram/R-HSA-109582.png" -o p53_diagram.png

# Export as PDF
curl "https://reactome.org/ContentService/exporter/document/event/R-HSA-109582.pdf" -o p53_pathway.pdf
```

### Method 7: Analysis Service Data

```bash
# Download interactors
wget https://reactome.org/download/current/interactors.db

# Gene expression data format specs
wget https://reactome.org/download/current/analysis_service_reference_data.zip
```

## File Inventory

### Core Files

| File | Size | Description |
|------|------|-------------|
| ReactomePathways.txt | ~5 MB | Pathway list |
| ReactomePathwaysRelation.txt | ~1 MB | Hierarchy |
| NCBI2Reactome.txt | ~20 MB | Gene mapping |
| UniProt2Reactome.txt | ~50 MB | Protein mapping |

### Pathway Representations

| File | Size | Description |
|------|------|-------------|
| biopax.zip | ~2 GB | BioPAX Level 3 |
| homo_sapiens.sbml.tar.bz2 | ~500 MB | SBML format |
| homo_sapiens.sbgn.tar.bz2 | ~300 MB | SBGN format |

### Gene Sets

| File | Size | Description |
|------|------|-------------|
| ReactomePathways.gmt.zip | ~10 MB | GSEA format |

### Database

| File | Size | Description |
|------|------|-------------|
| reactome.graphdb.dump | ~10 GB | Neo4j dump |
| interactors.db | ~1 GB | Interactors database |

## Post-Download Processing

```bash
# Parse pathway hierarchy
awk -F'\t' '{print $1"\t"$2"\t"$3}' ReactomePathways.txt | head

# Extract human pathways
grep "Homo sapiens" ReactomePathways.txt > human_pathways.txt

# Create gene-pathway mapping
awk -F'\t' '$6=="Homo sapiens" {print $1"\t"$2"\t"$4}' NCBI2Reactome.txt > gene_pathway_map.tsv

# Convert GMT to gene sets
python3 << 'EOF'
with open('ReactomePathways.gmt') as f:
    for line in f:
        parts = line.strip().split('\t')
        pathway_id = parts[0]
        pathway_name = parts[1]
        genes = parts[2:]
        print(f"{pathway_id}\t{pathway_name}\t{len(genes)} genes")
EOF

# Load BioPAX with paxtools
java -jar paxtools.jar toSIF Homo_sapiens.owl output.sif
```

## Verification

```bash
# Check file integrity
head -5 ReactomePathways.txt

# Count pathways by species
cut -f3 ReactomePathways.txt | sort | uniq -c | sort -rn | head

# Count gene-pathway associations (human)
grep "Homo sapiens" NCBI2Reactome.txt | wc -l
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major releases | Quarterly |
| Data updates | Monthly |
| Bug fixes | As needed |

## Common Issues

- **BioPAX parsing**: Use PaxTools or pybiopax library
- **Large files**: Use streaming for Neo4j import
- **Identifier mapping**: Map via UniProt or NCBI Gene ID
- **Pathway versioning**: Stable IDs have format R-HSA-NNNNNN
- **Species coverage**: Most complete for human; varies for others

## Pathway Hierarchy

| Level | Example |
|-------|---------|
| Top Level | Signal Transduction |
| Second Level | Signaling by Receptor Tyrosine Kinases |
| Third Level | Signaling by EGFR |
| Reaction | EGFR binds EGF |

## API Usage Examples

```python
import requests

# Search pathways
response = requests.get(
    "https://reactome.org/ContentService/search/query",
    params={"query": "apoptosis", "species": "Homo sapiens", "types": "Pathway"}
)
pathways = response.json()

# Get pathway participants
response = requests.get(
    "https://reactome.org/ContentService/data/participants/R-HSA-109582"
)
participants = response.json()

# Pathway analysis
genes = ["TP53", "BRCA1", "BRCA2", "ATM"]
response = requests.post(
    "https://reactome.org/AnalysisService/identifiers/projection",
    headers={"Content-Type": "text/plain"},
    data="\n".join(genes)
)
analysis = response.json()
```

## Related Resources

- [KEGG](../kegg/) - Alternative pathway database
- [Gene Ontology](../../03.diseases.phenotypes/3.1.disease.ontologies/) - Functional annotation
- [STRING](../../4.3.protein.protein.interactions/string/) - Protein interactions
