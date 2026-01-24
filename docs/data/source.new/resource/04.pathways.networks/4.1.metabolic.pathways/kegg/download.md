---
id: download-kegg
title: "KEGG Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# KEGG Download Instructions

## Quick Start

```bash
# Download pathway list via API
curl "https://rest.kegg.jp/list/pathway/hsa" -o hsa_pathways.txt

# Download single pathway (KGML format)
curl "https://rest.kegg.jp/get/hsa05200/kgml" -o hsa05200.kgml
```

## Prerequisites

- **curl** or **wget** for API access
- **Python** with biopython for parsing
- Understanding of KEGG licensing restrictions

## Important: Licensing Considerations

KEGG has specific licensing requirements:

| Use Type | Access |
|----------|--------|
| Academic (non-bulk) | Free via API/website |
| Academic (bulk/FTP) | Requires academic subscription |
| Commercial | Requires commercial license |

Bulk FTP downloads require a KEGG subscription. The methods below focus on API access.

## Download Methods

### Method 1: REST API (Recommended for Small-Scale)

```bash
# List all human pathways
curl "https://rest.kegg.jp/list/pathway/hsa" -o human_pathways.txt

# Get pathway info
curl "https://rest.kegg.jp/get/hsa05200" -o pathway_info.txt

# Get pathway image
curl "https://rest.kegg.jp/get/hsa05200/image" -o hsa05200.png

# Get KGML (pathway markup)
curl "https://rest.kegg.jp/get/hsa05200/kgml" -o hsa05200.kgml

# Get pathway genes
curl "https://rest.kegg.jp/link/hsa/pathway" -o pathway_genes.txt
```

### Method 2: Batch Downloads via API

```bash
# Download all human pathway KGMLs
curl "https://rest.kegg.jp/list/pathway/hsa" | cut -f1 | while read pathway; do
  echo "Downloading ${pathway}..."
  curl "https://rest.kegg.jp/get/${pathway}/kgml" -o "${pathway}.kgml"
  sleep 1  # Be respectful to the server
done

# Download compound information
curl "https://rest.kegg.jp/list/compound" -o compounds.txt
curl "https://rest.kegg.jp/list/compound" | head -100 | cut -f1 | while read cpd; do
  curl "https://rest.kegg.jp/get/${cpd}" -o "compounds/${cpd}.txt"
  sleep 0.5
done
```

### Method 3: Gene/Pathway Mappings

```bash
# Gene to pathway mapping
curl "https://rest.kegg.jp/link/pathway/hsa" -o gene_pathway_map.txt

# Pathway to gene mapping
curl "https://rest.kegg.jp/link/hsa/pathway" -o pathway_gene_map.txt

# Gene to KEGG ortholog
curl "https://rest.kegg.jp/link/ko/hsa" -o gene_ko_map.txt

# Disease to gene mapping
curl "https://rest.kegg.jp/link/hsa/disease" -o disease_gene_map.txt
```

### Method 4: Compound and Drug Data

```bash
# List all compounds
curl "https://rest.kegg.jp/list/compound" -o compounds.txt

# List all drugs
curl "https://rest.kegg.jp/list/drug" -o drugs.txt

# Get compound structure (MOL format)
curl "https://rest.kegg.jp/get/C00002/mol" -o C00002.mol

# Get drug info
curl "https://rest.kegg.jp/get/D00001" -o aspirin.txt
```

### Method 5: KEGG Orthology (KO)

```bash
# List all KO entries
curl "https://rest.kegg.jp/list/ko" -o ko_list.txt

# Get KO entry details
curl "https://rest.kegg.jp/get/K00001" -o K00001.txt

# KO to pathway mapping
curl "https://rest.kegg.jp/link/pathway/ko" -o ko_pathway_map.txt
```

### Method 6: FTP Access (Subscription Required)

```bash
# For subscribers only
ftp ftp.kegg.jp
# Login with subscription credentials
# cd /kegg/pathway/organisms/hsa
# mget *.kgml
```

## File Inventory

### Via API (Free)

| Data Type | Format | Size Estimate |
|-----------|--------|---------------|
| Pathway list | TSV | ~50 KB |
| KGML files | XML | ~50 KB each |
| Gene mapping | TSV | ~5 MB |
| Compound list | TSV | ~1 MB |
| KO list | TSV | ~2 MB |

### Via FTP (Subscription)

| Data Type | Size | Description |
|-----------|------|-------------|
| All pathways | ~500 MB | KGML format |
| Genes | ~100 MB | All organisms |
| Compounds | ~50 MB | Structures, info |
| Reactions | ~100 MB | Enzyme reactions |

## Post-Download Processing

```bash
# Parse pathway list
awk -F'\t' '{print $1"\t"$2}' human_pathways.txt

# Extract genes from KGML
python3 << 'EOF'
import xml.etree.ElementTree as ET

tree = ET.parse('hsa05200.kgml')
root = tree.getroot()

for entry in root.findall('.//entry[@type="gene"]'):
    gene_ids = entry.get('name').replace('hsa:', '').split()
    for gene_id in gene_ids:
        print(gene_id)
EOF

# Create pathway-gene matrix
python3 << 'EOF'
import pandas as pd

# Read gene-pathway mapping
df = pd.read_csv('gene_pathway_map.txt', sep='\t', header=None, names=['gene', 'pathway'])
df['gene'] = df['gene'].str.replace('hsa:', '')
df['pathway'] = df['pathway'].str.replace('path:', '')

# Pivot to matrix
matrix = pd.crosstab(df['gene'], df['pathway'])
matrix.to_csv('gene_pathway_matrix.csv')
EOF

# Convert KGML to networkx graph
python3 << 'EOF'
import networkx as nx
import xml.etree.ElementTree as ET

tree = ET.parse('hsa05200.kgml')
root = tree.getroot()

G = nx.DiGraph()
for relation in root.findall('.//relation'):
    entry1 = relation.get('entry1')
    entry2 = relation.get('entry2')
    rel_type = relation.find('subtype').get('name') if relation.find('subtype') is not None else 'unknown'
    G.add_edge(entry1, entry2, type=rel_type)

print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
EOF
```

## Verification

```bash
# Check pathway count
wc -l human_pathways.txt

# Verify KGML structure
head -20 hsa05200.kgml

# Check gene mapping
cut -f1 gene_pathway_map.txt | sort -u | wc -l
```

## Update Schedule

| Data Type | Frequency |
|-----------|-----------|
| Pathways | Monthly |
| Compounds | Monthly |
| Genes | Monthly |
| Diseases | Quarterly |

## Common Issues

- **Rate limiting**: Add delays between API calls (1 request/second recommended)
- **Bulk restrictions**: FTP requires subscription; API has implicit limits
- **ID formats**: KEGG uses specific prefixes (hsa:, path:, cpd:, etc.)
- **Data currency**: Some data may lag behind primary sources
- **Commercial use**: Requires explicit license

## KEGG Database Codes

| Code | Database |
|------|----------|
| path | Pathway |
| hsa | Human genes |
| ko | KEGG Orthology |
| cpd | Compound |
| dr | Drug |
| ds | Disease |
| br | BRITE hierarchy |
| md | Module |

## API Operations Reference

| Operation | URL Pattern |
|-----------|-------------|
| list | rest.kegg.jp/list/{db} |
| get | rest.kegg.jp/get/{id} |
| find | rest.kegg.jp/find/{db}/{query} |
| link | rest.kegg.jp/link/{db1}/{db2} |
| conv | rest.kegg.jp/conv/{db1}/{db2} |

## Alternative Free Resources

For unrestricted bulk downloads, consider:
- **Reactome** - Open pathway data
- **WikiPathways** - Community curated
- **Pathway Commons** - Aggregated pathways

## Related Resources

- [Reactome](../reactome/) - Alternative pathway database
- [GO](../../03.diseases.phenotypes/3.1.disease.ontologies/) - Gene Ontology
- [ChEMBL](../../02.compounds.molecules/2.2.pharmaceuticals/chembl/) - Bioactivity data
