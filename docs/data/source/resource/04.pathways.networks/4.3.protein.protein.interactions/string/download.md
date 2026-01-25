---
id: download-string
title: "STRING Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# STRING Download Instructions

## Quick Start

```bash
# Download human protein links
wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- 10-500GB disk space depending on species coverage

## No Registration Required

STRING data is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Species-Specific Downloads

```bash
# Human protein links (combined score)
wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz

# Human detailed interactions (with evidence channels)
wget https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz

# Human physical interactions only
wget https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz

# Human protein aliases (ID mapping)
wget https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz

# Human protein info
wget https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
```

### Method 2: All Species Bulk Download

```bash
# All protein links (LARGE: ~50GB compressed)
wget https://stringdb-downloads.org/download/protein.links.v12.0.txt.gz

# All protein info
wget https://stringdb-downloads.org/download/protein.info.v12.0.txt.gz

# Species list
wget https://stringdb-downloads.org/download/species.v12.0.txt
```

### Method 3: Evidence-Specific Downloads

```bash
# Experimental interactions
wget https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz

# Text mining interactions
wget https://stringdb-downloads.org/download/protein.textmining.v12.0/9606.protein.textmining.v12.0.txt.gz

# Homology-based transfers
wget https://stringdb-downloads.org/download/protein.homology.v12.0/9606.protein.homology.v12.0.txt.gz
```

### Method 4: Network Images

```bash
# Download network image for a protein
curl "https://string-db.org/api/image/network?identifiers=TP53&species=9606" -o tp53_network.png

# Higher resolution
curl "https://string-db.org/api/highres_image/network?identifiers=TP53&species=9606" -o tp53_network_hires.png
```

### Method 5: STRING API

```bash
# Get interactions for multiple proteins
curl "https://string-db.org/api/tsv/network?identifiers=TP53%0ABRCA1%0ABRCA2&species=9606" -o tp53_network.tsv

# Get interaction partners
curl "https://string-db.org/api/tsv/interaction_partners?identifiers=TP53&species=9606&limit=100" -o tp53_partners.tsv

# Get enrichment analysis
curl "https://string-db.org/api/tsv/enrichment?identifiers=TP53%0ABRCA1%0ABRCA2%0AATM&species=9606" -o enrichment.tsv

# Map identifiers
curl "https://string-db.org/api/tsv/get_string_ids?identifiers=p53%0Abrca1&species=9606" -o mapped_ids.tsv
```

### Method 6: Clusters/COG

```bash
# Clusters of Orthologous Groups
wget https://stringdb-downloads.org/download/COG.links.v12.0.txt.gz
wget https://stringdb-downloads.org/download/COG.mappings.v12.0.txt.gz
```

## File Inventory

### Per-Species Files

| File Pattern | Size (Human) | Description |
|--------------|--------------|-------------|
| *.protein.links.*.txt.gz | ~100 MB | Combined scores |
| *.protein.links.full.*.txt.gz | ~300 MB | All evidence channels |
| *.protein.physical.links.*.txt.gz | ~30 MB | Physical only |
| *.protein.aliases.*.txt.gz | ~50 MB | ID mapping |
| *.protein.info.*.txt.gz | ~10 MB | Protein metadata |

### Global Files

| File | Size | Description |
|------|------|-------------|
| protein.links.*.txt.gz | ~50 GB | All species |
| protein.info.*.txt.gz | ~5 GB | All proteins |
| species.*.txt | ~100 KB | Species list |

## Post-Download Processing

```bash
# Decompress
gunzip 9606.protein.links.v12.0.txt.gz

# Filter by score threshold (>700 = high confidence)
awk '$3 >= 700' 9606.protein.links.v12.0.txt > high_confidence_links.txt

# Convert STRING IDs to gene symbols
python3 << 'EOF'
import pandas as pd

# Load mapping
aliases = pd.read_csv('9606.protein.aliases.v12.0.txt.gz', sep='\t', compression='gzip')
gene_map = aliases[aliases['source'] == 'Ensembl_HGNC_symbol'][['#string_protein_id', 'alias']]
gene_map = gene_map.drop_duplicates(subset='#string_protein_id')
gene_map = dict(zip(gene_map['#string_protein_id'], gene_map['alias']))

# Load interactions
links = pd.read_csv('9606.protein.links.v12.0.txt', sep=' ')
links['protein1_symbol'] = links['protein1'].map(gene_map)
links['protein2_symbol'] = links['protein2'].map(gene_map)
links.to_csv('string_with_symbols.tsv', sep='\t', index=False)
EOF

# Create adjacency matrix
python3 << 'EOF'
import pandas as pd
import numpy as np

links = pd.read_csv('high_confidence_links.txt', sep=' ')
proteins = sorted(set(links['protein1'].tolist() + links['protein2'].tolist()))
n = len(proteins)
protein_idx = {p: i for i, p in enumerate(proteins)}

matrix = np.zeros((n, n), dtype=np.float16)
for _, row in links.iterrows():
    i, j = protein_idx[row['protein1']], protein_idx[row['protein2']]
    score = row['combined_score'] / 1000
    matrix[i, j] = score
    matrix[j, i] = score

# Save as sparse format for large networks
from scipy import sparse
sparse.save_npz('string_matrix.npz', sparse.csr_matrix(matrix))
EOF

# Extract subnetwork for gene list
grep -E "(ENSP00000269305|ENSP00000275493|ENSP00000418960)" 9606.protein.links.v12.0.txt
```

## Verification

```bash
# Check file structure
zcat 9606.protein.links.v12.0.txt.gz | head -5

# Count interactions
zcat 9606.protein.links.v12.0.txt.gz | wc -l

# Count unique proteins
zcat 9606.protein.links.v12.0.txt.gz | awk '{print $1"\n"$2}' | sort -u | wc -l

# Distribution of scores
zcat 9606.protein.links.v12.0.txt.gz | awk '{print int($3/100)*100}' | sort -n | uniq -c
```

---

## Dataset Versions

### Current Release: STRING v12.0

| Property | Value |
|----------|-------|
| Version | 12.0 |
| Release Date | 2023-08-21 |
| Total Size | ~50 GB (all species) |
| Organisms | 14,094 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| 9606.protein.links.v12.0.txt.gz | ~100 MB | ~12M | Human combined scores |
| 9606.protein.links.full.v12.0.txt.gz | ~300 MB | ~12M | Human all channels |
| 9606.protein.aliases.v12.0.txt.gz | ~50 MB | ~20K proteins | ID mapping |
| protein.links.v12.0.txt.gz | ~50 GB | All | All species |

### Previous Versions

| Version | Release | Organisms | Status |
|---------|---------|-----------|--------|
| 11.5 | 2021-08-12 | 14,094 | Archived |
| 11.0 | 2020-12-08 | 5,090 | Archived |
| 10.5 | 2017-05-14 | 2,031 | Archived |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| Base URL | `https://string-db.org/api` |
| Authentication | None required |
| Rate Limit | 1 request/second recommended |
| Response Format | TSV, JSON, PNG, SVG |

### API Endpoints

| Operation | Endpoint | Example |
|-----------|----------|---------|
| Network | `/tsv/network` | `?identifiers=TP53%0ABRCA1&species=9606` |
| Partners | `/tsv/interaction_partners` | `?identifiers=TP53&species=9606&limit=100` |
| Enrichment | `/tsv/enrichment` | `?identifiers=TP53%0ABRCA1%0ABRCA2&species=9606` |
| Map IDs | `/tsv/get_string_ids` | `?identifiers=p53&species=9606` |
| Image | `/image/network` | `?identifiers=TP53&species=9606` |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major versions | Every 1-2 years |
| Species updates | With major releases |
| Bug fixes | As needed |

## Common Issues

- **STRING IDs**: Use protein.aliases for mapping to gene symbols
- **Score interpretation**: 0-1000 scale; >700 is high confidence
- **Bidirectional edges**: Each interaction appears once (A-B, not B-A)
- **Self-loops**: Some proteins have self-interactions
- **Memory issues**: Full network is large; filter by species/score

## Score Channels

| Channel | Description |
|---------|-------------|
| neighborhood | Gene neighborhood |
| fusion | Gene fusion events |
| cooccurence | Phylogenetic co-occurrence |
| coexpression | Co-expression |
| experimental | Experimental data |
| database | Curated databases |
| textmining | Text mining |
| combined_score | Combined (main score) |

## Taxonomy IDs (Common Species)

| Species | NCBI Tax ID |
|---------|-------------|
| Human | 9606 |
| Mouse | 10090 |
| Rat | 10116 |
| Zebrafish | 7955 |
| D. melanogaster | 7227 |
| C. elegans | 6239 |
| S. cerevisiae | 4932 |
| E. coli K-12 | 511145 |

## Related Resources

- [BioGRID](../biogrid/) - Curated interactions
- [IntAct](../intact/) - Molecular interactions
- [Reactome](../../4.1.metabolic.pathways/reactome/) - Pathway context
