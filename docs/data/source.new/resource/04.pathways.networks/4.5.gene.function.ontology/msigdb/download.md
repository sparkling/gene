---
id: download-msigdb
title: "MSigDB Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# MSigDB - Download Documentation

## Overview

MSigDB provides gene sets via web download (requires registration) and a REST API. Gene sets are available in multiple formats and identifier types.

## Registration

Free registration required at: https://www.gsea-msigdb.org/gsea/register.jsp

## Web Downloads

### Download Portal

```
https://www.gsea-msigdb.org/gsea/downloads.jsp
```

### Available Files

| Collection | Symbols | Entrez |
|------------|---------|--------|
| Hallmark (H) | h.all.v2024.1.Hs.symbols.gmt | h.all.v2024.1.Hs.entrez.gmt |
| Curated (C2) | c2.all.v2024.1.Hs.symbols.gmt | c2.all.v2024.1.Hs.entrez.gmt |
| GO (C5) | c5.all.v2024.1.Hs.symbols.gmt | c5.all.v2024.1.Hs.entrez.gmt |
| Immunologic (C7) | c7.all.v2024.1.Hs.symbols.gmt | c7.all.v2024.1.Hs.entrez.gmt |
| Cell Type (C8) | c8.all.v2024.1.Hs.symbols.gmt | c8.all.v2024.1.Hs.entrez.gmt |
| All Collections | msigdb.v2024.1.Hs.symbols.gmt | msigdb.v2024.1.Hs.entrez.gmt |

### Subcollection Downloads

```
# KEGG pathways only
c2.cp.kegg.v2024.1.Hs.symbols.gmt

# Reactome pathways only
c2.cp.reactome.v2024.1.Hs.symbols.gmt

# GO Biological Process only
c5.go.bp.v2024.1.Hs.symbols.gmt

# TF targets
c3.tft.v2024.1.Hs.symbols.gmt
```

## REST API

### Base URL

```
https://www.gsea-msigdb.org/gsea/msigdb/api
```

### Search Gene Sets

```bash
# Search by keyword
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets?keyword=apoptosis"

# Search in specific collection
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets?keyword=apoptosis&collection=H"

# Search by gene
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets?gene=TP53"

# Pagination
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets?keyword=cancer&offset=0&limit=100"
```

### Get Gene Set Details

```bash
# By name
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets/HALLMARK_APOPTOSIS"

# By systematic name
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets/M5930"
```

### Get Genes in Gene Set

```bash
# Get gene list
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets/HALLMARK_APOPTOSIS/genes"

# With specific identifier type
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets/HALLMARK_APOPTOSIS/genes?idType=entrez"
```

### List Collections

```bash
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/collections"
```

### Download Gene Set (GMT)

```bash
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets/HALLMARK_APOPTOSIS/gmt"
```

## Python Examples

### Load GMT File

```python
def load_gmt(filepath):
    """Load gene sets from GMT file."""
    gene_sets = {}

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            description = parts[1]
            genes = parts[2:]

            gene_sets[name] = {
                'description': description,
                'genes': set(genes)
            }

    return gene_sets

# Example
hallmark = load_gmt('h.all.v2024.1.Hs.symbols.gmt')
print(f"Loaded {len(hallmark)} gene sets")
print(f"HALLMARK_APOPTOSIS has {len(hallmark['HALLMARK_APOPTOSIS']['genes'])} genes")
```

### API Access

```python
import requests

class MSigDBClient:
    def __init__(self):
        self.base_url = "https://www.gsea-msigdb.org/gsea/msigdb/api/v1"

    def search_gene_sets(self, keyword, collection=None, limit=100):
        """Search gene sets by keyword."""
        params = {"keyword": keyword, "limit": limit}
        if collection:
            params["collection"] = collection

        response = requests.get(f"{self.base_url}/gene_sets", params=params)
        return response.json()

    def get_gene_set(self, name):
        """Get gene set details."""
        response = requests.get(f"{self.base_url}/gene_sets/{name}")
        return response.json()

    def get_genes(self, name, id_type="symbol"):
        """Get genes in a gene set."""
        response = requests.get(
            f"{self.base_url}/gene_sets/{name}/genes",
            params={"idType": id_type}
        )
        return response.json()

# Example usage
client = MSigDBClient()
results = client.search_gene_sets("hypoxia", collection="H")
genes = client.get_genes("HALLMARK_HYPOXIA")
```

### Enrichment Analysis with GSEApy

```python
import gseapy as gp

# Preranked GSEA
rnk = pd.read_csv('gene_ranks.txt', sep='\t', header=None, names=['gene', 'rank'])

result = gp.prerank(
    rnk=rnk,
    gene_sets='h.all.v2024.1.Hs.symbols.gmt',
    outdir='gsea_results',
    min_size=15,
    max_size=500,
    permutation_num=1000
)

# View results
result.res2d.sort_values('NES', ascending=False).head(10)
```

### Over-Representation Analysis

```python
import gseapy as gp

# Gene list
my_genes = ['TP53', 'BRCA1', 'EGFR', 'KRAS', 'MYC', ...]

# ORA with MSigDB
result = gp.enrichr(
    gene_list=my_genes,
    gene_sets='MSigDB_Hallmark_2020',
    organism='Human',
    outdir='enrichr_results'
)
```

## GSEA Desktop Software

### Download

```
https://www.gsea-msigdb.org/gsea/downloads.jsp
```

### Requirements

- Java 11+
- 4GB RAM recommended

### Running GSEA

```bash
# Preranked analysis
java -Xmx4g -jar gsea-4.3.2.jar GSEA \
  -res expression.txt \
  -cls phenotype.cls \
  -gmx h.all.v2024.1.Hs.symbols.gmt \
  -out results/ \
  -rpt_label my_analysis
```

## File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| GMT | .gmt | Gene Matrix Transposed |
| GMX | .gmx | Gene Matrix (column-oriented) |
| GRP | .grp | Single gene list |
| XML | .xml | Full metadata |

## Species Support

| Species | Code | Available |
|---------|------|-----------|
| Human | Hs | Full database |
| Mouse | Mm | Orthologs |
| Rat | Rn | Orthologs |

## Versioning

Format: `{collection}.v{YYYY.N}.{Species}.{idtype}.gmt`

Example: `h.all.v2024.1.Hs.symbols.gmt`

- 2024.1 = First release of 2024
- Hs = Homo sapiens
- symbols = HGNC gene symbols

## See Also

- [Schema Documentation](./schema.md)
- [GSEA User Guide](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html)
