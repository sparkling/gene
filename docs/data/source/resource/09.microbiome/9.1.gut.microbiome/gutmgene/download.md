---
id: download-gutmgene
title: "gutMGene Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# gutMGene Download Instructions

## Quick Start

```bash
# Download via web interface
# Visit http://bio-annotation.cn/gutmgene and use Download section
```

## Prerequisites

- Web browser for download interface
- **Python** or **R** for data processing
- 100MB-1GB storage for full dataset

## No Registration Required

gutMGene data is freely available for academic use.

## Download Methods

### Method 1: Web Interface Download

1. Visit http://bio-annotation.cn/gutmgene
2. Navigate to "Download" section
3. Select data tables:
   - Microbe-gene associations
   - Microbe information
   - Gene information
4. Download TSV/Excel files

### Method 2: Query-Based Export

1. Visit http://bio-annotation.cn/gutmgene/search
2. Search by:
   - Microbe name/taxon
   - Gene symbol
   - Disease/condition
   - Tissue
3. Export search results to TSV

### Method 3: Programmatic Access (if available)

```python
import pandas as pd

# Note: Check website for current download URLs
# These are example patterns - verify actual URLs

# Download main association table
url = "http://bio-annotation.cn/gutmgene/data/associations.tsv"
associations = pd.read_csv(url, sep="\t")

# Download microbe table
microbes = pd.read_csv("http://bio-annotation.cn/gutmgene/data/microbes.tsv", sep="\t")

# Download gene table
genes = pd.read_csv("http://bio-annotation.cn/gutmgene/data/genes.tsv", sep="\t")
```

### Method 4: Manual Curation from Interface

```python
import requests
from bs4 import BeautifulSoup
import time

# Scrape search results (use responsibly)
def search_gutmgene(microbe=None, gene=None):
    """Search gutMGene database."""
    base_url = "http://bio-annotation.cn/gutmgene/search"
    params = {}
    if microbe:
        params['microbe'] = microbe
    if gene:
        params['gene'] = gene

    response = requests.get(base_url, params=params)
    # Parse results (implementation depends on page structure)
    return response.text

# Example: Search for Faecalibacterium associations
results = search_gutmgene(microbe="Faecalibacterium")
```

## File Inventory

### Expected Download Files

| File | Size (est.) | Description |
|------|-------------|-------------|
| associations.tsv | ~10 MB | All microbe-gene links |
| microbes.tsv | ~500 KB | Microbe taxonomy info |
| genes.tsv | ~2 MB | Host gene annotations |
| evidence.tsv | ~5 MB | Literature evidence |

### Data Fields (Associations Table)

| Column | Description |
|--------|-------------|
| microbe_taxid | NCBI Taxonomy ID |
| microbe_name | Species name |
| gene_symbol | Host gene symbol |
| gene_id | Entrez Gene ID |
| direction | up/down/altered |
| fold_change | Expression change |
| tissue | Affected tissue |
| pmid | PubMed reference |

## Post-Download Processing

```bash
# Convert to standard formats
python3 << 'EOF'
import pandas as pd

# Load associations
df = pd.read_csv("associations.tsv", sep="\t")

# Filter by direction
upregulated = df[df["direction"] == "up"]
downregulated = df[df["direction"] == "down"]

print(f"Total associations: {len(df)}")
print(f"Upregulated: {len(upregulated)}")
print(f"Downregulated: {len(downregulated)}")

# Group by microbe
microbe_counts = df.groupby("microbe_name").size().sort_values(ascending=False)
print("\nTop 10 microbes by association count:")
print(microbe_counts.head(10))

# Group by gene
gene_counts = df.groupby("gene_symbol").size().sort_values(ascending=False)
print("\nTop 10 genes by association count:")
print(gene_counts.head(10))

# Save filtered subsets
upregulated.to_csv("upregulated_associations.tsv", sep="\t", index=False)
downregulated.to_csv("downregulated_associations.tsv", sep="\t", index=False)
EOF
```

### Create Gene-Centric View

```python
import pandas as pd

# Load data
df = pd.read_csv("associations.tsv", sep="\t")

# Create gene-centric summary
gene_summary = df.groupby("gene_symbol").agg({
    "microbe_name": lambda x: ", ".join(sorted(set(x))),
    "direction": lambda x: list(x),
    "pmid": "count"
}).rename(columns={"pmid": "evidence_count"})

gene_summary.to_csv("gene_centric_summary.tsv", sep="\t")
```

### Create Microbe-Centric View

```python
import pandas as pd

# Create microbe-centric summary
microbe_summary = df.groupby("microbe_name").agg({
    "gene_symbol": lambda x: ", ".join(sorted(set(x))),
    "tissue": lambda x: ", ".join(sorted(set(x.dropna()))),
    "pmid": "nunique"
}).rename(columns={"pmid": "publication_count"})

microbe_summary.to_csv("microbe_centric_summary.tsv", sep="\t")
```

## Verification

```bash
# Check file integrity
wc -l associations.tsv

# Check column headers
head -1 associations.tsv

# Count unique entities
python3 << 'EOF'
import pandas as pd
df = pd.read_csv("associations.tsv", sep="\t")
print(f"Unique microbes: {df['microbe_taxid'].nunique()}")
print(f"Unique genes: {df['gene_id'].nunique()}")
print(f"Unique tissues: {df['tissue'].nunique()}")
print(f"Unique PMIDs: {df['pmid'].nunique()}")
EOF
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major updates | Annually |
| Literature curation | Continuous |

## Common Issues

- **No direct download link**: Use web interface
- **Rate limiting**: Add delays if scraping
- **Encoding issues**: Ensure UTF-8 handling
- **Missing data**: Some fields may be null
- **Commercial use**: Contact maintainers

## Integration with Other Databases

```python
# Link to GMrepo samples
import pandas as pd

gutmgene = pd.read_csv("associations.tsv", sep="\t")
# Get unique microbe taxids for GMrepo query
microbe_taxids = gutmgene["microbe_taxid"].unique()

# Link to KEGG pathways
# Use gene symbols to query KEGG for pathway enrichment

# Link to DisGeNET
# Use gene_ids to find disease associations
```

## Related Resources

- [GMrepo](../gmrepo/README.md) - Sample-level microbiome data
- [gutMDisorder](../../9.3.microbe.host.interactions/gutmdisorder/README.md) - Disease associations
- [VMH](../../9.3.microbe.host.interactions/vmh/README.md) - Metabolic models
