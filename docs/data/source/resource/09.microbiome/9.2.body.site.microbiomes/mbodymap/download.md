---
id: download-mbodymap
title: "mBodyMap Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# mBodyMap Download Instructions

## Quick Start

```bash
# Access via web interface
# Visit https://mbodymap.microbiome.cloud for interactive data access
```

## Prerequisites

- Web browser for interactive access
- **Python** with requests for programmatic access
- 1-50GB storage depending on scope

## No Registration Required

mBodyMap data is freely available for academic use.

## Download Methods

### Method 1: Web Interface

1. Visit https://mbodymap.microbiome.cloud
2. Navigate to body site of interest
3. Select samples or studies
4. Export data in desired format (TSV, BIOM)

### Method 2: API Access (if available)

```python
import requests
import pandas as pd

# Check API documentation at the website
base_url = "https://mbodymap.microbiome.cloud/api"

# Example: Get samples by body site
def get_samples_by_site(site_code):
    url = f"{base_url}/samples"
    params = {"site": site_code}
    response = requests.get(url, params=params)
    return response.json()

# Example: Get taxonomic profile
def get_profile(sample_id):
    url = f"{base_url}/profile/{sample_id}"
    response = requests.get(url)
    return response.json()

# Fetch gut samples
gut_samples = get_samples_by_site("GI_LI_COLON")
print(f"Found {len(gut_samples)} colon samples")
```

### Method 3: Download Pre-computed Datasets

```bash
# Check Downloads section of website for bulk data
# Common files may include:

# Abundance matrices by body site
wget https://mbodymap.microbiome.cloud/data/gut_abundance_matrix.tsv
wget https://mbodymap.microbiome.cloud/data/oral_abundance_matrix.tsv
wget https://mbodymap.microbiome.cloud/data/skin_abundance_matrix.tsv

# Metadata files
wget https://mbodymap.microbiome.cloud/data/sample_metadata.tsv
wget https://mbodymap.microbiome.cloud/data/body_site_reference.tsv

# Note: Actual URLs depend on current website structure
```

### Method 4: Query-Based Export

```python
import requests
import json

# Build query for specific criteria
query = {
    "body_site": "GI_LI_COLON",
    "health_status": "healthy",
    "age_range": [18, 65],
    "sequencing_method": "16S"
}

# Submit query and download results
response = requests.post(
    "https://mbodymap.microbiome.cloud/api/query",
    json=query
)

results = response.json()
print(f"Found {len(results['samples'])} matching samples")

# Save results
with open("query_results.json", "w") as f:
    json.dump(results, f, indent=2)
```

### Method 5: Integration with Other Resources

```bash
# mBodyMap may link to original studies in SRA/ENA
# Download raw data from primary sources

# Example: Get HMP data linked from mBodyMap
aws s3 sync s3://human-microbiome-project/HMASM/ ./hmp_data/ --no-sign-request

# Example: Get specific study data from ENA
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456/ERR123456.fastq.gz
```

## File Inventory

### Expected Data Files

| File Type | Size (est.) | Description |
|-----------|-------------|-------------|
| Abundance matrices | 100 MB - 1 GB | Per body site |
| Sample metadata | 10-50 MB | All samples |
| Reference profiles | 5-20 MB | Site-specific |
| Diversity metrics | 20-100 MB | Pre-computed |

### Data Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| TSV | .tsv | General analysis |
| BIOM | .biom | QIIME/microbiome tools |
| JSON | .json | API responses |
| CSV | .csv | Spreadsheet analysis |

## Post-Download Processing

```python
import pandas as pd
import numpy as np

# Load abundance matrix
abundance = pd.read_csv("gut_abundance_matrix.tsv", sep="\t", index_col=0)
metadata = pd.read_csv("sample_metadata.tsv", sep="\t", index_col=0)

print(f"Samples: {abundance.shape[1]}")
print(f"Taxa: {abundance.shape[0]}")

# Filter by body site
colon_samples = metadata[metadata["site_code"] == "GI_LI_COLON"].index
colon_abundance = abundance[colon_samples]

# Calculate diversity
from scipy.stats import entropy

def shannon_diversity(abundances):
    """Calculate Shannon diversity."""
    p = abundances[abundances > 0]
    return entropy(p, base=2)

diversity = colon_abundance.apply(shannon_diversity)
print(f"Mean Shannon diversity: {diversity.mean():.2f}")

# Compare body sites
site_diversity = {}
for site in metadata["site_code"].unique():
    samples = metadata[metadata["site_code"] == site].index
    site_abundance = abundance[abundance.columns.intersection(samples)]
    site_diversity[site] = site_abundance.apply(shannon_diversity).mean()

site_df = pd.DataFrame.from_dict(site_diversity, orient="index", columns=["shannon"])
print("\nDiversity by body site:")
print(site_df.sort_values("shannon", ascending=False).head(10))
```

### Create Body Site Comparison

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
abundance = pd.read_csv("abundance_matrix.tsv", sep="\t", index_col=0)
metadata = pd.read_csv("sample_metadata.tsv", sep="\t", index_col=0)

# Calculate site-level means
site_means = {}
for site in metadata["site_code"].unique():
    samples = metadata[metadata["site_code"] == site].index
    site_abundance = abundance[abundance.columns.intersection(samples)]
    site_means[site] = site_abundance.mean(axis=1)

site_df = pd.DataFrame(site_means)

# Top taxa across sites
top_taxa = abundance.mean(axis=1).nlargest(20).index
heatmap_data = site_df.loc[top_taxa]

# Plot heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data, cmap="YlOrRd", annot=True, fmt=".2f")
plt.title("Top 20 Taxa Across Body Sites")
plt.tight_layout()
plt.savefig("body_site_comparison.png", dpi=150)
```

## Verification

```bash
# Check file integrity
wc -l gut_abundance_matrix.tsv

# Verify column counts match metadata
python3 << 'EOF'
import pandas as pd

abundance = pd.read_csv("gut_abundance_matrix.tsv", sep="\t", index_col=0)
metadata = pd.read_csv("sample_metadata.tsv", sep="\t", index_col=0)

# Check overlap
overlap = set(abundance.columns) & set(metadata.index)
print(f"Samples in abundance: {len(abundance.columns)}")
print(f"Samples in metadata: {len(metadata)}")
print(f"Overlapping: {len(overlap)}")
EOF
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New studies | Quarterly |
| Database refresh | Semi-annual |

## Common Issues

- **Large files**: Download in chunks if needed
- **Format variations**: Check delimiter and header rows
- **Sample ID mapping**: Match across different files
- **Missing data**: Some samples may lack full metadata
- **Version tracking**: Note database version for reproducibility

## Integration Examples

```python
# Link to HMP samples
def link_to_hmp(mbodymap_samples, hmp_metadata):
    """Find HMP equivalents for mBodyMap samples."""
    links = []
    for sample in mbodymap_samples:
        # Match by study/subject if available
        pass
    return links

# Link to GMrepo for gut samples
def link_to_gmrepo(sample_metadata):
    """Find GMrepo projects for gut samples."""
    gut_samples = sample_metadata[
        sample_metadata["site_code"].str.startswith("GI")
    ]
    # Query GMrepo API for matching projects
    pass

# Convert to QIIME2 format
def export_to_qiime(abundance_df, metadata_df, output_prefix):
    """Export for QIIME2 analysis."""
    # Feature table
    from biom import Table
    table = Table(abundance_df.values,
                  abundance_df.index.tolist(),
                  abundance_df.columns.tolist())
    with biom_open(f"{output_prefix}_table.biom", "w") as f:
        table.to_hdf5(f, "mBodyMap export")

    # Metadata
    metadata_df.to_csv(f"{output_prefix}_metadata.tsv", sep="\t")
```

## Related Resources

- [HMP](../hmp/README.md) - US reference microbiome data
- [HOMD](../homd/README.md) - Oral microbiome focus
- [GMrepo](../../9.1.gut.microbiome/gmrepo/README.md) - Gut microbiome repository
