---
id: download-gmrepo
title: "GMrepo Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# GMrepo Download Instructions

## Quick Start

```bash
# Fetch samples by phenotype via API
curl "https://gmrepo.humangut.info/api/getRunsByPhenotype?phenotype=Type%202%20diabetes" \
  -o t2d_samples.json
```

## Prerequisites

- **curl** or **wget** for downloads
- **Python** with requests library for batch access
- **jq** for JSON parsing (optional)
- 1-10GB storage for processed data

## No Registration Required

GMrepo data is freely available for academic use (CC BY-NC 3.0).

## Download Methods

### Method 1: REST API Queries

```bash
# Get samples by phenotype
curl "https://gmrepo.humangut.info/api/getRunsByPhenotype?phenotype=Type%202%20diabetes" \
  -o t2d_samples.json

# Get samples by project
curl "https://gmrepo.humangut.info/api/getRunsByProjectID?projectID=PRJEB6070" \
  -o project_samples.json

# Get species abundance for a sample
curl "https://gmrepo.humangut.info/api/getSpeciesAbundance?runID=SRR1234567" \
  -o sample_abundance.json

# Get marker taxa for a disease
curl "https://gmrepo.humangut.info/api/getMarkerTaxa?phenotype=Type%202%20diabetes" \
  -o t2d_markers.json

# Get all phenotypes
curl "https://gmrepo.humangut.info/api/getPhenotypeList" \
  -o phenotypes.json

# Get project list
curl "https://gmrepo.humangut.info/api/getProjectList" \
  -o projects.json

# Get statistics
curl "https://gmrepo.humangut.info/api/getStatistics" \
  -o stats.json
```

### Method 2: Python Batch Download

```python
import requests
import json
import time

base_url = "https://gmrepo.humangut.info/api"

def get_samples_by_phenotype(phenotype):
    """Fetch all samples for a phenotype."""
    url = f"{base_url}/getRunsByPhenotype"
    response = requests.get(url, params={"phenotype": phenotype})
    return response.json()

def get_abundance(run_id):
    """Get species abundance for a sample."""
    url = f"{base_url}/getSpeciesAbundance"
    response = requests.get(url, params={"runID": run_id})
    return response.json()

def get_marker_taxa(phenotype):
    """Get disease marker taxa."""
    url = f"{base_url}/getMarkerTaxa"
    response = requests.get(url, params={"phenotype": phenotype})
    return response.json()

# Example: Download T2D data
phenotype = "Type 2 diabetes"

# Get samples
samples = get_samples_by_phenotype(phenotype)
print(f"Found {len(samples)} samples for {phenotype}")

# Save samples
with open("t2d_samples.json", "w") as f:
    json.dump(samples, f, indent=2)

# Get marker taxa
markers = get_marker_taxa(phenotype)
with open("t2d_markers.json", "w") as f:
    json.dump(markers, f, indent=2)

# Get abundance for first 10 samples
for sample in samples[:10]:
    run_id = sample["run_id"]
    abundance = get_abundance(run_id)
    with open(f"abundance_{run_id}.json", "w") as f:
        json.dump(abundance, f, indent=2)
    time.sleep(0.5)  # Be polite
    print(f"Downloaded {run_id}")
```

### Method 3: Bulk Data Download (Processed Tables)

```bash
# Download processed abundance tables (if available)
# Check Downloads page for bulk data
curl "https://gmrepo.humangut.info/Downloads" -o downloads_page.html

# Download from the downloads page (URLs may change)
# wget "https://gmrepo.humangut.info/data/abundance_table.tsv.gz"
```

### Method 4: Web Interface Export

1. Visit https://gmrepo.humangut.info/search/f
2. Apply filters (phenotype, project, organism)
3. Select samples of interest
4. Export to CSV/TSV

### Method 5: Programmatic Access to Raw Sequences

```bash
# GMrepo links to SRA/ENA for raw sequences
# Get sample metadata including SRA accessions
curl "https://gmrepo.humangut.info/api/getRunsByPhenotype?phenotype=Healthy" \
  | jq -r '.[].run_id' > sra_accessions.txt

# Download raw sequences from SRA (requires SRA toolkit)
while read acc; do
  fastq-dump --split-files $acc
done < sra_accessions.txt
```

## File Inventory

### API Response Files

| Endpoint | Size | Description |
|----------|------|-------------|
| getPhenotypeList | ~50 KB | All phenotypes |
| getProjectList | ~500 KB | All projects |
| getRunsByPhenotype | Variable | Samples per disease |
| getSpeciesAbundance | ~10-50 KB | Per-sample abundance |
| getMarkerTaxa | ~100 KB | Disease markers |

### Estimated Data Sizes

| Data Type | Size | Description |
|-----------|------|-------------|
| Sample metadata | ~50 MB | All 118K samples |
| Abundance profiles | ~5-10 GB | All samples |
| Marker taxa | ~10 MB | All phenotypes |
| Raw sequences | ~100-400 GB | Via SRA links |

## Post-Download Processing

```bash
# Parse JSON to TSV
python3 << 'EOF'
import json
import csv

# Convert samples to TSV
with open("t2d_samples.json") as f:
    samples = json.load(f)

with open("t2d_samples.tsv", "w", newline="") as f:
    if samples:
        writer = csv.DictWriter(f, fieldnames=samples[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(samples)
EOF

# Parse abundance data
python3 << 'EOF'
import json
import pandas as pd
import glob

# Combine abundance files into matrix
all_abundances = {}
for filename in glob.glob("abundance_*.json"):
    with open(filename) as f:
        data = json.load(f)
        run_id = filename.replace("abundance_", "").replace(".json", "")
        for species in data.get("species", []):
            taxon = species["taxon_name"]
            abundance = species["relative_abundance"]
            if taxon not in all_abundances:
                all_abundances[taxon] = {}
            all_abundances[taxon][run_id] = abundance

df = pd.DataFrame(all_abundances).T.fillna(0)
df.to_csv("abundance_matrix.tsv", sep="\t")
print(f"Created matrix: {df.shape[0]} taxa x {df.shape[1]} samples")
EOF

# Create BIOM format (for QIIME compatibility)
python3 << 'EOF'
import pandas as pd
from biom import Table
from biom.util import biom_open

df = pd.read_csv("abundance_matrix.tsv", sep="\t", index_col=0)
table = Table(df.values, df.index.tolist(), df.columns.tolist())

with biom_open("abundance.biom", "w") as f:
    table.to_hdf5(f, "gmrepo_export")
EOF
```

## Verification

```bash
# Check JSON validity
cat t2d_samples.json | jq . > /dev/null && echo "Valid JSON"

# Count samples
cat t2d_samples.json | jq 'length'

# Check phenotypes
cat phenotypes.json | jq '.[].phenotype' | head -20

# Verify marker taxa structure
cat t2d_markers.json | jq 'keys'
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New projects | Continuous |
| Database refresh | Quarterly |
| Version releases | Annually |

## Common Issues

- **Rate limiting**: Add delays between API calls
- **Large responses**: Some phenotypes have thousands of samples
- **JSON parsing**: Use streaming for large files
- **Raw sequences**: Not hosted by GMrepo; access via SRA/ENA
- **Commercial use**: Requires permission (CC BY-NC 3.0)

## API Rate Limits

```python
import time

# Recommended: 1 request per second
def rate_limited_request(url, params=None):
    response = requests.get(url, params=params)
    time.sleep(1)  # Wait 1 second between requests
    return response
```

## Integration with SRA

```bash
# GMrepo provides links to raw sequences in SRA
# Download raw sequences using SRA toolkit

# Install SRA toolkit
# conda install -c bioconda sra-tools

# Download FASTQ files
fastq-dump --split-files --gzip SRR1234567

# Or use faster fasterq-dump
fasterq-dump SRR1234567
```

## Related Resources

- [HMP](../../9.2.body.site.microbiomes/hmp/_index.md) - US reference microbiomes
- [MetaHIT](../metahit/_index.md) - European metagenomes
- [gutMGene](../gutmgene/_index.md) - Microbe-gene links
