---
id: download-hmp
title: "HMP Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Human Microbiome Project (HMP) Download Instructions

## Quick Start

```bash
# Download via AWS S3 (no account required)
aws s3 sync s3://human-microbiome-project/HMASM/ ./HMASM/ --no-sign-request
```

## Prerequisites

- **AWS CLI** for S3 downloads (recommended)
- **wget** or **curl** for direct downloads
- **Aspera** for high-speed transfers (optional)
- 50GB-50TB storage depending on data scope

## No Registration Required

HMP data is openly available under NIH data sharing policy.

## Download Methods

### Method 1: AWS S3 (Recommended for Large Downloads)

```bash
# Install AWS CLI if needed
pip install awscli

# List available data
aws s3 ls s3://human-microbiome-project/ --no-sign-request

# Download entire HMASM (processed data)
aws s3 sync s3://human-microbiome-project/HMASM/ ./HMASM/ --no-sign-request

# Download specific study (IBDMDB)
aws s3 sync s3://human-microbiome-project/HMASM/IBDMDB/ ./IBDMDB/ --no-sign-request

# Download only 16S data
aws s3 sync s3://human-microbiome-project/HMASM/IBDMDB/16S/ ./IBDMDB_16S/ --no-sign-request

# Download WGS data
aws s3 sync s3://human-microbiome-project/HMASM/IBDMDB/WGS/ ./IBDMDB_WGS/ --no-sign-request

# Download specific file types
aws s3 sync s3://human-microbiome-project/HMASM/ ./data/ \
  --exclude "*" --include "*.biom" --no-sign-request
```

### Method 2: Portal Data Browser

1. Visit https://portal.hmpdacc.org
2. Navigate to search: https://portal.hmpdacc.org/search/f
3. Apply filters:
   - Study (HMP1, IBDMDB, T2D, MOMS-PI)
   - Data type (16S, WGS, Proteomics, etc.)
   - Body site
4. Add files to cart
5. Download manifest or files directly

### Method 3: Direct HTTP Download

```bash
# Download from HMPDACC directly
wget https://downloads.hmpdacc.org/dacc/hmp1/otu_table.biom.gz

# Download specific dataset
wget -r -np -nH --cut-dirs=3 \
  https://downloads.hmpdacc.org/dacc/hmp1/16S/

# Download metadata
wget https://downloads.hmpdacc.org/dacc/metadata/hmp1_sample_metadata.tsv
```

### Method 4: Aspera High-Speed Transfer

```bash
# Install Aspera Connect
# Download from: https://www.ibm.com/aspera/connect/

# High-speed download
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -k1 -T -l 300m \
  fasp-beta@downloads.hmpdacc.org:/dacc/hmp1/16S ./hmp1_16S/
```

### Method 5: Using Manifest File

```bash
# Download manifest from portal cart
# Manifest contains file URLs

# Parse and download
while IFS=$'\t' read -r url md5 size filename; do
  echo "Downloading $filename..."
  wget -O "$filename" "$url"
  echo "$md5  $filename" | md5sum -c -
done < manifest.tsv
```

### Method 6: OSDF API (Programmatic)

```bash
# Query OSDF for specific data
curl "https://osdf.igs.umaryland.edu/api/query?node_type=16s_raw_seq_set&study=IBDMDB" \
  -o ibdmdb_16s_metadata.json

# Parse URLs from response
cat ibdmdb_16s_metadata.json | jq -r '.results[].urls[]' > download_urls.txt

# Download files
wget -i download_urls.txt
```

## File Inventory

### HMP1 (Healthy Cohort)

| Data Type | Size | Description |
|-----------|------|-------------|
| 16S profiles | ~5 GB | OTU tables, BIOM |
| WGS | ~10 TB | Shotgun sequences |
| Metadata | ~100 MB | Sample annotations |

### iHMP Studies

| Study | 16S | WGS | Other -omics | Total |
|-------|-----|-----|--------------|-------|
| IBDMDB | ~2 TB | ~5 TB | ~2 TB | ~9 TB |
| T2D | ~1 TB | ~3 TB | ~1 TB | ~5 TB |
| MOMS-PI | ~2 TB | ~4 TB | ~2 TB | ~8 TB |

### File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| FASTQ | .fastq.gz | Raw sequences |
| FASTA | .fasta.gz | Processed sequences |
| BIOM | .biom | OTU/ASV tables |
| TSV | .tsv | Metadata, abundance |
| JSON | .json | OSDF schema files |

## Post-Download Processing

```bash
# Decompress files
gunzip *.fastq.gz

# Convert BIOM to TSV
biom convert -i otu_table.biom -o otu_table.tsv --to-tsv

# Load into QIIME2
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path otu_table.biom \
  --output-path otu_table.qza

# Parse metadata
python3 << 'EOF'
import pandas as pd

# Load sample metadata
metadata = pd.read_csv("hmp1_sample_metadata.tsv", sep="\t")

# Filter by body site
stool = metadata[metadata["body_site"] == "stool"]
print(f"Stool samples: {len(stool)}")

# Get subject info
subjects = metadata.groupby("subject_id").agg({
    "sample_id": "count",
    "body_site": lambda x: ", ".join(sorted(set(x)))
})
subjects.to_csv("subject_summary.tsv", sep="\t")
EOF
```

### Process Multi-omics Data

```python
import pandas as pd
import os

# Load different data types
def load_ibdmdb_data(data_dir):
    """Load IBDMDB multi-omics data."""

    data = {}

    # 16S taxonomic profiles
    data['16s'] = pd.read_csv(
        os.path.join(data_dir, "16S", "taxonomic_profiles.tsv"),
        sep="\t", index_col=0
    )

    # WGS functional profiles
    data['wgs_func'] = pd.read_csv(
        os.path.join(data_dir, "WGS", "humann_pathabundance.tsv"),
        sep="\t", index_col=0
    )

    # Metabolomics
    data['metabolome'] = pd.read_csv(
        os.path.join(data_dir, "Metabolomics", "metabolites.tsv"),
        sep="\t", index_col=0
    )

    return data
```

## Verification

```bash
# Verify MD5 checksums
md5sum -c checksums.md5

# Check BIOM file
biom summarize-table -i otu_table.biom

# Validate FASTQ
seqkit stats sample.fastq.gz

# Check sample counts
aws s3 ls s3://human-microbiome-project/HMASM/IBDMDB/16S/ --no-sign-request | wc -l
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New data | As studies complete |
| Portal updates | Monthly |
| Major releases | Project milestones |

## Common Issues

- **Large data volume**: Use AWS S3 for resumable transfers
- **Storage space**: Plan for 10-50TB for complete iHMP
- **Network speed**: Consider Aspera for faster downloads
- **File organization**: Follow HMP directory structure
- **Metadata matching**: Use OSDF schema for linking files

## Study-Specific Notes

### IBDMDB (IBD Multi'omics)

```bash
# Download IBDMDB specifically
aws s3 sync s3://human-microbiome-project/HMASM/IBDMDB/ ./IBDMDB/ --no-sign-request

# Key data types
# - 16S rRNA profiles
# - Metagenomics (taxonomic + functional)
# - Metatranscriptomics
# - Proteomics
# - Metabolomics
# - Host transcriptomics
```

### T2D (Type 2 Diabetes)

```bash
# Download T2D study
aws s3 sync s3://human-microbiome-project/HMASM/T2D/ ./T2D/ --no-sign-request

# Longitudinal samples tracking diabetes progression
```

### MOMS-PI (Pregnancy)

```bash
# Download pregnancy study
aws s3 sync s3://human-microbiome-project/HMASM/MOMSPI/ ./MOMSPI/ --no-sign-request

# Multiple body sites during pregnancy
```

## Related Resources

- [GMrepo](../gmrepo/README.md) - Curated gut microbiome data
- [MetaHIT](../metahit/README.md) - European reference
- [HOMD](../../9.2.body.site.microbiomes/homd/README.md) - Oral microbiome
