---
id: download-hmp
title: "HMP Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# HMP Download Instructions (Body Site Microbiomes)

## Quick Start

```bash
# Download all body site data via AWS S3 (no authentication required)
aws s3 sync s3://human-microbiome-project/ ./hmp/ --no-sign-request
```

## Prerequisites

- **AWS CLI** (recommended for large downloads)
- **wget** or **curl** for direct HTTP downloads
- **Aspera** (optional, for high-speed transfers)
- Storage: 50GB-50TB depending on data scope

## No Registration Required

HMP data is openly available under NIH public domain data sharing policy. No account or data use agreement required for download.

---

## Download Methods

### Method 1: AWS S3 (Recommended)

```bash
# Install AWS CLI if needed
pip install awscli

# List available data
aws s3 ls s3://human-microbiome-project/ --no-sign-request

# Download complete dataset
aws s3 sync s3://human-microbiome-project/ ./hmp/ --no-sign-request

# Download specific body site data (reference set)
aws s3 sync s3://human-microbiome-project/HMASM/ ./HMASM/ --no-sign-request

# Download reference genomes
aws s3 sync s3://human-microbiome-project/HMP_REFERENCE_GENOMES/ ./genomes/ --no-sign-request

# Filter by file type (16S abundance tables)
aws s3 sync s3://human-microbiome-project/HMASM/ ./data/ \
  --exclude "*" --include "*.biom" --no-sign-request

# Download specific study
aws s3 sync s3://human-microbiome-project/HMASM/HMP1/ ./HMP1/ --no-sign-request
```

### Method 2: HMP Portal Data Browser

1. Visit https://portal.hmpdacc.org
2. Navigate to search: https://portal.hmpdacc.org/search/f
3. Apply body site filters:
   - **Supersite**: Airways, Gastrointestinal, Oral, Skin, Urogenital
   - **Body Site**: Specific anatomical location
   - **Data Type**: 16S, WGS, Proteomics, Metabolomics, etc.
4. Add files to cart
5. Download manifest or files directly

### Method 3: Direct HTTP Download

```bash
# Download from HMPDACC directly
wget https://downloads.hmpdacc.org/dacc/hmp1/otu_table.biom.gz

# Download metadata
wget https://downloads.hmpdacc.org/dacc/metadata/hmp1_sample_metadata.tsv

# Download body site-specific data
wget -r -np -nH --cut-dirs=3 \
  https://downloads.hmpdacc.org/dacc/hmp1/oral/

# Download reference genomes
wget -r -np -nH --cut-dirs=3 \
  https://downloads.hmpdacc.org/reference_genomes/
```

### Method 4: Aspera High-Speed Transfer

```bash
# Install Aspera Connect from IBM
# https://www.ibm.com/aspera/connect/

# High-speed download
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -k1 -T -l 300m \
  fasp-beta@downloads.hmpdacc.org:/dacc/hmp1/ ./hmp1/
```

### Method 5: NCBI SRA Access

```bash
# Using SRA Toolkit
# Install: https://github.com/ncbi/sra-tools

# Download specific sample by SRA accession
prefetch SRR346657

# Convert to FASTQ
fastq-dump --split-files SRR346657

# Bulk download by BioProject
prefetch --option-file hmp_sra_accessions.txt
```

### Method 6: OSDF API (Programmatic)

```bash
# Query specific body site samples
curl "https://osdf.igs.umaryland.edu/api/query?body_site=Tongue%20dorsum&node_type=sample" \
  -o tongue_samples.json

# Get file URLs
cat tongue_samples.json | jq -r '.results[].urls[]' > download_urls.txt

# Download files
wget -i download_urls.txt
```

---

## File Inventory by Body Site

### HMP1 Reference Set

| Body Site Category | Sample Count | 16S Data | WGS Data |
|--------------------|--------------|----------|----------|
| **Oral** (9 sites) | ~4,000 | ~3 TB | ~5 TB |
| **Stool** | ~3,500 | ~2 TB | ~3 TB |
| **Skin** (4 sites) | ~1,000 | ~800 GB | ~1 TB |
| **Nasal** | ~2,000 | ~1 TB | ~1.5 TB |
| **Vaginal** (3 sites) | ~500 | ~400 GB | ~600 GB |

### iHMP Studies (Multi-site Longitudinal)

| Study | Body Sites | Samples | Total Size |
|-------|------------|---------|------------|
| IBDMDB | Stool, Biopsy | ~2,000 | ~9 TB |
| T2D | Stool, Nasal | ~1,500 | ~5 TB |
| MOMS-PI | Vaginal, Oral, Skin | ~3,000 | ~8 TB |

### File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| FASTQ | .fastq.gz | Raw sequences |
| FASTA | .fasta.gz | Processed sequences |
| BIOM | .biom | OTU/ASV tables |
| TSV | .tsv | Metadata, abundance |
| JSON | .json | OSDF metadata |

---

## Body Site-Specific Downloads

### Oral Microbiome Data

```bash
# All oral samples (9 subsites)
aws s3 sync s3://human-microbiome-project/HMASM/ ./oral/ \
  --exclude "*" --include "*oral*" --no-sign-request

# Specific oral sites
# Tongue dorsum, Buccal mucosa, Hard palate, etc.
curl "https://portal.hmpdacc.org/api/search?supersite=Oral" \
  -o oral_manifest.json
```

### Skin Microbiome Data

```bash
# All skin samples (4 subsites)
aws s3 sync s3://human-microbiome-project/HMASM/ ./skin/ \
  --exclude "*" --include "*skin*" --no-sign-request

# Retroauricular crease, Antecubital fossa
```

### Nasal/Airways Data

```bash
# Anterior nares samples
curl "https://portal.hmpdacc.org/api/search?body_site=Anterior%20nares" \
  -o nasal_manifest.json
```

### Vaginal Microbiome Data

```bash
# All vaginal subsites
aws s3 sync s3://human-microbiome-project/HMASM/ ./vaginal/ \
  --exclude "*" --include "*vaginal*" --no-sign-request

# Posterior fornix, Mid vagina, Vaginal introitus
```

---

## Post-Download Processing

### Decompress and Organize

```bash
# Decompress files
gunzip -r ./hmp/

# Organize by body site
mkdir -p organized/{oral,skin,nasal,gut,vaginal}
# Use metadata to sort files
```

### Convert Formats

```bash
# BIOM to TSV
biom convert -i otu_table.biom -o otu_table.tsv --to-tsv

# Load into QIIME2
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path otu_table.biom \
  --output-path otu_table.qza
```

### Filter by Body Site

```python
import pandas as pd

# Load sample metadata
metadata = pd.read_csv("hmp_sample_metadata.tsv", sep="\t")

# Filter by supersite
oral_samples = metadata[metadata["supersite"] == "Oral"]
skin_samples = metadata[metadata["supersite"] == "Skin"]

# Filter by specific body site
tongue = metadata[metadata["body_site"] == "Tongue dorsum"]

# Get sample counts per site
site_counts = metadata.groupby(["supersite", "body_site"]).size()
print(site_counts)
```

### Cross-Site Analysis Setup

```python
import pandas as pd
from biom import load_table

# Load abundance matrix
otu_table = load_table("hmp_otu_table.biom")

# Load metadata with body site annotations
metadata = pd.read_csv("hmp_sample_metadata.tsv", sep="\t", index_col=0)

# Subset by body site for comparison
def get_site_samples(otu_table, metadata, body_site):
    """Get OTU table for specific body site."""
    site_samples = metadata[metadata["body_site"] == body_site].index
    return otu_table.filter(site_samples, axis="sample")

# Get body site-specific tables
tongue_otu = get_site_samples(otu_table, metadata, "Tongue dorsum")
stool_otu = get_site_samples(otu_table, metadata, "Stool")
skin_otu = get_site_samples(otu_table, metadata, "Retroauricular crease")
```

---

## Verification

```bash
# Verify MD5 checksums
md5sum -c checksums.md5

# Check BIOM file integrity
biom validate-table -i otu_table.biom

# Verify FASTQ
seqkit stats sample.fastq.gz

# Check download completeness
aws s3 ls s3://human-microbiome-project/HMASM/ --no-sign-request --recursive | wc -l
```

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Data corrections | As needed |
| Portal updates | Monthly |
| Major releases | Project milestones |
| New iHMP data | As studies complete |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Large data volume | Use AWS S3 sync with `--exclude` filters |
| Slow download | Use Aspera for faster transfers |
| Storage space | Plan for 10-50TB for complete HMP |
| Body site matching | Use OSDF metadata for sample-file linking |
| Cross-site comparison | Ensure consistent processing pipeline |

---

## Cloud Analysis Options

### AWS

```bash
# Data already on AWS Open Data Registry
# Launch EC2 in same region (us-east-1) for free data transfer
# Use s3:// paths directly in analysis tools
```

### Google Cloud

```bash
# Transfer to GCS
gsutil -m cp -r s3://human-microbiome-project/ gs://your-bucket/hmp/
```

### Terra/FireCloud

HMP data available as workspace datasets for cloud-based analysis.

---

## Related Resources

- [HMP Portal](https://portal.hmpdacc.org) - Main data access
- [AWS Open Data](https://registry.opendata.aws/human-microbiome-project/) - Cloud access
- [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA28331) - Project metadata
- [HOMD](../homd/README.md) - Oral microbiome detail
- [mBodyMap](../mbodymap/README.md) - Body site atlas
- [HMP in Gut Microbiome](../../9.1.gut.microbiome/hmp/README.md) - Stool-focused data
