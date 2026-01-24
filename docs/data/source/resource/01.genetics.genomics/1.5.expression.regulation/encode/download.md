---
id: download-encode
title: "ENCODE Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# ENCODE Download Instructions

## Quick Start

```bash
# Download specific experiment files using ENCODE portal
curl -O "https://www.encodeproject.org/files/ENCFF123ABC/@@download/ENCFF123ABC.bam"

# Or use the metadata-based batch download
xargs -L 1 curl -O -J -L < files.txt
```

## Prerequisites

- **curl** or **wget** for downloads
- **jq** for JSON parsing (optional)
- **samtools** for BAM processing
- **bedtools** for BED/BigBed manipulation
- 1-100TB storage depending on data scope

## No Registration Required

ENCODE data is freely available under ENCODE Data Use Policy.

## Download Methods

### Method 1: ENCODE Portal Search + Download

1. Navigate to https://www.encodeproject.org
2. Use search filters (assay type, biosample, etc.)
3. Select experiments/files
4. Click "Download" to get file manifest

```bash
# Download files from manifest
cut -f1 files.txt | xargs -I {} curl -O -L "https://www.encodeproject.org/files/{}/@@download/{}.{ext}"
```

### Method 2: REST API Batch Download

```bash
# Search for ChIP-seq experiments in K562
curl -s "https://www.encodeproject.org/search/?type=Experiment\
&assay_title=ChIP-seq&biosample_ontology.term_name=K562\
&format=json&limit=all" | jq '.["@graph"][] | .accession'

# Get file list for an experiment
curl -s "https://www.encodeproject.org/experiments/ENCSR000EAD/?format=json" | \
  jq -r '.files[] | select(.file_format=="bam") | .href'

# Batch download all files from search
curl -s "https://www.encodeproject.org/batch_download/?type=Experiment\
&assay_title=ChIP-seq&biosample_ontology.term_name=K562" > files.txt
xargs -L 1 curl -O -J -L < files.txt
```

### Method 3: Metadata TSV Download

```bash
# Get metadata for all experiments
curl -O "https://www.encodeproject.org/metadata/?type=Experiment&format=tsv"

# Filter and download specific files
awk -F'\t' '$1=="ChIP-seq" && $2=="GRCh38" {print $NF}' metadata.tsv | \
  xargs -I {} curl -O -L {}
```

### Method 4: AWS S3 Mirror (Fastest for Large Downloads)

```bash
# ENCODE data is mirrored on AWS Open Data
# Use AWS CLI (no credentials required)
aws s3 ls --no-sign-request s3://encode-public/

# Download specific file
aws s3 cp --no-sign-request \
  s3://encode-public/2023/01/01/ENCFF123ABC.bam ./

# Sync entire dataset
aws s3 sync --no-sign-request \
  s3://encode-public/2023/01/01/ ./encode_data/ \
  --exclude "*" --include "*.bigWig"
```

### Method 5: Specific Data Type Downloads

```bash
# ChIP-seq peaks (BED narrowPeak)
curl -s "https://www.encodeproject.org/search/?type=File\
&file_format=bed&output_type=peaks&assay_title=ChIP-seq\
&assembly=GRCh38&format=json&limit=100" | \
  jq -r '.["@graph"][] | .href' | \
  xargs -L 1 curl -O -L

# RNA-seq gene quantifications
curl -s "https://www.encodeproject.org/search/?type=File\
&file_format=tsv&output_type=gene+quantifications\
&assay_title=RNA-seq&format=json" | \
  jq -r '.["@graph"][] | .href' | head -10 | \
  xargs -L 1 curl -O -J -L

# ATAC-seq signal tracks
curl -s "https://www.encodeproject.org/search/?type=File\
&file_format=bigWig&assay_title=ATAC-seq\
&assembly=GRCh38&format=json" | \
  jq -r '.["@graph"][] | .href' | \
  xargs -L 1 curl -O -J -L
```

## File Inventory

### By Assay Type

| Assay | Typical Files | Size Range |
|-------|---------------|------------|
| ChIP-seq | BAM, bigWig, narrowPeak | 1-50 GB/experiment |
| RNA-seq | BAM, bigWig, TSV quantifications | 5-100 GB/experiment |
| ATAC-seq | BAM, bigWig, peaks | 2-30 GB/experiment |
| Hi-C | PAIRS, HIC, cool | 50-500 GB/experiment |
| WGBS | BAM, bigWig, BED | 50-200 GB/experiment |

### By File Format

| Format | Description | Tools |
|--------|-------------|-------|
| BAM | Aligned reads | samtools |
| bigWig | Signal tracks | UCSC tools |
| narrowPeak | ChIP-seq peaks | bedtools |
| broadPeak | Broad peaks | bedtools |
| TSV | Quantifications | awk, pandas |
| HIC | Hi-C matrices | juicer |

### Annotation Files

| File | Description |
|------|-------------|
| ENCFF762MJQ | GRCh38 gene annotations |
| ENCFF285DRD | GENCODE v29 GTF |
| Blacklist | ENCODE blacklist regions |

## Post-Download Processing

```bash
# Index BAM files
samtools index ENCFF123ABC.bam

# Convert bigWig to bedGraph
bigWigToBedGraph ENCFF123ABC.bigWig ENCFF123ABC.bedGraph

# Merge peaks from replicates
bedtools intersect -a rep1_peaks.bed -b rep2_peaks.bed > reproducible_peaks.bed

# Extract reads from specific region
samtools view -b ENCFF123ABC.bam chr17:43044295-43125483 > brca1_region.bam
```

## Verification

```bash
# Check BAM integrity
samtools quickcheck ENCFF123ABC.bam

# Verify file against ENCODE metadata
curl -s "https://www.encodeproject.org/files/ENCFF123ABC/?format=json" | \
  jq '{md5sum, file_size}'

# Compare local checksum
md5sum ENCFF123ABC.bam
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New experiments | Rolling (daily) |
| Assembly updates | Major releases |
| Reprocessing | Periodic pipeline updates |

## Common Issues

- **Large downloads**: Use AWS mirror for faster transfers
- **Paired files**: Some analyses require matching files (e.g., R1/R2 for fastq)
- **Replicate handling**: Check biological vs technical replicates
- **Assembly versions**: Ensure GRCh38 vs hg38 naming consistency
- **Status filtering**: Filter for "released" status; avoid "revoked" files

## ENCODE File Statuses

| Status | Description |
|--------|-------------|
| released | Publicly available |
| archived | Superseded version |
| revoked | Quality issues |
| in progress | Not yet released |

## Quality Metrics

```bash
# Get quality metrics for experiment
curl -s "https://www.encodeproject.org/experiments/ENCSR000EAD/?format=json" | \
  jq '.audit'

# Filter for high-quality data
curl -s "https://www.encodeproject.org/search/?type=Experiment\
&audit.NOT_COMPLIANT.path=*&audit.ERROR.path=*"
```

## Recommended Downloads by Analysis

### For ChIP-seq Analysis
- Fold-change over control (bigWig)
- Optimal IDR peaks (narrowPeak)
- Alignment files (BAM)

### For RNA-seq Analysis
- Gene quantifications (TSV)
- Transcript quantifications (TSV)
- Signal tracks (bigWig)

### For Chromatin Accessibility
- ATAC-seq peaks (narrowPeak)
- DNase-seq peaks (narrowPeak)
- Signal tracks (bigWig)
