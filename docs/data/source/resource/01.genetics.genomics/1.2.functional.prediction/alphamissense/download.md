---
id: download-alphamissense
title: "AlphaMissense Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# AlphaMissense Download Instructions

## Quick Start

```bash
# Download hg38 predictions from Google Cloud
gsutil cp gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz .
```

## Prerequisites

- gsutil (Google Cloud SDK) or wget
- ~700 MB disk space for hg38 file
- gzip for decompression

## Download Methods

### Primary: Google Cloud Storage

```bash
# List available files
gsutil ls gs://dm_alphamissense/

# Download GRCh38 predictions
gsutil cp gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz .

# Download GRCh37 predictions (if needed)
gsutil cp gs://dm_alphamissense/AlphaMissense_hg19.tsv.gz .

# Download gene-level summary
gsutil cp gs://dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz .
```

### Alternative: Zenodo

```bash
# Download from Zenodo mirror
wget https://zenodo.org/records/8360242/files/AlphaMissense_hg38.tsv.gz
```

### Alternative: Direct URL

```bash
# Using curl
curl -O https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| AlphaMissense_hg38.tsv.gz | 643 MB | All variants, GRCh38 |
| AlphaMissense_hg19.tsv.gz | 643 MB | All variants, GRCh37 |
| AlphaMissense_gene_hg38.tsv.gz | 2 MB | Gene-level summary |
| AlphaMissense_aa_substitutions.tsv.gz | 6.9 GB | All AA substitutions |

## Post-Download Processing

```bash
# Decompress
gunzip AlphaMissense_hg38.tsv.gz

# View header
head -1 AlphaMissense_hg38.tsv

# Count variants
wc -l AlphaMissense_hg38.tsv

# Filter for a specific gene (example: BRAF)
grep "BRAF" AlphaMissense_hg38.tsv > braf_variants.tsv
```

## Verification

```bash
# Check file integrity
md5sum AlphaMissense_hg38.tsv.gz

# Verify expected variant count
zcat AlphaMissense_hg38.tsv.gz | wc -l
# Expected: ~71,000,000 lines
```

## Dataset Versions

### Current Release: v1.0.0

| Property | Value |
|----------|-------|
| Version | 1.0.0 |
| Release Date | 2023-09-22 |
| Total Size | ~7.5 GB (all files) |
| Variants | 71M missense predictions |
| Genes | 19,233 protein-coding |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| AlphaMissense_hg38.tsv.gz | 643 MB | 71M | All variants (GRCh38) |
| AlphaMissense_hg19.tsv.gz | 643 MB | 71M | All variants (GRCh37) |
| AlphaMissense_gene_hg38.tsv.gz | 2 MB | 19K | Gene-level summary |
| AlphaMissense_aa_substitutions.tsv.gz | 6.9 GB | 216M | All AA substitutions |

### Classification Thresholds

| Classification | Score Range | Count |
|----------------|-------------|-------|
| Likely Pathogenic | >= 0.564 | 32% |
| Ambiguous | 0.340 - 0.564 | 23% |
| Likely Benign | < 0.340 | 45% |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | gs://dm_alphamissense/ (GCS) |
| Rate Limit | N/A (bulk download) |
| Auth Required | No |
| Response Format | TSV (gzipped) |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major version | As published |
| Data updates | N/A (static dataset) |

## Integration Notes

- VEP Plugin available: `--plugin AlphaMissense,file=AlphaMissense_hg38.tsv.gz`
- Pre-integrated in dbNSFP v5.0+
- Can be loaded into annotation databases (Annovar, SnpEff)
