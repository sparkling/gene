---
id: download-gnomad
title: "gnomAD Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# gnomAD Download Instructions

## Quick Start

```bash
# Download gnomAD v4.1 exomes VCF (chr1) - requires gsutil
gsutil -m cp gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr1.vcf.bgz* .
```

## Prerequisites

- **gsutil** (Google Cloud SDK) or **AWS CLI** for cloud downloads
- **bcftools** or **tabix** for VCF indexing and subsetting
- **Hail** (optional) for Hail Table format
- Minimum 500GB-2TB disk space depending on data scope

## No Registration Required

gnomAD data is publicly available without registration. Data is released under ODC-ODbL 1.0 license.

## Download Methods

### Method 1: Google Cloud Storage (Recommended)

```bash
# Install gsutil if needed
pip install gsutil

# List available files
gsutil ls gs://gcp-public-data--gnomad/release/

# Download gnomAD v4.1 exomes (sites-only VCF, all chromosomes)
for chr in {1..22} X Y; do
  gsutil -m cp \
    gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz \
    gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz.tbi \
    ./gnomad_exomes/
done

# Download gnomAD v4.1 genomes
for chr in {1..22} X Y; do
  gsutil -m cp \
    gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz \
    gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz.tbi \
    ./gnomad_genomes/
done
```

### Method 2: AWS Open Data

```bash
# Using AWS CLI (no credentials required for public data)
aws s3 ls --no-sign-request s3://gnomad-public-us-east-1/release/

# Download specific files
aws s3 cp --no-sign-request \
  s3://gnomad-public-us-east-1/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr1.vcf.bgz \
  ./gnomad_exomes/
```

### Method 3: Hail Tables (For Large-Scale Analysis)

```python
import hail as hl

# Initialize Hail
hl.init()

# Load gnomAD v4 Hail Table directly from cloud
gnomad_ht = hl.read_table(
    'gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht'
)

# Filter and export subset
gnomad_ht.filter(gnomad_ht.freq[0].AF > 0.01).export('common_variants.tsv')
```

### Method 4: Download Portal

Navigate to https://gnomad.broadinstitute.org/downloads for interactive file selection.

## File Inventory

### gnomAD v4.1 Exomes

| File Pattern | Size (per chr) | Description |
|--------------|----------------|-------------|
| gnomad.exomes.v4.1.sites.chr*.vcf.bgz | 2-15 GB | Sites-only VCF |
| gnomad.exomes.v4.1.sites.chr*.vcf.bgz.tbi | 2-5 MB | Tabix index |
| gnomad.exomes.v4.1.sites.ht | ~250 GB total | Hail Table format |

### gnomAD v4.1 Genomes

| File Pattern | Size (per chr) | Description |
|--------------|----------------|-------------|
| gnomad.genomes.v4.1.sites.chr*.vcf.bgz | 5-30 GB | Sites-only VCF |
| gnomad.genomes.v4.1.sites.chr*.vcf.bgz.tbi | 3-8 MB | Tabix index |
| gnomad.genomes.v4.1.sites.ht | ~400 GB total | Hail Table format |

### Supplementary Files

| File | Size | Description |
|------|------|-------------|
| gnomad.v4.1.constraint_metrics.tsv | ~50 MB | Gene constraint scores (pLI, LOEUF) |
| gnomad.v4.1.sv.sites.vcf.bgz | ~5 GB | Structural variants |
| gnomad.v4.1.coverage.summary.tsv.bgz | ~100 MB | Coverage statistics |

## Post-Download Processing

```bash
# Verify file integrity
md5sum -c gnomad.exomes.v4.1.sites.chr1.vcf.bgz.md5

# Index VCF if needed
tabix -p vcf gnomad.exomes.v4.1.sites.chr1.vcf.bgz

# Extract specific region
bcftools view -r chr17:43044295-43125483 \
  gnomad.exomes.v4.1.sites.chr17.vcf.bgz \
  -o brca1_variants.vcf.gz -Oz

# Convert to TSV for simpler processing
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' \
  gnomad.exomes.v4.1.sites.chr1.vcf.bgz > chr1_allele_freq.tsv
```

## Verification

```bash
# Check VCF header
bcftools view -h gnomad.exomes.v4.1.sites.chr1.vcf.bgz | head -50

# Count variants
bcftools view -H gnomad.exomes.v4.1.sites.chr1.vcf.bgz | wc -l

# Verify index
tabix -l gnomad.exomes.v4.1.sites.chr1.vcf.bgz
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major versions (v4, v5) | Every 2-3 years |
| Minor updates (v4.1, v4.2) | Every 6-12 months |
| Constraint metrics | Updated with major releases |

## Common Issues

- **Insufficient disk space**: gnomAD full download requires 1-2TB; download chromosome by chromosome
- **Slow download speeds**: Use regional mirrors (us-east-1 for US, europe-west1 for EU)
- **Memory errors with Hail**: Increase executor memory or use cloud-based Hail clusters
- **Missing index files**: Always download .tbi files alongside VCFs
- **Version mismatch**: Ensure VCF version matches your reference genome (GRCh38 for v4+)

## Storage Optimization

```bash
# Download only specific populations (using bcftools after download)
bcftools view -i 'INFO/AF_nfe > 0' gnomad.exomes.v4.1.sites.chr1.vcf.bgz \
  -o nfe_variants.vcf.gz -Oz

# Keep only common variants
bcftools view -i 'INFO/AF > 0.01' gnomad.exomes.v4.1.sites.chr1.vcf.bgz \
  -o common_variants.vcf.gz -Oz
```
