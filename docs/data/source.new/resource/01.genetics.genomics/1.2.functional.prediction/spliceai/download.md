---
id: download-spliceai
title: "SpliceAI Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# SpliceAI Download Instructions

## Quick Start

```bash
# Download pre-computed scores from Illumina
wget https://basespace.illumina.com/s/otSPW8hnhaZR/download -O spliceai_scores.raw.snv.hg38.vcf.gz
```

## Prerequisites

- wget or web browser (BaseSpace account may be required)
- ~20 GB disk space
- tabix for indexed queries
- bcftools (optional)

## Download Methods

### Primary: Illumina BaseSpace

```bash
# Navigate to Illumina BaseSpace and download
# https://basespace.illumina.com/s/otSPW8hnhaZR

# Files available:
# - spliceai_scores.raw.snv.hg38.vcf.gz
# - spliceai_scores.raw.snv.hg19.vcf.gz
# - spliceai_scores.raw.indel.hg38.vcf.gz
# - spliceai_scores.raw.indel.hg19.vcf.gz
```

### Alternative: Broad Institute Lookup API

```bash
# Query individual variants via API
curl "https://spliceailookup-api.broadinstitute.org/spliceai/?hg=38&variant=2-179446218-A-G"
```

### Alternative: GitHub (Model Only)

```bash
# Clone model for scoring novel variants
git clone https://github.com/Illumina/SpliceAI.git
pip install spliceai
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| spliceai_scores.raw.snv.hg38.vcf.gz | ~15 GB | SNVs, GRCh38 |
| spliceai_scores.raw.snv.hg19.vcf.gz | ~15 GB | SNVs, GRCh37 |
| spliceai_scores.raw.indel.hg38.vcf.gz | ~2 GB | Indels, GRCh38 |
| spliceai_scores.raw.indel.hg19.vcf.gz | ~2 GB | Indels, GRCh37 |

## Post-Download Processing

```bash
# Create tabix index
tabix -p vcf spliceai_scores.raw.snv.hg38.vcf.gz

# Query specific region
tabix spliceai_scores.raw.snv.hg38.vcf.gz 2:179446000-179447000

# Extract high-scoring variants
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SpliceAI\n' \
  spliceai_scores.raw.snv.hg38.vcf.gz | \
  awk -F'|' '$3>=0.5 || $4>=0.5 || $5>=0.5 || $6>=0.5' > high_spliceai.tsv
```

## Verification

```bash
# Check VCF header
bcftools view -h spliceai_scores.raw.snv.hg38.vcf.gz | head

# Test tabix query
tabix spliceai_scores.raw.snv.hg38.vcf.gz 1:100000-100100
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major version | Infrequent |
| Model updates | As published |

## Integration Notes

- VEP Plugin available: `--plugin SpliceAI,snv=spliceai_scores.raw.snv.hg38.vcf.gz,indel=spliceai_scores.raw.indel.hg38.vcf.gz`
- Integrated in dbNSFP
- Can run model locally for novel variants: `spliceai -I input.vcf -O output.vcf -R genome.fa -A grch38`
