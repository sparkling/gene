---
id: download-cadd
title: "CADD Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# CADD Download Instructions

## Quick Start

```bash
# Download pre-scored SNVs for GRCh38 (large file ~80GB)
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```

## Prerequisites

- wget or curl
- ~100 GB disk space for full genome
- tabix for indexed queries
- bgzip (optional, for compression)

## Download Methods

### Primary: CADD Download Server

```bash
# GRCh38 pre-scored SNVs (all possible)
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi

# GRCh37 pre-scored SNVs
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi

# Pre-scored indels (scored from gnomAD, ClinVar, etc.)
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz
```

### Alternative: Web Scoring

```bash
# For novel variants or small batches, use web scoring
# https://cadd.gs.washington.edu/score
# Upload VCF or paste variants
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| whole_genome_SNVs.tsv.gz | ~80 GB | All possible SNVs |
| whole_genome_SNVs.tsv.gz.tbi | ~3 MB | Tabix index |
| gnomad.genomes.r4.0.indel.tsv.gz | ~2 GB | Pre-scored indels |
| InDels_inclAnno.tsv.gz | Varies | Annotated indels |

## Post-Download Processing

```bash
# Query specific region using tabix
tabix whole_genome_SNVs.tsv.gz 7:140753336-140753336

# Query gene region (BRAF)
tabix whole_genome_SNVs.tsv.gz 7:140719327-140924929 > braf_cadd.tsv

# Extract high-scoring variants
zcat whole_genome_SNVs.tsv.gz | awk '$6 >= 20' > high_cadd_variants.tsv
```

## Verification

```bash
# Check tabix index works
tabix whole_genome_SNVs.tsv.gz 1:100000-100100

# Verify file integrity
md5sum whole_genome_SNVs.tsv.gz
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major version | Every 1-2 years |
| Annotation updates | As integrated data updates |

## Integration Notes

- VEP Plugin available: `--plugin CADD,whole_genome_SNVs.tsv.gz`
- Integrated into dbNSFP
- Can be used with ANNOVAR, SnpSift
