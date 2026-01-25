---
id: download-dbvar
title: "dbVar Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# dbVar Download Instructions

## Quick Start

```bash
# Download non-redundant deletions (GRCh38)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.nr_deletions.vcf.gz
```

## Prerequisites

- wget, curl, or FTP client
- ~20 GB disk space for comprehensive data
- tabix for VCF indexing
- bcftools (recommended)

## Download Methods

### Primary: NCBI FTP

```bash
# Non-redundant variant sets (recommended for most uses)
# Deletions
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.nr_deletions.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.nr_deletions.vcf.gz.tbi

# Duplications
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.nr_duplications.vcf.gz

# Insertions
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.nr_insertions.vcf.gz

# All SV types
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_region.vcf.gz
```

### Alternative: Study-Specific Downloads

```bash
# List available studies
curl ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/

# Download specific study (e.g., gnomAD SVs - nstd166)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/nstd166/vcf/nstd166.GRCh38.variant_call.vcf.gz
```

### Alternative: BED Format

```bash
# BED files (simpler format)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/bed/GRCh38.nr_deletions.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/bed/GRCh38.nr_duplications.bed.gz
```

### Alternative: Web Interface

```bash
# Search and download via web
# https://www.ncbi.nlm.nih.gov/dbvar/

# Advanced search allows filtering by:
# - Variant type
# - Size range
# - Chromosome/region
# - Clinical significance
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| GRCh38.nr_deletions.vcf.gz | ~500 MB | Non-redundant deletions |
| GRCh38.nr_duplications.vcf.gz | ~200 MB | Non-redundant duplications |
| GRCh38.nr_insertions.vcf.gz | ~300 MB | Non-redundant insertions |
| GRCh38.variant_region.vcf.gz | ~2 GB | All variant regions |
| nstd166.*.vcf.gz | ~1 GB | gnomAD SV study |

## Post-Download Processing

```bash
# Index VCF files
tabix -p vcf GRCh38.nr_deletions.vcf.gz

# Query specific region
tabix GRCh38.nr_deletions.vcf.gz 1:1000000-2000000

# Filter by size (>10kb deletions)
bcftools view -i 'abs(SVLEN)>10000' GRCh38.nr_deletions.vcf.gz > large_dels.vcf

# Extract specific SV type
bcftools view -i 'SVTYPE="DEL"' GRCh38.variant_region.vcf.gz > deletions_only.vcf

# Convert to BED
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' GRCh38.nr_deletions.vcf.gz > deletions.bed
```

## Verification

```bash
# Count variants
bcftools view -H GRCh38.nr_deletions.vcf.gz | wc -l

# Check SV types present
bcftools query -f '%INFO/SVTYPE\n' GRCh38.variant_region.vcf.gz | sort | uniq -c
```

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | Continuous (2026-01) |
| Release Date | Daily updates |
| Total Size | ~3 GB (all VCFs) |
| SV Records | 7M+ structural variants |
| Studies | 500+ studies |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| GRCh38.nr_deletions.vcf.gz | ~500 MB | 2M | Non-redundant deletions |
| GRCh38.nr_duplications.vcf.gz | ~200 MB | 800K | Non-redundant duplications |
| GRCh38.nr_insertions.vcf.gz | ~300 MB | 1M | Non-redundant insertions |
| GRCh38.variant_region.vcf.gz | ~2 GB | 7M | All variant regions |

### Major Studies

| Study ID | Description | Variants |
|----------|-------------|----------|
| nstd166 | gnomAD-SV | 400K |
| nstd102 | 1000G Phase 3 SVs | 68K |
| nstd186 | TOPMed SVs | 200K |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://www.ncbi.nlm.nih.gov/dbvar/ |
| Rate Limit | 3 req/sec (10 with API key) |
| Auth Required | No |
| Response Format | XML, VCF, BED |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Data updates | Continuous |
| FTP refresh | Daily |
| Major studies | As submitted |

## Integration Notes

- Synchronized with European Variation Archive (EVA/DGVa)
- ClinVar SVs are forwarded to dbVar
- gnomAD-SV (nstd166) provides population frequencies
- Compatible with AnnotSV for clinical annotation
- Use non-redundant files for annotation to avoid double-counting
