---
id: download-dbnsfp
title: "dbNSFP Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# dbNSFP Download Instructions

## Quick Start

```bash
# Download from Amazon S3
wget https://dbnsfp.s3.amazonaws.com/dbNSFP5.3.1a.zip
```

## Prerequisites

- wget or curl
- ~50 GB disk space
- unzip utility
- Java (for search utility, optional)

## Download Methods

### Primary: Amazon S3

```bash
# Full database (GRCh38 + GRCh37)
wget https://dbnsfp.s3.amazonaws.com/dbNSFP5.3.1a.zip

# Gene annotations only (smaller)
wget https://dbnsfp.s3.amazonaws.com/dbNSFP5.3_gene.complete.gz
```

### Alternative: Academic Institutions

```bash
# Some mirrors available; check website
# https://www.dbnsfp.org/
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| dbNSFP5.3.1a.zip | ~40 GB | Full database |
| dbNSFP5.3_gene.complete.gz | ~100 MB | Gene-level only |
| dbNSFP5.3_variant.chr*.gz | Varies | Per-chromosome |

## Post-Download Processing

```bash
# Unzip
unzip dbNSFP5.3.1a.zip -d dbNSFP5.3

# Navigate to directory
cd dbNSFP5.3

# Concatenate chromosomes (optional, creates single file)
zcat dbNSFP5.3_variant.chr1.gz > dbNSFP5.3_variant.txt
for i in {2..22} X Y M; do
  zcat dbNSFP5.3_variant.chr${i}.gz | tail -n +2 >> dbNSFP5.3_variant.txt
done

# Build search index (requires Java)
java search_dbNSFP53 -b
```

## Verification

```bash
# Check file count
ls -la dbNSFP5.3/dbNSFP5.3_variant.chr*.gz | wc -l
# Expected: 24 files (1-22, X, Y)

# Verify column count
zcat dbNSFP5.3_variant.chr1.gz | head -1 | tr '\t' '\n' | wc -l
# Expected: 600+ columns
```

## Dataset Versions

### Current Release: v5.3

| Property | Value |
|----------|-------|
| Version | 5.3.1a |
| Release Date | 2025-10 |
| Total Size | ~40 GB (zipped) |
| nsSNVs | 83M non-synonymous |
| ssSNVs | 2.4M splice-site |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| dbNSFP5.3.1a.zip | ~40 GB | 85M | Full database |
| dbNSFP5.3_gene.complete.gz | ~100 MB | 19K | Gene-level annotations |
| dbNSFP5.3_variant.chr*.gz | varies | per chr | Per-chromosome files |

### Prediction Algorithms (36 total)

| Category | Algorithms |
|----------|------------|
| Conservation | SIFT, SIFT4G, PROVEAN, phyloP, phastCons |
| Functional | PolyPhen2, MutationTaster, FATHMM-XF |
| Ensemble | CADD, REVEL, MetaLR, MetaRNN, MetaSVM |
| Deep Learning | ESM1b, AlphaMissense, MutFormer |

### Previous Versions

| Version | Release | Size | Status |
|---------|---------|------|--------|
| v5.2 | 2025-06 | ~38 GB | Archived |
| v5.1 | 2025-02 | ~36 GB | Archived |
| v5.0 | 2024-08 | ~35 GB | Archived |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://dbnsfp.s3.amazonaws.com/ |
| Rate Limit | N/A (bulk download) |
| Auth Required | No |
| Response Format | TSV (zipped) |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major version | Every 6-12 months |
| Minor updates | As new predictors added |

## Integration Notes

- VEP Plugin: dbNSFP plugin available
- SnpSift: `SnpSift dbnsfp` command
- ANNOVAR: Use dbnsfp protocol
- Some commercial-use restrictions on individual predictors
