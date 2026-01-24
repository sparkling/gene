---
id: download-dbnsfp
title: "dbNSFP Download Instructions"
type: download
parent: _index.md
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
