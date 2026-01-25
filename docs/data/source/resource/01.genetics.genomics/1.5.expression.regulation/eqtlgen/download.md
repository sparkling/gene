---
id: download-eqtlgen
title: "eQTLGen Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# eQTLGen Download Instructions

## Quick Start

```bash
# Download significant cis-eQTLs
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonijValsRemoved.txt.gz
```

## Prerequisites

- wget or curl
- ~10 GB disk space for full data
- gzip for decompression
- R or Python for analysis (recommended)

## Download Methods

### Primary: eQTLGen Downloads

```bash
# Significant cis-eQTLs (FDR < 0.05)
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonijValsRemoved.txt.gz

# Full cis-eQTL results (all tested)
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonijValsRemoved.txt.gz

# Trans-eQTLs (significant)
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonigValsRemoved.txt.gz

# Full trans-eQTL results
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonijValsRemoved.txt.gz
```

### Alternative: Web Downloads

```bash
# Navigate to eQTLGen download page
# https://www.eqtlgen.org/cis-eqtls.html
# https://www.eqtlgen.org/trans-eqtls.html

# Select specific datasets:
# - Full summary statistics
# - Significant results only
# - SMR-ready formats
```

### Alternative: API Lookup

```bash
# Query specific SNP-gene pair
curl "https://www.eqtlgen.org/api/v1/cis-eqtl?snp=rs12345&gene=ENSG00000012048"
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| cis-eQTLsFDR0.05*.txt.gz | ~200 MB | Significant cis-eQTLs |
| cis-eQTLsFDR*.txt.gz | ~5 GB | All cis-eQTL tests |
| trans-eQTLsFDR0.05*.txt.gz | ~50 MB | Significant trans-eQTLs |
| trans-eQTLsFDR*.txt.gz | ~2 GB | All trans-eQTL tests |

## Post-Download Processing

```bash
# Decompress
gunzip 2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonijValsRemoved.txt.gz

# View header
head -1 2019-12-11-cis-eQTLsFDR0.05*.txt

# Filter for specific gene
grep "BRCA1" 2019-12-11-cis-eQTLsFDR0.05*.txt > brca1_eqtls.txt

# Filter for genome-wide significant
awk '$13 < 5e-8' 2019-12-11-cis-eQTLsFDR0.05*.txt > gwas_significant_eqtls.txt
```

## Verification

```bash
# Count significant cis-eQTLs
zcat 2019-12-11-cis-eQTLsFDR0.05*.txt.gz | wc -l
# Expected: ~16,987 significant cis-eQTLs

# Check column headers
zcat 2019-12-11-cis-eQTLsFDR0.05*.txt.gz | head -1
```

## Dataset Versions

### Current Release: Phase I

| Property | Value |
|----------|-------|
| Version | Phase I (2019-12-11) |
| Release Date | 2019-12-11 |
| Total Size | ~7 GB |
| Samples | 31,684 individuals |
| Cohorts | 37 cohorts |
| Reference | GRCh37/hg19 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| cis-eQTLsFDR0.05*.txt.gz | ~200 MB | 16,987 | Significant cis-eQTLs |
| cis-eQTLsFDR*.txt.gz | ~5 GB | 10M+ | All cis-eQTL tests |
| trans-eQTLsFDR0.05*.txt.gz | ~50 MB | 3,853 | Significant trans-eQTLs |
| trans-eQTLsFDR*.txt.gz | ~2 GB | 100M+ | All trans-eQTL tests |

### eQTL Statistics

| Category | Count | Description |
|----------|-------|-------------|
| Tested genes | 19,250 | Protein-coding genes |
| Significant cis-eGenes | 16,987 | FDR < 0.05 |
| Significant trans-eGenes | 3,853 | FDR < 0.05 |
| Tested SNPs | 10M+ | Common variants |

### Tissue/Cell Type

| Type | Description |
|------|-------------|
| Blood | Whole blood eQTLs |
| PBMC | Peripheral blood cells |
| Monocytes | Subset available |
| CD4+ T cells | Subset available |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://www.eqtlgen.org/api/v1/ |
| Rate Limit | 10 req/sec |
| Auth Required | No |
| Response Format | JSON, TSV |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Phase I | Released 2019 |
| Phase II | In progress (single-cell) |
| Updates | As new analyses completed |

## Integration Notes

- Coordinates are hg19/GRCh37 (liftover may be needed for GRCh38)
- Compatible with SMR (Summary-based Mendelian Randomization)
- Can be used for colocalization analyses (coloc, eCAVIAR)
- Useful for GWAS signal interpretation and fine-mapping
