---
id: download-cancer-gene-census
title: "Cancer Gene Census Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Cancer Gene Census Download Instructions

## Quick Start

```bash
# Register at COSMIC, then download via authenticated request
curl -H "Authorization: Basic YOUR_TOKEN" \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/cancer_gene_census.csv"
```

## Prerequisites

- COSMIC account (free registration)
- curl or web browser
- Minimal disk space (<10 MB)

## Download Methods

### Primary: COSMIC Website (Registered)

```bash
# 1. Register at https://cancer.sanger.ac.uk/cosmic/register

# 2. Login and navigate to Downloads
# https://cancer.sanger.ac.uk/cosmic/download

# 3. Select "Cancer Gene Census" from available files

# 4. Download directly or via command line with authentication
# Get your token from account settings

# Download CGC
curl -H "Authorization: Basic $(echo -n 'email:password' | base64)" \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/cancer_gene_census.csv" \
  -o cancer_gene_census.csv
```

### Alternative: Web Download

```bash
# Navigate to COSMIC downloads
# https://cancer.sanger.ac.uk/cosmic/download

# Select:
# - Cancer Gene Census (CGC)
# - Click download button
# - File downloads as CSV
```

### Alternative: API Query

```bash
# COSMIC API for gene lookup
# Requires authentication

# Get gene information
curl -H "Authorization: Basic YOUR_TOKEN" \
  "https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=TP53"
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| cancer_gene_census.csv | ~500 KB | Full CGC list |
| cancer_gene_census.tsv | ~500 KB | Tab-separated version |
| hallmarks_of_cancer.csv | ~100 KB | Hallmark annotations |

## Post-Download Processing

```bash
# View file structure
head -1 cancer_gene_census.csv

# Count genes
wc -l cancer_gene_census.csv

# Extract Tier 1 genes
awk -F',' '$5 == 1' cancer_gene_census.csv > tier1_genes.csv

# Extract oncogenes
grep -i "oncogene" cancer_gene_census.csv > oncogenes.csv

# Extract tumor suppressors
grep -i "TSG" cancer_gene_census.csv > tsgs.csv

# Get genes for specific cancer type
grep -i "breast" cancer_gene_census.csv > breast_cancer_genes.csv
```

## Verification

```bash
# Count total genes
tail -n +2 cancer_gene_census.csv | wc -l
# Expected: ~736

# Count by tier
cut -d',' -f5 cancer_gene_census.csv | sort | uniq -c

# Verify specific gene present
grep "TP53" cancer_gene_census.csv
```

## Dataset Versions

### Current Release: COSMIC v103

| Property | Value |
|----------|-------|
| Version | v103 (CGC) |
| Release Date | 2025-11-18 |
| Total Size | ~500 KB |
| Tier 1 Genes | 350+ |
| Tier 2 Genes | 400+ |
| Total CGC Genes | 736 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| cancer_gene_census.csv | ~500 KB | 736 | Full CGC list |
| hallmarks_of_cancer.csv | ~100 KB | 736 | Hallmark annotations |
| cancer_gene_census.tsv | ~500 KB | 736 | Tab-separated version |

### Gene Classification

| Category | Count | Description |
|----------|-------|-------------|
| Tier 1 | 350+ | Strong evidence, well-characterized |
| Tier 2 | 400+ | Emerging evidence |
| Oncogenes | 350+ | Gain-of-function drivers |
| TSGs | 300+ | Tumor suppressors |
| Fusion genes | 150+ | Recurrent fusions |

### Hallmarks Coverage

| Hallmark | Genes | Description |
|----------|-------|-------------|
| Proliferation | 200+ | Cell cycle, growth signals |
| Apoptosis evasion | 150+ | Death resistance |
| Invasion/metastasis | 100+ | Tissue spread |
| Angiogenesis | 50+ | Blood vessel formation |
| Immune evasion | 80+ | Immune escape |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://cancer.sanger.ac.uk/cosmic/api/ |
| Rate Limit | 10 req/sec |
| Auth Required | Yes (COSMIC account) |
| Response Format | CSV, JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| CGC updates | Quarterly |
| COSMIC releases | Quarterly (v103, v104, etc.) |
| New genes | As evidence accumulated |

## Integration Notes

- CGC is part of COSMIC database
- Free for academic use with registration
- Commercial use requires COSMIC license
- Gene symbols match HGNC nomenclature
- Cross-reference with COSMIC for mutation details
- Can be used for cancer panel design and variant interpretation
