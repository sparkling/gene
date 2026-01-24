---
id: download-pharmvar
title: "PharmVar Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# PharmVar Download Instructions

## Quick Start

```bash
# Download CYP2D6 allele definitions
curl "https://www.pharmvar.org/api-service/alleles?gene=CYP2D6" \
  -H "Accept: application/json" > cyp2d6_alleles.json
```

## Prerequisites

- curl or web browser
- Minimal disk space
- JSON parser (jq recommended)

## Download Methods

### Primary: PharmVar API

```bash
# List all genes
curl "https://www.pharmvar.org/api-service/genes" \
  -H "Accept: application/json"

# Get alleles for specific gene
curl "https://www.pharmvar.org/api-service/alleles?gene=CYP2D6" \
  -H "Accept: application/json" > cyp2d6_alleles.json

curl "https://www.pharmvar.org/api-service/alleles?gene=CYP2C19" \
  -H "Accept: application/json" > cyp2c19_alleles.json

# Get specific allele details
curl "https://www.pharmvar.org/api-service/allele?name=CYP2D6*4" \
  -H "Accept: application/json"
```

### Alternative: Web Downloads

```bash
# Navigate to PharmVar download section
# https://www.pharmvar.org/download

# Available files per gene:
# - Allele definition tables (TSV)
# - Reference sequences (FASTA)
# - Functional status tables

# Example direct download
wget "https://www.pharmvar.org/get-download-file?name=CYP2D6&refSeq=NG_008376.4&alleleDef=true&type=json"
```

### Alternative: Gene-Specific Pages

```bash
# Each gene has dedicated page with downloads
# https://www.pharmvar.org/gene/CYP2D6

# Click "Download" tab for:
# - TSV allele definitions
# - FASTA sequences
# - VCF files (where available)
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| Gene_alleles.json | 50-500 KB | API response per gene |
| Gene_allele_definitions.tsv | 20-200 KB | Tabular definitions |
| Gene_reference.fasta | 10-100 KB | Reference sequence |
| Gene.vcf | Varies | Variant positions |

## Post-Download Processing

```bash
# Parse JSON for allele names and functions
cat cyp2d6_alleles.json | jq '.[] | {name: .alleleName, function: .function}'

# Extract defining variants
cat cyp2d6_alleles.json | jq '.[] | select(.alleleName == "*4") | .definingVariants'

# Convert to simple table
cat cyp2d6_alleles.json | jq -r '.[] | [.alleleName, .function // "unknown", .activityValue // "N/A"] | @tsv' > cyp2d6_table.tsv
```

## Verification

```bash
# Count alleles per gene
cat cyp2d6_alleles.json | jq 'length'

# Verify specific allele exists
cat cyp2d6_alleles.json | jq '.[] | select(.alleleName == "*4")'
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New alleles | As discovered/submitted |
| Database updates | Continuous |
| Major releases | Quarterly |

## Integration Notes

- PharmVar is the authoritative source for star allele nomenclature
- Used by CPIC, DPWG, and clinical laboratories
- Star allele callers (Stargazer, Aldy, Cyrius) reference PharmVar
- Free for academic and commercial use with citation
- Allele assignments require matching all defining variants
