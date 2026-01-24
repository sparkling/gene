---
id: download-cpic
title: "CPIC Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# CPIC Download Instructions

## Quick Start

```bash
# Access CPIC data via PharmGKB API
curl "https://api.pharmgkb.org/v1/data/clinicalAnnotation?view=max" -H "Accept: application/json"
```

## Prerequisites

- curl or web browser
- PharmGKB account (optional, for API)
- Minimal disk space (guidelines are documents)

## Download Methods

### Primary: CPIC Website

```bash
# Navigate to CPIC website for guidelines
# https://cpicpgx.org/guidelines/

# Each guideline page includes:
# - Full text PDF
# - Supplemental Excel files
# - Allele functionality tables
# - Diplotype-phenotype tables
```

### Alternative: PharmGKB API

```bash
# Get all CPIC guidelines
curl "https://api.pharmgkb.org/v1/data/guideline?source=cpic" \
  -H "Accept: application/json" > cpic_guidelines.json

# Get specific gene-drug annotation
curl "https://api.pharmgkb.org/v1/data/clinicalAnnotation?location.genes.symbol=CYP2D6&relatedChemicals.name=codeine" \
  -H "Accept: application/json"

# Get allele functionality
curl "https://api.pharmgkb.org/v1/data/gene/PA128/alleles" \
  -H "Accept: application/json"
```

### Alternative: PharmGKB Downloads

```bash
# Download curated tables from PharmGKB
# https://www.pharmgkb.org/downloads

# Clinical annotations
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip

# Allele definitions
wget https://api.pharmgkb.org/v1/download/file/data/alleleDefinitions.zip

# Allele functionality
wget https://api.pharmgkb.org/v1/download/file/data/alleleFunctionality.zip
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| Guideline PDFs | 1-5 MB each | Full text documents |
| Supplemental Excel | 100-500 KB | Data tables |
| clinicalAnnotations.zip | ~5 MB | All annotations |
| alleleDefinitions.zip | ~2 MB | Star allele definitions |

## Post-Download Processing

```bash
# Parse PharmGKB JSON
cat cpic_guidelines.json | jq '.data[] | {gene: .gene.symbol, drug: .relatedChemicals[0].name}'

# Extract specific gene
cat cpic_guidelines.json | jq '.data[] | select(.gene.symbol == "CYP2D6")'

# Convert to TSV (using jq)
cat cpic_guidelines.json | jq -r '.data[] | [.gene.symbol, .relatedChemicals[0].name, .recommendation.text] | @tsv'
```

## Verification

```bash
# Count guidelines
cat cpic_guidelines.json | jq '.data | length'

# Check specific guideline exists
cat cpic_guidelines.json | jq '.data[] | select(.relatedChemicals[0].name == "codeine")'
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New guidelines | As published |
| Guideline updates | Continuous |
| PharmGKB sync | Real-time |

## Integration Notes

- CPIC content is free and open access
- Citation required for each guideline
- PharmGKB provides structured data access
- Clinical decision support implementations should use official tables
- Some EHR vendors have pre-built CPIC integrations
