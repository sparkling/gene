---
id: download-cpic
title: "CPIC Download Instructions"
type: download
parent: README.md
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

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | 2025-Q4 |
| Release Date | Continuous updates |
| Total Size | ~10 MB |
| Guidelines | 27 published |
| Gene-Drug Pairs | 135+ |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| Guideline PDFs | 1-5 MB each | 27 | Full text documents |
| Supplemental Excel | 100-500 KB | 27 | Data tables |
| clinicalAnnotations.zip | ~5 MB | 750+ | PharmGKB format |
| alleleDefinitions.zip | ~2 MB | 25 genes | Star allele definitions |
| alleleFunctionality.zip | ~1 MB | 25 genes | Function assignments |

### Published Guidelines

| Category | Genes | Examples |
|----------|-------|----------|
| Analgesics | CYP2D6 | codeine, tramadol |
| Anticoagulants | CYP2C9, VKORC1 | warfarin |
| Antidepressants | CYP2D6, CYP2C19 | SSRIs, TCAs |
| Oncology | DPYD, TPMT, UGT1A1 | fluoropyrimidines |
| Cardiology | CYP2C19, SLCO1B1 | clopidogrel, statins |

### Guideline Status

| Status | Count | Description |
|--------|-------|-------------|
| Published | 27 | Full CPIC guideline |
| In Progress | 5+ | Under development |
| Update Pending | 3 | Revision in progress |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://api.pharmgkb.org/v1/data/guideline |
| Rate Limit | 10 req/sec |
| Auth Required | No |
| Response Format | JSON |

---

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
