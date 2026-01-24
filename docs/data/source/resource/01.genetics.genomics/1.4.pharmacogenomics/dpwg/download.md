---
id: download-dpwg
title: "DPWG Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# DPWG Download Instructions

## Quick Start

```bash
# Access DPWG guidelines via PharmGKB API
curl "https://api.pharmgkb.org/v1/data/guideline?source=dpwg" \
  -H "Accept: application/json" > dpwg_guidelines.json
```

## Prerequisites

- curl or web browser
- PharmGKB account (optional)
- Minimal disk space

## Download Methods

### Primary: PharmGKB (International Access)

```bash
# Get all DPWG guidelines
curl "https://api.pharmgkb.org/v1/data/guideline?source=dpwg" \
  -H "Accept: application/json" > dpwg_guidelines.json

# Get specific gene annotations
curl "https://api.pharmgkb.org/v1/data/clinicalAnnotation?source=DPWG&location.genes.symbol=CYP2C19" \
  -H "Accept: application/json"

# Download structured data files
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip
```

### Alternative: KNMP (Dutch Source)

```bash
# KNMP provides guidelines in Dutch via G-Standaard
# https://www.knmp.nl/

# Requires KNMP membership for full access
# Academic access may be available upon request
```

### Alternative: PharmGKB Web Interface

```bash
# Browse DPWG guidelines online
# https://www.pharmgkb.org/page/dpwg

# Individual guideline pages include:
# - Recommendation text
# - Gene-drug interaction details
# - Comparison with CPIC guidelines
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| dpwg_guidelines.json | ~2 MB | All DPWG via PharmGKB |
| clinicalAnnotations.zip | ~5 MB | All sources including DPWG |
| G-Standaard export | Varies | Dutch national (restricted) |

## Post-Download Processing

```bash
# Parse PharmGKB JSON
cat dpwg_guidelines.json | jq '.data[] | {gene: .gene.symbol, drug: .relatedChemicals[0].name}'

# Filter urgent recommendations
cat dpwg_guidelines.json | jq '.data[] | select(.annotation.urgency == "Urgent")'

# Compare DPWG with CPIC for same gene-drug
# Download both and join on gene+drug fields
```

## Verification

```bash
# Count DPWG guidelines
cat dpwg_guidelines.json | jq '.data | length'

# Verify specific guideline
cat dpwg_guidelines.json | jq '.data[] | select(.relatedChemicals[0].name == "clopidogrel")'
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| DPWG updates | Quarterly to G-Standaard |
| PharmGKB sync | Monthly |
| New guidelines | As published |

## Integration Notes

- PharmGKB provides English translations
- CPIC-DPWG harmonization ongoing for key gene-drug pairs
- G-Standaard integration for Dutch EHR systems
- CC BY-SA 4.0 license for PharmGKB-accessed content
- Contact KNMP for commercial use of G-Standaard
