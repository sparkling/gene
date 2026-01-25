---
id: download-dpwg
title: "DPWG Download Instructions"
type: download
parent: README.md
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

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | 2025-Q4 |
| Release Date | 2025-12 |
| Total Size | ~5 MB |
| Guidelines | 100+ gene-drug pairs |
| Languages | Dutch, English (PharmGKB) |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| dpwg_guidelines.json (PharmGKB) | ~2 MB | 100+ | All DPWG guidelines |
| clinicalAnnotations.zip | ~5 MB | 1000+ | Including DPWG annotations |
| G-Standaard (Dutch) | varies | 100+ | Dutch pharmacy system |

### Gene-Drug Coverage

| Category | Examples | Count |
|----------|----------|-------|
| CYP2D6 drugs | codeine, tramadol, tamoxifen | 25+ |
| CYP2C19 drugs | clopidogrel, omeprazole | 15+ |
| CYP2C9 drugs | warfarin, phenytoin | 10+ |
| Other genes | DPYD, TPMT, VKORC1 | 50+ |

### Comparison with CPIC

| Aspect | DPWG | CPIC |
|--------|------|------|
| Origin | Netherlands | USA |
| Focus | European practice | Global |
| Integration | G-Standaard | EMR/EHR |
| Harmonization | Ongoing | Ongoing |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://api.pharmgkb.org/v1/data/guideline |
| Rate Limit | 10 req/sec |
| Auth Required | No (PharmGKB account recommended) |
| Response Format | JSON |

---

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
