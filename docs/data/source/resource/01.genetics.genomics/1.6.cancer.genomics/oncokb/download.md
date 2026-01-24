---
id: download-oncokb
title: "OncoKB Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# OncoKB Download Instructions

## Quick Start

```bash
# Annotate a mutation via API (requires token)
curl -X GET "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=BRAF&alteration=V600E" \
  -H "Authorization: Bearer YOUR_API_TOKEN"
```

## Prerequisites

- OncoKB account (free for academic)
- API token from account settings
- curl or HTTP client
- Minimal disk space for API responses

## Download Methods

### Primary: OncoKB API

```bash
# Get API token from https://www.oncokb.org/account/settings

# Annotate single mutation
curl -X GET "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=BRAF&alteration=V600E" \
  -H "Authorization: Bearer YOUR_API_TOKEN" > braf_v600e.json

# Annotate by genomic coordinates
curl -X GET "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange?genomicLocation=7,140453136,140453136,T,A&referenceGenome=GRCh38" \
  -H "Authorization: Bearer YOUR_API_TOKEN"

# Annotate CNA
curl -X GET "https://www.oncokb.org/api/v1/annotate/copyNumberAlterations?hugoSymbol=ERBB2&copyNameAlterationType=AMPLIFICATION" \
  -H "Authorization: Bearer YOUR_API_TOKEN"

# Annotate fusion
curl -X GET "https://www.oncokb.org/api/v1/annotate/structuralVariants?hugoSymbolA=BCR&hugoSymbolB=ABL1&structuralVariantType=FUSION" \
  -H "Authorization: Bearer YOUR_API_TOKEN"

# Get all gene info
curl -X GET "https://www.oncokb.org/api/v1/genes" \
  -H "Authorization: Bearer YOUR_API_TOKEN" > all_genes.json

# Get evidence levels
curl -X GET "https://www.oncokb.org/api/v1/levels" \
  -H "Authorization: Bearer YOUR_API_TOKEN"
```

### Alternative: Batch Annotation

```bash
# Batch annotate mutations (POST)
curl -X POST "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange" \
  -H "Authorization: Bearer YOUR_API_TOKEN" \
  -H "Content-Type: application/json" \
  -d '[
    {"hugoSymbol": "BRAF", "alteration": "V600E"},
    {"hugoSymbol": "KRAS", "alteration": "G12D"},
    {"hugoSymbol": "EGFR", "alteration": "L858R"}
  ]' > batch_annotations.json
```

### Alternative: Data Files (Registered Users)

```bash
# Login to OncoKB and navigate to Downloads
# https://www.oncokb.org/account/settings

# Available downloads (academic license):
# - All annotated variants
# - Therapeutic implications
# - Gene summaries

# Files delivered via download link after request approval
```

### Alternative: cBioPortal Integration

```bash
# OncoKB annotations are available in cBioPortal
# When querying cBioPortal, enable OncoKB column

# Via cBioPortal API
curl "https://www.cbioportal.org/api/molecular-profiles/brca_tcga_mutations/mutations?sampleListId=brca_tcga_all" \
  -H "Accept: application/json" > brca_with_oncokb.json
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| API responses | Varies | Per-query JSON |
| all_genes.json | ~1 MB | All annotated genes |
| data_files (registered) | ~10 MB | Complete annotations |

## Post-Download Processing

```bash
# Parse API response
cat braf_v600e.json | jq '.oncogenic, .treatments[0].level'

# Extract actionable variants
cat batch_annotations.json | jq '.[] | select(.treatments != null) | {gene: .query.hugoSymbol, alteration: .query.alteration, level: .treatments[0].level}'

# Filter for Level 1 evidence
cat batch_annotations.json | jq '.[] | select(.treatments[0].level == "1")'
```

## Verification

```bash
# Test API connection
curl -I -X GET "https://www.oncokb.org/api/v1/info" \
  -H "Authorization: Bearer YOUR_API_TOKEN"

# Verify BRAF V600E returns expected result
cat braf_v600e.json | jq '.oncogenic'
# Expected: "Oncogenic"
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Data updates | Weekly |
| Evidence review | Continuous |
| FDA label updates | As approved |

## Integration Notes

- Free for academic, license required for commercial
- API rate limits apply (check documentation)
- FDA-recognized tumor mutation database
- Integrated into cBioPortal, various commercial platforms
- VICC harmonization with CIViC and other knowledgebases
- Contact oncokb@oncokb.org for commercial licensing
