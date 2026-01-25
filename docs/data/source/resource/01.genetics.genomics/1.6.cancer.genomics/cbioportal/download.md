---
id: download-cbioportal
title: "cBioPortal Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# cBioPortal Download Instructions

## Quick Start

```bash
# Download mutations for a study via API
curl "https://www.cbioportal.org/api/molecular-profiles/brca_tcga_mutations/mutations?sampleListId=brca_tcga_all" \
  -H "Accept: application/json" > brca_tcga_mutations.json
```

## Prerequisites

- curl or web browser
- Varies by study (100 MB to 10+ GB)
- JSON parser (jq recommended)
- R or Python for analysis

## Download Methods

### Primary: cBioPortal API

```bash
# List all studies
curl "https://www.cbioportal.org/api/studies" \
  -H "Accept: application/json" > all_studies.json

# Get study details
curl "https://www.cbioportal.org/api/studies/brca_tcga" \
  -H "Accept: application/json"

# Get samples in study
curl "https://www.cbioportal.org/api/studies/brca_tcga/samples" \
  -H "Accept: application/json" > brca_samples.json

# Get mutations
curl "https://www.cbioportal.org/api/molecular-profiles/brca_tcga_mutations/mutations?sampleListId=brca_tcga_all" \
  -H "Accept: application/json" > brca_mutations.json

# Get clinical data
curl "https://www.cbioportal.org/api/studies/brca_tcga/clinical-data" \
  -H "Accept: application/json" > brca_clinical.json
```

### Alternative: Web Download

```bash
# Navigate to cBioPortal
# https://www.cbioportal.org/

# Select study → Download → Download tab
# Available formats: TSV, MAF, Seg files

# Direct study download links
# https://www.cbioportal.org/study/summary?id=brca_tcga
```

### Alternative: Datahub (Bulk Download)

```bash
# Clone cBioPortal datahub for bulk data
git clone https://github.com/cBioPortal/datahub.git

# Navigate to specific study
cd datahub/public/
ls brca_tcga/

# Files include:
# - data_mutations.txt (MAF format)
# - data_CNA.txt (copy number)
# - data_clinical_patient.txt
# - data_clinical_sample.txt
```

### Alternative: R Package (cgdsr)

```r
# Install package
install.packages("cgdsr")
library(cgdsr)

# Connect
mycgds <- CGDS("https://www.cbioportal.org/")

# Get studies
studies <- getCancerStudies(mycgds)

# Get mutation data
mutations <- getMutationData(mycgds, "brca_tcga_all", "brca_tcga_mutations", "TP53")
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| data_mutations.txt | 10-500 MB | MAF-format mutations |
| data_CNA.txt | 1-50 MB | Copy number data |
| data_clinical_patient.txt | 1-10 MB | Patient clinical data |
| data_clinical_sample.txt | 1-10 MB | Sample clinical data |
| data_mrna_seq_v2_rsem.txt | 10-100 MB | Expression data |

## Post-Download Processing

```bash
# Parse API JSON response
cat brca_mutations.json | jq '.[].hugoGeneSymbol' | sort | uniq -c | sort -rn | head

# Extract specific gene mutations
cat brca_mutations.json | jq '.[] | select(.hugoGeneSymbol == "TP53")'

# Convert to TSV
cat brca_mutations.json | jq -r '.[] | [.hugoGeneSymbol, .proteinChange, .mutationType, .sampleId] | @tsv' > mutations.tsv

# Filter for specific mutation types
grep "Missense" data_mutations.txt > missense_only.txt
```

## Verification

```bash
# Count mutations per gene (API data)
cat brca_mutations.json | jq '.[].hugoGeneSymbol' | sort | uniq -c | sort -rn | head -20

# Count samples
cat brca_mutations.json | jq '.[].sampleId' | sort -u | wc -l
```

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | Continuous (2026-01) |
| Release Date | Daily updates |
| Total Size | ~500 GB (all studies) |
| Studies | 400+ |
| Samples | 200,000+ |
| Genes | 20,000+ |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| data_mutations.txt (per study) | 10-500 MB | varies | MAF-format mutations |
| data_CNA.txt (per study) | 1-50 MB | varies | Copy number data |
| data_clinical_patient.txt | 1-10 MB | varies | Patient clinical data |
| data_clinical_sample.txt | 1-10 MB | varies | Sample clinical data |
| data_mrna_seq_v2_rsem.txt | 10-100 MB | varies | Expression data |

### Major Studies

| Study | Samples | Cancer Type |
|-------|---------|-------------|
| TCGA Pan-Cancer | 11,000+ | 33 cancer types |
| METABRIC | 2,500+ | Breast cancer |
| MSK-IMPACT | 50,000+ | Various cancers |
| GENIE | 100,000+ | Multi-institution |
| PCAWG | 2,800+ | Whole genomes |

### Data Types Available

| Type | Coverage | Description |
|------|----------|-------------|
| Mutations | 98% studies | SNVs, indels |
| CNA | 90% studies | Copy number |
| Expression | 60% studies | mRNA levels |
| Methylation | 30% studies | DNA methylation |
| Fusions | 50% studies | Structural variants |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://www.cbioportal.org/api/ |
| Rate Limit | 100 req/min |
| Auth Required | No |
| Response Format | JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New studies | Continuous |
| Data updates | As published |
| API version | Stable |

## Integration Notes

- Public API, no authentication required
- ODbL license allows commercial use with attribution
- Rate limiting applies to API (use delays between requests)
- Python client: `pip install cbioportalR`
- Compatible with maftools for mutation analysis
- OncoKB annotations available in portal
