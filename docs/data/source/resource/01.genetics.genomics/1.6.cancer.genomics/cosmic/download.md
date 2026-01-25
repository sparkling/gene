---
id: download-cosmic
title: "COSMIC Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# COSMIC Download Instructions

## Quick Start

```bash
# After registration, download with authentication
curl -H "Authorization: Basic $(echo -n 'email:password' | base64)" \
  -o CosmicMutantExport.tsv.gz \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/CosmicMutantExport.tsv.gz"
```

## Prerequisites

- **COSMIC account** (free registration required)
- **curl** or **wget** for downloads
- **gunzip** for decompression
- 50-200GB disk space depending on files downloaded

## Registration Process

### Step 1: Create Account

1. Navigate to https://cancer.sanger.ac.uk/cosmic/register
2. Provide email and institution details
3. Accept terms of use (free for academic/non-commercial)
4. Verify email address

### Step 2: Access Downloads

1. Log in at https://cancer.sanger.ac.uk/cosmic
2. Navigate to https://cancer.sanger.ac.uk/cosmic/download
3. Select genome version (GRCh38 recommended)
4. Generate authentication token for command-line access

## Download Methods

### Method 1: Web Interface

1. Log in to COSMIC
2. Go to Data Downloads page
3. Select file category (Mutations, Genes, etc.)
4. Click download link

### Method 2: Command Line (with Authentication)

```bash
# Set credentials
COSMIC_EMAIL="your@email.com"
COSMIC_PASSWORD="your_password"
AUTH=$(echo -n "${COSMIC_EMAIL}:${COSMIC_PASSWORD}" | base64)

# Download complete mutation export
curl -H "Authorization: Basic ${AUTH}" \
  -o CosmicMutantExport.tsv.gz \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/CosmicMutantExport.tsv.gz"

# Download coding mutations
curl -H "Authorization: Basic ${AUTH}" \
  -o CosmicCodingMuts.vcf.gz \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/VCF/CosmicCodingMuts.vcf.gz"

# Download gene census
curl -H "Authorization: Basic ${AUTH}" \
  -o cancer_gene_census.csv \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/cancer_gene_census.csv"
```

### Method 3: SFTP Access (Large Transfers)

```bash
# Connect via SFTP
sftp cosmic@sftp-cancer.sanger.ac.uk

# Navigate and download
cd /cosmic/grch38/cosmic/v99/
get CosmicMutantExport.tsv.gz

# Batch download
mget *.tsv.gz
```

### Method 4: Specific Data Downloads

```bash
# Gene census (cancer genes)
curl -H "Authorization: Basic ${AUTH}" \
  -o cancer_gene_census.csv \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/cancer_gene_census.csv"

# Fusions
curl -H "Authorization: Basic ${AUTH}" \
  -o CosmicFusionExport.tsv.gz \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/CosmicFusionExport.tsv.gz"

# CNV data
curl -H "Authorization: Basic ${AUTH}" \
  -o CosmicCompleteCNA.tsv.gz \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/CosmicCompleteCNA.tsv.gz"

# Resistance mutations
curl -H "Authorization: Basic ${AUTH}" \
  -o CosmicResistanceMutations.tsv.gz \
  "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v99/CosmicResistanceMutations.tsv.gz"
```

## File Inventory

### Core Mutation Data

| File | Size | Description |
|------|------|-------------|
| CosmicMutantExport.tsv.gz | ~5 GB | Complete mutation export |
| CosmicCodingMuts.vcf.gz | ~2 GB | Coding mutations (VCF) |
| CosmicNonCodingVariants.vcf.gz | ~3 GB | Non-coding variants |
| CosmicMutantExportCensus.tsv.gz | ~500 MB | Census genes mutations only |

### Gene and Target Data

| File | Size | Description |
|------|------|-------------|
| cancer_gene_census.csv | ~500 KB | Cancer Gene Census |
| CosmicHGNC.tsv.gz | ~1 MB | HGNC gene mapping |
| All_COSMIC_Genes.fasta.gz | ~50 MB | Gene sequences |

### Structural Variants

| File | Size | Description |
|------|------|-------------|
| CosmicFusionExport.tsv.gz | ~100 MB | Gene fusions |
| CosmicCompleteCNA.tsv.gz | ~2 GB | Copy number alterations |
| CosmicBreakpointsExport.tsv.gz | ~200 MB | Translocation breakpoints |

### Clinical Data

| File | Size | Description |
|------|------|-------------|
| CosmicSample.tsv.gz | ~500 MB | Sample metadata |
| CosmicResistanceMutations.tsv.gz | ~10 MB | Drug resistance mutations |
| classification.csv | ~1 MB | Tumor classifications |

## Post-Download Processing

```bash
# Decompress files
gunzip CosmicMutantExport.tsv.gz

# Extract specific columns
cut -f1,2,3,17,18,19 CosmicMutantExport.tsv > mutations_basic.tsv

# Filter for specific gene
grep "BRAF" CosmicMutantExport.tsv > braf_mutations.tsv

# Convert VCF to BED
bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' CosmicCodingMuts.vcf.gz > cosmic.bed

# Index VCF
tabix -p vcf CosmicCodingMuts.vcf.gz
```

## Verification

```bash
# Check file integrity
gzip -t CosmicMutantExport.tsv.gz

# Count mutations
zcat CosmicMutantExport.tsv.gz | wc -l

# Preview file structure
zcat CosmicMutantExport.tsv.gz | head -5

# Check gene census
wc -l cancer_gene_census.csv
```

## Dataset Versions

### Current Release: v103

| Property | Value |
|----------|-------|
| Version | v103 |
| Release Date | 2025-11-18 |
| Total Size | ~20 GB |
| Coding Mutations | 8M+ |
| Non-coding Mutations | 10M+ |
| Tumor Samples | 1.5M+ |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| CosmicMutantExport.tsv.gz | ~5 GB | 8M | Complete mutation export |
| CosmicCodingMuts.vcf.gz | ~2 GB | 8M | Coding mutations (VCF) |
| CosmicNonCodingVariants.vcf.gz | ~3 GB | 10M | Non-coding variants |
| cancer_gene_census.csv | ~500 KB | 736 | Cancer Gene Census |
| CosmicFusionExport.tsv.gz | ~100 MB | 30K | Gene fusions |
| CosmicCompleteCNA.tsv.gz | ~2 GB | 5M | Copy number data |

### Version History

| Version | Release | New Features | Status |
|---------|---------|--------------|--------|
| v103 | 2025-11 | Pipeline upgrade | Current |
| v102 | 2025-08 | 8 new hallmark genes | Archived |
| v101 | 2025-05 | Rare cancer focus | Archived |
| v100 | 2024-11 | 100th release milestone | Archived |

### Data Categories

| Category | Records | Description |
|----------|---------|-------------|
| Somatic mutations | 8M+ | Point mutations |
| Gene Census | 736 | Driver genes (Tier 1+2) |
| Fusions | 30K+ | Gene fusion events |
| CNV | 5M+ | Copy number alterations |
| Resistance | 5K+ | Drug resistance mutations |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://cancer.sanger.ac.uk/cosmic/api/ |
| Rate Limit | 10 req/sec |
| Auth Required | Yes (COSMIC account) |
| Response Format | JSON, TSV, VCF |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major versions | Quarterly |
| Bug fixes | As needed |
| Gene Census | With major releases |

## Common Issues

- **Authentication errors**: Ensure email/password are URL-encoded if special characters
- **Download timeouts**: Use SFTP for large files; supports resume
- **Missing columns**: Schema changes between versions; check release notes
- **GRCh37 vs GRCh38**: Ensure consistent genome build; lift-over available
- **License restrictions**: Commercial use requires separate license

## License Tiers

| Tier | Access | Cost |
|------|--------|------|
| Academic | Full download access | Free |
| Non-commercial | Full download access | Free |
| Commercial | Full access + redistribution | Contact Wellcome Sanger |

## Data Categories

| Category | Description |
|----------|-------------|
| Somatic mutations | Point mutations in cancer samples |
| Gene Census | Tier 1 & 2 cancer driver genes |
| Fusions | Gene fusion events |
| CNV | Copy number variations |
| Resistance | Drug resistance mutations |
| Expression | Gene expression in cancer |
| Methylation | DNA methylation patterns |

## API Alternative

```bash
# COSMIC REST API (limited queries)
curl "https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=BRAF" \
  -H "Accept: application/json"
```
