---
id: download-uk-biobank
title: "UK Biobank Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# UK Biobank Download Instructions

## Quick Start

```bash
# After approval, use ukbfetch to download data
./ukbfetch -a<application_id> -b<basket_id> -o<output_directory>
```

## Prerequisites

- **Approved UK Biobank application** (required - typically 3-6 months process)
- **ukbfetch/ukbconv tools** from UK Biobank
- **Data Showcase access** credentials
- Minimum 5-50TB storage depending on data types requested
- High-bandwidth internet connection

## Registration Process

### Step 1: Create Account

1. Navigate to https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access
2. Register as a researcher
3. Complete the registration questionnaire

### Step 2: Submit Application

1. Create a research application describing your project
2. Specify data fields required (phenotypes, genetics, imaging)
3. Declare institution and funding
4. Submit for review (allow 3-6 months)

### Step 3: Data Access Agreement

1. Sign the Material Transfer Agreement (MTA)
2. Complete ethics documentation
3. Pay access fees if applicable
4. Receive application ID and encryption keys

## Download Methods

### Method 1: UKB Bulk Download (ukbfetch)

```bash
# Download ukbfetch tool
wget -nd biobank.ndph.ox.ac.uk/showcase/util/ukbfetch

# Make executable
chmod +x ukbfetch

# Download phenotype data
./ukbfetch -a12345 -b67890 -d./phenotypes

# Download specific dataset files
./ukbfetch -a12345 -e12345_1234 -d./data

# Resume interrupted download
./ukbfetch -a12345 -b67890 -d./phenotypes -r
```

### Method 2: UK Biobank Research Analysis Platform (RAP)

```bash
# Use DNAnexus CLI for cloud-based access
dx login
dx select project-xxx

# List available files
dx ls /Bulk/

# Download files
dx download /Bulk/Genotype\ Results/PLINK\ format\ genetic\ data/
```

### Method 3: Genetic Data via Application

```bash
# Imputed genotypes (BGEN format)
./ukbfetch -a12345 -s1 -m  # Download chromosome 1 markers
./ukbfetch -a12345 -s1 -c  # Download chromosome 1 calls

# Exome sequencing data (PLINK format)
./ukbfetch -a12345 -dexome -f./exome_data

# WGS data (CRAM format) - RAP only
# Access via DNAnexus platform
```

### Method 4: Phenotype Data

```bash
# Download phenotype file (encoded)
./ukbfetch -a12345 -b67890 -o./phenotypes/ukb12345.enc

# Decrypt the file
./ukbunpack ukb12345.enc keyfile

# Convert to usable format
./ukbconv ukb12345.enc_ukb docs  # Generate documentation
./ukbconv ukb12345.enc_ukb r     # Convert to R format
./ukbconv ukb12345.enc_ukb csv   # Convert to CSV
```

## File Inventory

### Genetic Data

| File Type | Size | Description |
|-----------|------|-------------|
| Imputed genotypes (BGEN) | ~2TB total | 97M variants, ~500K samples |
| Exome sequencing (PLINK) | ~500GB | ~200K samples |
| WGS (CRAM) | ~50TB+ | ~200K samples |
| CNV calls | ~50GB | Copy number variants |

### Phenotype Data

| File Type | Size | Description |
|-----------|------|-------------|
| Main dataset | ~20GB | ~7,000 phenotype fields |
| Hospital episode data | ~100GB | HES records |
| Primary care data | ~50GB | GP records |
| Death registry | ~5GB | Mortality data |

### Imaging Data

| File Type | Size | Description |
|-----------|------|-------------|
| Brain MRI | ~15TB | ~50K participants |
| Cardiac MRI | ~10TB | ~50K participants |
| Abdominal MRI | ~5TB | DXA and other scans |

## Post-Download Processing

```bash
# Convert BGEN to PLINK format
plink2 --bgen ukb_imp_chr1_v3.bgen ref-first \
  --sample ukb12345_imp_chr1_v3_s487395.sample \
  --make-bed --out ukb_chr1

# Extract specific samples
plink2 --bgen ukb_imp_chr1_v3.bgen \
  --sample ukb12345_imp_chr1_v3.sample \
  --keep samples_of_interest.txt \
  --make-bed --out subset_chr1

# Convert phenotype file
./ukbconv ukb12345.enc_ukb r --ifieldid 21001,21002,31,34
```

## Verification

```bash
# Verify downloaded file integrity
md5sum -c ukb12345.md5

# Check BGEN file
bgenix -g ukb_imp_chr1_v3.bgen -list

# Validate sample file
head ukb12345_imp_chr1_v3.sample
wc -l ukb12345_imp_chr1_v3.sample
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Phenotype updates | Quarterly |
| Genetic data releases | Annually |
| Imaging data | Ongoing recruitment |
| Primary care linkage | Annually |

## Common Issues

- **Application rejection**: Ensure clear scientific rationale and data minimization
- **Encryption key errors**: Keys are application-specific; do not share between projects
- **Download timeouts**: Use ukbfetch with -r flag to resume; consider overnight downloads
- **Storage limitations**: Plan for 10-50TB; cloud storage recommended
- **Consent restrictions**: Some data fields have additional restrictions; check Data Showcase
- **RAP billing**: Research Analysis Platform has compute costs; budget accordingly

## Access Tiers

| Tier | Description | Data Included |
|------|-------------|---------------|
| Standard | Basic phenotypes + imputed genetics | Most researchers |
| Return of Results | Includes genetic results | Additional ethics |
| Imaging | Brain/Cardiac/Abdominal MRI | Subset applications |
| Primary Care | GP records linkage | Additional approval |
| WGS | Whole genome sequencing | RAP only |

## Important Notes

- Data must remain within approved secure computing environment
- Results must be returned to UK Biobank before publication
- Data cannot be transferred to third parties
- Annual renewal required for continued access
