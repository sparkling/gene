---
id: download-pharmgkb
title: "PharmGKB Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# PharmGKB Download Instructions

## Quick Start

```bash
# Download core annotation files (requires registration)
curl -u "username:password" -O https://api.pharmgkb.org/v1/download/file/data/annotations.zip
```

## Prerequisites

- **Free PharmGKB account** (required for downloads)
- **curl** or **wget** for command-line access
- **unzip** for extraction
- Approximately 2GB disk space

## Registration Process

### Step 1: Create Account

1. Navigate to https://www.pharmgkb.org/register
2. Provide email and institution details
3. Accept Terms of Use
4. Verify email address

### Step 2: Generate API Access

1. Log in to PharmGKB
2. Go to Account Settings
3. Generate API key for programmatic access

## Download Methods

### Method 1: Web Interface

1. Log in at https://www.pharmgkb.org
2. Navigate to Downloads: https://www.pharmgkb.org/downloads
3. Select data type
4. Click download

### Method 2: REST API

```bash
# Set credentials
PHARMGKB_USER="your_email"
PHARMGKB_PASS="your_password"

# Download clinical annotations
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o clinical_annotations.zip \
  "https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip"

# Download drug labels
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o drug_labels.zip \
  "https://api.pharmgkb.org/v1/download/file/data/drugLabels.zip"

# Download variant annotations
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o var_annotations.zip \
  "https://api.pharmgkb.org/v1/download/file/data/variantAnnotations.zip"
```

### Method 3: Bulk Data Downloads

```bash
# All annotations package
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o annotations.zip \
  "https://api.pharmgkb.org/v1/download/file/data/annotations.zip"

# Genes data
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o genes.zip \
  "https://api.pharmgkb.org/v1/download/file/data/genes.zip"

# Drugs data
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o drugs.zip \
  "https://api.pharmgkb.org/v1/download/file/data/drugs.zip"

# Relationships
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o relationships.zip \
  "https://api.pharmgkb.org/v1/download/file/data/relationships.zip"
```

### Method 4: Specific Data Categories

```bash
# CPIC guidelines
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o cpic_guidelines.zip \
  "https://api.pharmgkb.org/v1/download/file/data/cpicPairs.zip"

# FDA drug labels with PGx
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o fda_labels.zip \
  "https://api.pharmgkb.org/v1/download/file/data/drugLabels.zip"

# Pathways
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o pathways.zip \
  "https://api.pharmgkb.org/v1/download/file/data/pathways-tsv.zip"

# Automated annotations
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o automated_annotations.zip \
  "https://api.pharmgkb.org/v1/download/file/data/automatedAnnotations.zip"
```

### Method 5: VCF Annotations

```bash
# VCF-formatted PGx variants
curl -u "${PHARMGKB_USER}:${PHARMGKB_PASS}" \
  -o pharmgkb_variants.vcf.gz \
  "https://api.pharmgkb.org/v1/download/file/data/pharmgkb_vcf.zip"
```

## File Inventory

### Core Annotation Files

| File | Size | Description |
|------|------|-------------|
| clinical_annotations.tsv | ~20 MB | Clinical annotations |
| var_drug_ann.tsv | ~5 MB | Variant-drug annotations |
| var_pheno_ann.tsv | ~3 MB | Variant-phenotype annotations |

### Drug Data

| File | Size | Description |
|------|------|-------------|
| drugs.tsv | ~5 MB | Drug information |
| drug_labels.tsv | ~10 MB | FDA drug labels |
| drugLabels.tsv | ~15 MB | All regulatory labels |

### Gene Data

| File | Size | Description |
|------|------|-------------|
| genes.tsv | ~1 MB | Gene information |
| gene_relationships.tsv | ~500 KB | Gene-drug relationships |

### Guideline Data

| File | Size | Description |
|------|------|-------------|
| cpicPairs.tsv | ~500 KB | CPIC gene-drug pairs |
| dosingGuidelines.tsv | ~2 MB | Dosing guidelines |

### Pathway Data

| File | Size | Description |
|------|------|-------------|
| pathways-tsv/*.tsv | ~50 MB total | PK/PD pathways |

## Post-Download Processing

```bash
# Extract archives
unzip clinical_annotations.zip -d clinical_annotations/

# Parse clinical annotations
head -1 clinical_annotations/clinical_annotations.tsv | tr '\t' '\n' | nl

# Filter for specific gene (CYP2D6)
awk -F'\t' '$2 ~ /CYP2D6/' clinical_annotations.tsv > cyp2d6_annotations.tsv

# Get high-evidence annotations
awk -F'\t' '$11 ~ /1A|1B/' clinical_annotations.tsv > high_evidence.tsv

# Create variant lookup table
cut -f1,2,3 var_drug_ann.tsv | sort -u > variant_lookup.tsv

# Convert to VCF-friendly format
awk -F'\t' 'NR>1 {print $3"\t"$4"\t.\t"$5"\t"$6"\t.\tPASS\tGENE="$2}' \
  var_drug_ann.tsv > pgx_variants.vcf
```

## Verification

```bash
# Check file structure
head -3 clinical_annotations.tsv

# Count annotations by evidence level
cut -f11 clinical_annotations.tsv | sort | uniq -c

# Count unique gene-drug pairs
cut -f2,4 var_drug_ann.tsv | sort -u | wc -l
```

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | Continuous (2026-01) |
| Release Date | Weekly updates |
| Total Size | ~2 GB (all downloads) |
| Drugs | 750+ |
| Genes | 1,800+ |
| Annotations | 9,000+ curated |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| clinicalAnnotations.zip | ~20 MB | 2000+ | Clinical annotations |
| variantAnnotations.zip | ~5 MB | 5000+ | Variant annotations |
| drugLabels.zip | ~15 MB | 800+ | FDA/EMA drug labels |
| relationships.zip | ~10 MB | 20K+ | Gene-drug relationships |
| pathways-tsv.zip | ~50 MB | 200+ | PK/PD pathways |

### Evidence Level Distribution

| Level | Description | Count |
|-------|-------------|-------|
| 1A | CPIC/FDA guideline | 150+ |
| 1B | Published clinical | 300+ |
| 2A | VIP gene evidence | 200+ |
| 2B | Moderate evidence | 500+ |
| 3 | Low evidence | 800+ |
| 4 | Functional only | 500+ |

### Gene Coverage

| Category | Genes | Description |
|----------|-------|-------------|
| VIP Genes | 55 | Very Important Pharmacogenes |
| CYP450 | 12 | Drug metabolism enzymes |
| Transporters | 20+ | Drug transporters |
| Other | 1700+ | Additional pharmacogenes |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://api.pharmgkb.org/v1/ |
| Rate Limit | 10 req/sec |
| Auth Required | Yes (free account) |
| Response Format | JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Annotations | Weekly |
| Guidelines | As published |
| Drug labels | Monthly |
| Full database | Quarterly |

## Common Issues

- **Authentication errors**: Ensure credentials are URL-encoded
- **File format changes**: Check release notes for schema updates
- **Incomplete downloads**: Use curl with -C - to resume
- **CPIC vs PharmGKB**: CPIC guidelines are subset; PharmGKB has more annotations
- **Evidence levels**: Filter by evidence level for clinical use

## Evidence Levels

| Level | Description |
|-------|-------------|
| 1A | CPIC guideline or FDA label |
| 1B | Published clinical study |
| 2A | VIP gene with clinical evidence |
| 2B | Moderate evidence |
| 3 | Low evidence |
| 4 | Functional only |

## API Alternative

```bash
# Get gene information via API
curl "https://api.pharmgkb.org/v1/data/gene/PA124"

# Search for variants
curl "https://api.pharmgkb.org/v1/data/variant?symbol=rs1065852"

# Get drug-gene relationships
curl "https://api.pharmgkb.org/v1/data/drug/PA449015/relatedGenes"
```

## Integration with Other Resources

```bash
# Map PharmGKB IDs to other databases
# genes.tsv contains NCBI Gene ID, HGNC ID, Ensembl ID

# Annotate VCF with PharmGKB
bcftools annotate -a pharmgkb_variants.vcf.gz \
  -c CHROM,POS,REF,ALT,INFO/PGKB your_variants.vcf.gz

# Cross-reference with ClinVar
# Use variant IDs (rsIDs) for matching
```

## Related Resources

- [CPIC](../cpic/README.md) - Clinical Pharmacogenetics Implementation
- [FDA Labels](../../08.literature.knowledge/8.4.regulatory.legal/) - Regulatory information
- [DrugBank](../../02.compounds.molecules/2.2.pharmaceuticals/drugbank/) - Drug data
