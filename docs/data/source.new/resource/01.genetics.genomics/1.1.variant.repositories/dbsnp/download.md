---
id: download-dbsnp
title: "dbSNP Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# dbSNP Download Instructions

## Quick Start

```bash
# Download latest dbSNP VCF (GRCh38)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi
```

## Prerequisites

- **wget** or **curl** for downloads
- **bcftools/tabix** for VCF manipulation
- **rsync** for efficient bulk transfer
- 100-500GB disk space for full database

## No Registration Required

dbSNP data is public domain and freely available.

## Download Methods

### Method 1: Latest Release VCF (Recommended)

```bash
# GRCh38 (RefSeq chromosome names)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

# GRCh38 (UCSC chromosome names: chr1, chr2, etc.)
# Note: Convert RefSeq to UCSC names after download
bcftools annotate --rename-chrs chr_name_conv.txt \
  GCF_000001405.40.gz -Oz -o dbsnp_grch38_ucsc.vcf.gz

# GRCh37 (legacy)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
```

### Method 2: Per-Chromosome Downloads

```bash
# List available chromosome files
curl -s https://ftp.ncbi.nih.gov/snp/latest_release/VCF/ | grep -o 'GCF.*chr[0-9XYM]*\.gz'

# Download specific chromosome
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40_chr1.gz
```

### Method 3: ALFA (Allele Frequency Aggregator)

```bash
# Download ALFA frequencies
wget https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/freq.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/freq.vcf.gz.tbi

# ALFA by population
wget https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/SAMN10492695.freq.vcf.gz  # European
wget https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/SAMN10492696.freq.vcf.gz  # African
```

### Method 4: JSON Batch Downloads

```bash
# Download JSON format (alternative to VCF)
wget https://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-merged.json.bz2
wget https://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chr1.json.bz2
```

### Method 5: rsync for Full Mirror

```bash
# Mirror entire latest release
rsync -av rsync://ftp.ncbi.nih.gov/snp/latest_release/ ./dbsnp_latest/

# Mirror VCF only
rsync -av rsync://ftp.ncbi.nih.gov/snp/latest_release/VCF/ ./dbsnp_vcf/
```

### Method 6: E-utilities API

```bash
# Query specific RS numbers
efetch -db snp -id 334 -format xml > rs334.xml

# Batch query
echo -e "rs334\nrs1805007\nrs429358" | \
  while read rs; do
    efetch -db snp -id ${rs#rs} -format json >> snp_batch.json
  done
```

## File Inventory

### VCF Files

| File | Size | Description |
|------|------|-------------|
| GCF_000001405.40.gz | ~50 GB | All dbSNP (GRCh38) |
| GCF_000001405.25.gz | ~45 GB | All dbSNP (GRCh37) |
| GCF_*_chr*.gz | 1-5 GB each | Per-chromosome files |

### ALFA Frequency Files

| File | Size | Description |
|------|------|-------------|
| freq.vcf.gz | ~15 GB | All ALFA frequencies |
| SAMN*.freq.vcf.gz | ~2-5 GB | Population-specific |

### JSON Files

| File | Size | Description |
|------|------|-------------|
| refsnp-merged.json.bz2 | ~5 GB | Merged RS records |
| refsnp-chr*.json.bz2 | 500MB-2GB | Per-chromosome JSON |

### Supplementary Files

| File | Size | Description |
|------|------|-------------|
| RsMergeArch.bcp.gz | ~200 MB | RS merge history |
| SNPChrPosOnRef_*.bcp.gz | ~5 GB | Position mapping |

## Post-Download Processing

```bash
# Index VCF if needed
tabix -p vcf GCF_000001405.40.gz

# Convert RefSeq to UCSC chromosome names
cat > chr_conv.txt << 'EOF'
NC_000001.11 chr1
NC_000002.12 chr2
NC_000003.12 chr3
# ... add all chromosomes
NC_000024.10 chrY
NC_012920.1 chrM
EOF

bcftools annotate --rename-chrs chr_conv.txt GCF_000001405.40.gz \
  -Oz -o dbsnp_b156_GRCh38.vcf.gz

# Extract common variants (AF > 0.01)
bcftools view -i 'INFO/COMMON=1' GCF_000001405.40.gz \
  -Oz -o dbsnp_common.vcf.gz

# Extract specific region
bcftools view -r NC_000017.11:43044295-43125483 GCF_000001405.40.gz \
  -Oz -o brca1_snps.vcf.gz

# Create rsID lookup table
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' GCF_000001405.40.gz \
  > rsid_positions.tsv
```

## Verification

```bash
# Check VCF header
bcftools view -h GCF_000001405.40.gz | head -30

# Count variants
bcftools view -H GCF_000001405.40.gz | wc -l

# Verify specific RS
bcftools view -i 'ID="rs334"' GCF_000001405.40.gz

# Check ALFA populations
bcftools view -h freq.vcf.gz | grep "##SAMPLE"
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Build releases | Every 2-3 months |
| ALFA updates | Quarterly |
| Merge updates | Weekly |

## Common Issues

- **Chromosome naming**: RefSeq uses NC_* names; convert to chr* for most tools
- **Large file size**: Download by chromosome for memory-constrained systems
- **RS merging**: Check RsMergeArch for historical RS ID changes
- **Missing MAF**: Not all SNPs have ALFA frequencies; some have only TOPMED or gnomAD
- **Multiallelic sites**: Some RS IDs map to multiple alleles at same position

## dbSNP INFO Fields

| Field | Description |
|-------|-------------|
| RS | dbSNP RS identifier |
| GENEINFO | Gene symbol and ID |
| COMMON | Common variant flag (AF>0.01) |
| FREQ | Allele frequencies by study |
| CLNVI | ClinVar variation ID |
| VC | Variant class |

## ALFA Population Codes

| Code | Population |
|------|------------|
| SAMN10492695 | European |
| SAMN10492696 | African |
| SAMN10492697 | East Asian |
| SAMN10492698 | South Asian |
| SAMN10492699 | Latin American |
| SAMN10492703 | Other |

## Integration with Other Resources

```bash
# Annotate your VCF with dbSNP RS IDs
bcftools annotate -a GCF_000001405.40.gz \
  -c ID your_variants.vcf.gz -Oz -o annotated.vcf.gz

# Extract ALFA frequencies
bcftools +fill-tags GCF_000001405.40.gz -Oz -o with_af.vcf.gz -- -t AF

# Cross-reference with ClinVar
bcftools isec -n=2 dbsnp.vcf.gz clinvar.vcf.gz -o overlap/
```
