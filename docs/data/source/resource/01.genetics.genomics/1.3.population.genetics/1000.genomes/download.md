---
id: download-1000-genomes
title: "1000 Genomes Project Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# 1000 Genomes Project Download Instructions

## Quick Start

```bash
# Download Phase 3 VCF for chromosome 1 (GRCh38)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

## Prerequisites

- wget, curl, or aspera
- ~500 GB disk space for full dataset
- tabix for indexed queries
- bcftools (recommended)

## Download Methods

### Primary: EBI FTP

```bash
# Phase 3 high-coverage (30x) - GRCh38
# Per-chromosome VCFs
for chr in {1..22} X; do
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
done
```

### Alternative: AWS S3

```bash
# Using AWS CLI (no AWS account needed for public data)
aws s3 ls --no-sign-request s3://1000genomes/

# Download specific file
aws s3 cp --no-sign-request \
  s3://1000genomes/1000G_2504_high_coverage/working/20220422_3202_phased/chr1.vcf.gz .
```

### Alternative: Aspera (Faster)

```bash
# Using Aspera Connect
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -Tr -Q -l 300M -P33001 \
  fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/data_collections/ .
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| chr*.vcf.gz | ~5-20 GB each | Per-chromosome phased VCFs |
| chr*.vcf.gz.tbi | ~1-5 MB each | Tabix indices |
| integrated_call_samples_v3.20130502.ALL.panel | 50 KB | Sample metadata |
| 20131219.populations.tsv | 10 KB | Population descriptions |

## Post-Download Processing

```bash
# Create local index if needed
tabix -p vcf 1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

# Extract specific region
tabix 1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz 1:100000-200000

# Extract population-specific samples
bcftools view -S GBR_samples.txt input.vcf.gz -Oz -o GBR_only.vcf.gz

# Calculate allele frequencies
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' input.vcf.gz
```

## Verification

```bash
# Check sample count
bcftools query -l input.vcf.gz | wc -l
# Expected: 2504 (Phase 3)

# Check variant count
bcftools view -H input.vcf.gz | wc -l
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major release | Completed (Phase 3) |
| Reference updates | GRCh38 remapping available |
| 30x high-coverage | Released 2020 |

## Integration Notes

- Widely used as imputation reference panel
- Compatible with Beagle, IMPUTE2, Minimac
- gnomAD includes 1000 Genomes samples
- Sample IDs match across releases
