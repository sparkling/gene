---
id: download-clinvar
title: "ClinVar Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# ClinVar Download Instructions

## Quick Start

```bash
# Download latest ClinVar VCF (GRCh38)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

## Prerequisites

- **wget** or **curl** for downloads
- **tabix/bcftools** for VCF manipulation
- **Entrez Direct** (optional) for E-utilities access
- Approximately 5-10GB disk space

## No Registration Required

ClinVar data is in the public domain and freely available without registration.

## Download Methods

### Method 1: FTP VCF Files (Recommended)

```bash
# GRCh38 VCF (current standard)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# GRCh37 VCF (legacy)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi

# Dated archive (for reproducibility)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20240120.vcf.gz
```

### Method 2: XML Full Archive

```bash
# Complete XML archive (most detailed)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

# Variation archive (VCV records)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_00-latest.xml.gz
```

### Method 3: Tab-Delimited Summary Files

```bash
# Variant summary (most commonly used)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# Submission summary
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz

# Gene-condition associations
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_condition_source_id

# Variation-allele mapping
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variation_allele.txt.gz
```

### Method 4: E-utilities API

```bash
# Install Entrez Direct
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Search ClinVar
esearch -db clinvar -query "BRCA1[gene] AND pathogenic[clinical_significance]" | \
  efetch -format xml > brca1_pathogenic.xml

# Get specific VCV record
efetch -db clinvar -id VCV000000123 -format xml > vcv000000123.xml
```

### Method 5: Specific Data Subsets

```bash
# Pathogenic variants only
bcftools view -i 'CLNSIG="Pathogenic"' clinvar.vcf.gz -o pathogenic.vcf.gz -Oz

# Variants with review status
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt

# Gene-specific downloads via rsync
rsync -av rsync://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ ./clinvar_vcf/
```

## File Inventory

### VCF Files

| File | Size | Description |
|------|------|-------------|
| clinvar.vcf.gz | ~150 MB | All variants (VCF) |
| clinvar.vcf.gz.tbi | ~1 MB | Tabix index |
| clinvar_papu.vcf.gz | ~80 MB | Pathogenic/likely pathogenic only |

### XML Files

| File | Size | Description |
|------|------|-------------|
| ClinVarFullRelease_*.xml.gz | ~2 GB | Complete XML database |
| ClinVarVariationRelease_*.xml.gz | ~1.5 GB | Variation-centric XML |
| RCVAccession_*.xml.gz | ~3 GB | RCV records |

### Tab-Delimited Files

| File | Size | Description |
|------|------|-------------|
| variant_summary.txt.gz | ~200 MB | Variant summary table |
| submission_summary.txt.gz | ~500 MB | All submissions |
| gene_condition_source_id | ~5 MB | Gene-disease mapping |
| variation_allele.txt.gz | ~100 MB | Variation-allele mapping |
| hgvs4variation.txt.gz | ~300 MB | HGVS nomenclature |

## Post-Download Processing

```bash
# Decompress and inspect
gunzip -k clinvar.vcf.gz
bcftools view -h clinvar.vcf | tail -50

# Extract pathogenic variants
bcftools view -i 'INFO/CLNSIG ~ "Pathogenic"' clinvar.vcf.gz | \
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/GENEINFO\n' \
  > pathogenic_variants.tsv

# Filter by review status (expert panel or higher)
bcftools view -i 'INFO/CLNREVSTAT ~ "reviewed_by_expert_panel|practice_guideline"' \
  clinvar.vcf.gz -o expert_reviewed.vcf.gz -Oz

# Convert to BED format
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\t%INFO/CLNSIG\n' clinvar.vcf.gz > clinvar.bed

# Parse variant_summary
zcat variant_summary.txt.gz | \
  awk -F'\t' '$7=="GRCh38" && $6~/Pathogenic/' > pathogenic_grch38.tsv
```

## Verification

```bash
# Check VCF integrity
bcftools view -h clinvar.vcf.gz | grep "##fileDate"

# Count variants
bcftools view -H clinvar.vcf.gz | wc -l

# Verify tabix index
tabix -l clinvar.vcf.gz

# Check variant_summary columns
zcat variant_summary.txt.gz | head -1 | tr '\t' '\n' | nl
```

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | 2026-01 (Weekly) |
| Release Date | 2026-01-20 |
| Total Size | ~2.5 GB (XML), ~177 MB (VCF) |
| Variant Records | 2.5M+ VCV records |
| Submissions | 4M+ SCV records |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| clinvar.vcf.gz | ~150 MB | 2.5M | All variants (VCF format) |
| ClinVarFullRelease.xml.gz | ~2.5 GB | 2.5M | Complete XML database |
| variant_summary.txt.gz | ~200 MB | 2.5M | Tab-delimited summary |
| submission_summary.txt.gz | ~500 MB | 4M | All submissions |

### Previous Versions

| Version | Release | Size | Status |
|---------|---------|------|--------|
| 2025-12 | 2025-12-16 | ~2.4 GB | Archived |
| 2025-11 | 2025-11-18 | ~2.4 GB | Archived |
| 2025-10 | 2025-10-21 | ~2.3 GB | Archived |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ |
| Rate Limit | 3 req/sec (10 with API key) |
| Auth Required | No (API key recommended) |
| Response Format | XML, JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Weekly VCF | Every Monday |
| Monthly archive | First of month |
| XML release | Weekly |

## Common Issues

- **Version mismatch**: Ensure GRCh37/GRCh38 matches your pipeline
- **CLNSIG encoding**: Pathogenic significance may have multiple values separated by "|"
- **Missing variants**: Some variants only in XML (complex alleles, haplotypes)
- **Star alleles**: Pharmacogenomic star alleles may have special handling
- **Conflicting interpretations**: Check CLNREVSTAT for conflicts

## CLNSIG Values Reference

| Value | Description |
|-------|-------------|
| Pathogenic | Disease-causing |
| Likely_pathogenic | Probably disease-causing |
| Uncertain_significance | VUS |
| Likely_benign | Probably not disease-causing |
| Benign | Not disease-causing |
| Conflicting_interpretations | Multiple labs disagree |

## CLNREVSTAT Values Reference

| Value | Stars | Description |
|-------|-------|-------------|
| practice_guideline | 4 | Practice guideline |
| reviewed_by_expert_panel | 3 | Expert panel review |
| criteria_provided,_multiple_submitters | 2 | Multiple concordant submitters |
| criteria_provided,_single_submitter | 1 | Single submitter with criteria |
| no_assertion_criteria_provided | 0 | No criteria |

## Integration with Other Resources

```bash
# Add gnomAD frequencies
bcftools annotate -a gnomad.vcf.gz -c AF clinvar.vcf.gz -o clinvar_annotated.vcf.gz -Oz

# Cross-reference with dbSNP
bcftools query -f '%ID\n' clinvar.vcf.gz | grep "^rs" > clinvar_rsids.txt
```
