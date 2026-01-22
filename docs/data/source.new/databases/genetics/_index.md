---
title: "Genetics Databases"
parent: ../_index.md
category: genetics
last_updated: 2026-01-22
status: draft
---

# Genetics Databases

Databases for genetic variation, clinical annotations, and pharmacogenomics.

## Database Catalog

| Database | Type | Tier | Update Frequency | Access Method | Size |
|----------|------|------|------------------|---------------|------|
| **dbSNP** | SNP/Variant Repository | 1 | Monthly | REST API, Bulk | 1.5+ billion variants |
| **ClinVar** | Clinical Variant Annotation | 1 | Monthly | REST API, FTP | 2.8+ million variants |
| **gnomAD** | Population Frequencies | 1 | Quarterly | Browser, VCF | 125,748 exomes, 71,702 genomes |
| **PharmGKB** | Pharmacogenomics | 1 | Monthly | REST API, Download | 3,000+ genes, 700+ drugs |
| **dbNSFP** | Functional Prediction | 2 | Annual | Download | 84+ million variants |
| **ExAC** | Exome Aggregation | 2 | Static (superseded by gnomAD) | Browser | 60,706 exomes |
| **GTEx** | Gene Expression | 2 | Biannual | Portal, API | 54 tissues, 948 donors |
| **UK Biobank** | Population Genomics | 3 | Ongoing | Approved Access | 500,000+ participants |
| **1000 Genomes** | Population Variation | 2 | Static (complete) | FTP | 2,504 individuals |
| **OMIM** | Genetic Disorders | 2 | Daily | Web, API (subscription) | 25,000+ entries |
| **GWAS Catalog** | Association Studies | 2 | Weekly | REST API, Download | 420,000+ associations |
| **CIViC** | Clinical Interpretations | 2 | Continuous | API, Download | 3,700+ evidence items |

## Primary Use Cases

### Tier 1 (MVP) Focus

#### dbSNP
- **Purpose**: Master variant registry
- **Use**: Map user variants to rs IDs
- **Integration**: Links to all other databases
- **Data**: Position, alleles, population frequencies

#### ClinVar
- **Purpose**: Clinical significance of variants
- **Use**: Pathogenicity classification
- **Integration**: dbSNP rs IDs, OMIM, ClinGen
- **Data**: Pathogenic, benign, VUS classifications

#### gnomAD
- **Purpose**: Population allele frequencies
- **Use**: Determine variant rarity
- **Integration**: dbSNP rs IDs
- **Data**: Allele frequencies across populations

#### PharmGKB
- **Purpose**: Drug-gene interactions
- **Use**: Pharmacogenomic recommendations
- **Integration**: dbSNP, drug databases
- **Data**: Gene-drug pairs, dosing guidelines, clinical annotations

## Data Integration Workflow

```
User Variant (VCF)
    ↓
dbSNP (rs ID mapping)
    ↓
┌─────────────┬─────────────┬─────────────┐
│             │             │             │
ClinVar       gnomAD        PharmGKB
│             │             │
Pathogenicity Frequency     Drug Response
│             │             │
└─────────────┴─────────────┴─────────────┘
                   ↓
            Combined Report
```

## Access Methods

### dbSNP
```bash
# REST API
curl "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/rs429358"

# Bulk download
wget ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/
```

### ClinVar
```bash
# E-utilities
esearch -db clinvar -query "APOE[gene]" | efetch -format vcv

# FTP download
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
```

### gnomAD
```bash
# Browser: https://gnomad.broadinstitute.org
# VCF download
gsutil -m cp -r gs://gcp-public-data--gnomad/release/4.1/ .
```

### PharmGKB
```bash
# API (requires key)
curl "https://api.pharmgkb.org/v1/data/clinicalAnnotation?symbol=APOE"

# Bulk download
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip
```

## Tier 2 Databases

### dbNSFP
- **Purpose**: Functional prediction scores (SIFT, PolyPhen, CADD, etc.)
- **Size**: 84 million variants
- **Access**: Download only (large file)
- **Use**: Predict impact of missense variants

### GTEx
- **Purpose**: Tissue-specific gene expression
- **Size**: 54 tissues, 948 donors
- **Access**: Portal, API
- **Use**: eQTL analysis, expression profiling

### GWAS Catalog
- **Purpose**: Genome-wide association studies
- **Size**: 420,000+ associations
- **Access**: REST API, download
- **Use**: Disease-variant associations

## Update Strategy

### Automated Updates (Tier 1)
- dbSNP: Monthly sync via FTP
- ClinVar: Monthly incremental updates
- gnomAD: Quarterly major releases
- PharmGKB: Monthly API refresh

### Manual Updates (Tier 2)
- dbNSFP: Annual major releases
- GTEx: Check for new releases biannually
- GWAS Catalog: Weekly updates available

## Storage Estimates

| Database | Storage Required | Format |
|----------|------------------|--------|
| dbSNP (GRCh38) | 100-200 GB | VCF.gz, JSON |
| ClinVar | 5-10 GB | VCF.gz, XML |
| gnomAD | 500 GB - 2 TB | VCF.gz |
| PharmGKB | 100 MB | JSON, TSV |
| dbNSFP | 50-100 GB | TSV.gz |
| GTEx | 50-100 GB | TSV, BED |

## Navigation

- **Parent**: [Database Sources](../_index.md)
- **Related**: [Pathways](../pathways/_index.md), [Compounds](../compounds/_index.md)
