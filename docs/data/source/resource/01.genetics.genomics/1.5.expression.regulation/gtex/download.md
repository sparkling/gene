---
id: download-gtex
title: "GTEx Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# GTEx Download Instructions

## Quick Start

```bash
# Download GTEx v8 gene expression matrix (open access)
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
```

## Prerequisites

- **wget** or **curl** for open access downloads
- **dbGaP approved access** for protected data (individual-level)
- **SRA toolkit** for sequence data
- 50GB-10TB storage depending on data types

## Access Tiers

### Open Access (No Registration)
- Summary statistics
- Gene expression matrices (aggregated)
- eQTL results
- Sample annotations

### Protected Access (dbGaP Required)
- Individual-level genotypes
- Raw sequence data (RNA-seq, WGS)
- Detailed phenotype data

## Download Methods

### Method 1: GTEx Portal (Open Access)

```bash
# Gene expression (TPM)
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

# Gene expression (read counts)
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz

# Sample annotations
wget https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

# Subject phenotypes
wget https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
```

### Method 2: eQTL Data

```bash
# Single-tissue cis-eQTLs (all tissues)
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar

# Significant variant-gene pairs (example: whole blood)
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz

# All tested pairs (for specific tissue)
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz

# sQTLs (splicing QTLs)
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_sQTL.tar
```

### Method 3: dbGaP Access (Protected Data)

```bash
# 1. Apply for access at dbGaP (phs000424)
# 2. Download repository key from dbGaP
# 3. Use SRA toolkit to download

# Configure SRA toolkit with dbGaP key
vdb-config --import ~/prj_12345.ngc

# Download specific samples
prefetch SRR1234567
fasterq-dump SRR1234567

# Download genotypes (VCF)
# Available after dbGaP approval via AnVIL or direct download
```

### Method 4: Google Cloud Storage (Bulk)

```bash
# List all GTEx files
gsutil ls gs://adult-gtex/

# Download specific directory
gsutil -m cp -r gs://adult-gtex/bulk-gex/v8/rna-seq/ ./gtex_expression/

# Download eQTL data
gsutil -m cp -r gs://adult-gtex/bulk-gex/v8/single-tissue-cis-qtl/ ./gtex_eqtl/
```

### Method 5: AnVIL (Terra Platform)

For protected data and cloud-based analysis:
1. Request access via dbGaP
2. Link to AnVIL at https://anvil.terra.bio
3. Access GTEx workspace
4. Download or analyze in cloud

## File Inventory

### Expression Data (Open)

| File | Size | Description |
|------|------|-------------|
| gene_tpm.gct.gz | ~1.5 GB | Gene TPM matrix |
| gene_reads.gct.gz | ~1.5 GB | Gene read counts |
| transcript_tpm.gct.gz | ~4 GB | Transcript TPM |
| exon_reads.gct.gz | ~5 GB | Exon-level counts |

### eQTL Data (Open)

| File | Size | Description |
|------|------|-------------|
| GTEx_Analysis_v8_eQTL.tar | ~50 GB | All cis-eQTL results |
| *.signif_variant_gene_pairs.txt.gz | 5-50 MB each | Significant pairs |
| *.egenes.txt.gz | 1-10 MB each | eGene-level results |

### Metadata (Open)

| File | Size | Description |
|------|------|-------------|
| SampleAttributesDS.txt | ~20 MB | Sample metadata |
| SubjectPhenotypesDS.txt | ~1 MB | Donor phenotypes |
| tissue_colors.txt | ~5 KB | Tissue color codes |

### Protected Data (dbGaP)

| Data Type | Size | Description |
|-----------|------|-------------|
| WGS VCF | ~500 GB | Whole genome genotypes |
| RNA-seq FASTQ | ~10 TB | Raw sequence data |
| Phenotypes | ~100 MB | Detailed clinical |

## Post-Download Processing

```bash
# Decompress expression matrix
gunzip GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

# Parse GCT format (skip header lines)
tail -n +3 gene_tpm.gct | head

# Extract specific tissue from eQTL
zcat Whole_Blood.v8.signif_variant_gene_pairs.txt.gz | \
  awk -F'\t' '$7<5e-8' > highly_sig_eqtls.txt

# Convert eQTL to BED format
zcat *.signif_variant_gene_pairs.txt.gz | \
  awk -F'\t' 'NR>1 {split($1,a,"_"); print a[1]"\t"a[2]-1"\t"a[2]"\t"$2}' \
  > eqtl_snps.bed

# Merge expression with metadata
python3 << 'EOF'
import pandas as pd

expr = pd.read_csv('gene_tpm.gct', sep='\t', skiprows=2, index_col=0)
meta = pd.read_csv('SampleAttributesDS.txt', sep='\t')
# Filter for specific tissue
blood_samples = meta[meta['SMTSD']=='Whole Blood']['SAMPID']
blood_expr = expr[blood_samples.tolist()]
blood_expr.to_csv('blood_expression.tsv', sep='\t')
EOF
```

## Verification

```bash
# Check GCT file structure
head -5 gene_tpm.gct

# Count samples and genes
awk 'NR==3 {print "Genes:", NF-2}' gene_tpm.gct
wc -l gene_tpm.gct

# Verify eQTL format
zcat Whole_Blood.v8.signif_variant_gene_pairs.txt.gz | head -5
```

## Dataset Versions

### Current Release: V10

| Property | Value |
|----------|-------|
| Version | V10 |
| Release Date | 2024-11 |
| Total Size | ~12 TB (all data types) |
| Donors | 948 |
| Tissues | 54 |
| Samples | 17,382 RNA-seq |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| gene_tpm.gct.gz | ~1.5 GB | 56K genes | TPM matrix |
| gene_reads.gct.gz | ~1.5 GB | 56K genes | Read counts |
| GTEx_Analysis_v10_eQTL.tar | ~50 GB | varies | All cis-eQTL results |
| SampleAttributesDS.txt | ~20 MB | 17K | Sample metadata |
| SubjectPhenotypesDS.txt | ~1 MB | 948 | Donor phenotypes |

### Tissue Coverage (54 tissues)

| Category | Tissues | Samples |
|----------|---------|---------|
| Brain | 13 regions | 2,500+ |
| Cardiovascular | 5 (heart, arteries) | 800+ |
| Digestive | 8 (GI tract, liver) | 1,200+ |
| Reproductive | 6 (various) | 600+ |
| Other | 22 | 12,000+ |

### Version History

| Version | Release | Donors | Tissues | Status |
|---------|---------|--------|---------|--------|
| V10 | 2024-11 | 948 | 54 | Current |
| V8 | 2020-06 | 838 | 54 | Available |
| V9 | 2022 | - | - | Metadata only |
| V7 | 2017-09 | 714 | 53 | Legacy |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://gtexportal.org/api/v2/ |
| Rate Limit | 10 req/sec |
| Auth Required | No (open access data) |
| Response Format | JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major versions | Every 2-3 years |
| Data freezes | With major releases |
| Portal updates | Continuous |

## Common Issues

- **GCT format**: First two lines are metadata; skip when parsing
- **Tissue names**: Use exact spelling from tissue list
- **Sample IDs**: Format is GTEX-DONOR-TISSUE-ALIQUOT
- **dbGaP delays**: Protected access approval takes 2-4 weeks
- **Memory issues**: Expression matrices are large; use chunked reading

## Tissue List (49 Tissues)

| Category | Tissues |
|----------|---------|
| Brain | 13 regions |
| Cardiovascular | Heart (2), Artery (3) |
| Digestive | Esophagus (3), Stomach, Colon (2), Small Intestine |
| Other | Lung, Liver, Skin (2), Adipose (2), etc. |

## dbGaP Application Tips

1. Request access at https://dbgap.ncbi.nlm.nih.gov
2. Study accession: phs000424 (GTEx)
3. Provide detailed research plan
4. Institutional signing official required
5. Allow 2-4 weeks for approval
