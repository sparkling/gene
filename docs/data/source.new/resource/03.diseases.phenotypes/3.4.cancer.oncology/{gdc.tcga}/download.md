---
id: download-gdc-tcga
title: "GDC/TCGA Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# GDC/TCGA Download Instructions

## Quick Start

```bash
# Install GDC Data Transfer Tool
wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip

# Download open-access mutation data
./gdc-client download --manifest manifest.txt
```

## Prerequisites

- **GDC Data Transfer Tool** for bulk downloads
- **dbGaP authorization** for controlled-access data
- **eRA Commons account** for controlled data
- Approximately 1TB+ for full TCGA dataset

## Registration Requirements

| Data Type | Access Level | Requirements |
|-----------|--------------|--------------|
| Somatic mutations (MAF) | Open | None |
| Gene expression | Open | None |
| Clinical data | Open | None |
| Aligned reads (BAM) | Controlled | dbGaP approval |
| Raw sequences (FASTQ) | Controlled | dbGaP approval |

## Download Methods

### Method 1: GDC Portal Web Interface

```bash
# 1. Go to https://portal.gdc.cancer.gov/
# 2. Use Exploration or Repository tabs
# 3. Add files to cart
# 4. Download manifest or individual files

# Download manifest from cart
# Click "Download" > "Manifest"

# Use manifest with gdc-client
./gdc-client download -m gdc_manifest.txt
```

### Method 2: GDC Data Transfer Tool

```bash
# Download and install
wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip
chmod +x gdc-client

# Download by UUID (open access)
./gdc-client download <file_uuid>

# Download from manifest
./gdc-client download -m manifest.txt -d ./downloads/

# Download with token (controlled access)
./gdc-client download -m manifest.txt -t token.txt -d ./downloads/

# Parallel downloads
./gdc-client download -m manifest.txt -n 4 -d ./downloads/
```

### Method 3: GDC API

```bash
# Files endpoint
curl "https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D&format=json&size=10" \
  -o tcga_brca_files.json

# Cases endpoint
curl "https://api.gdc.cancer.gov/cases?filters=%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D&format=json&size=10" \
  -o tcga_brca_cases.json

# Get mutations (open access)
curl "https://api.gdc.cancer.gov/ssm_occurrences?filters=%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22case.project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D&format=json&size=100" \
  -o tcga_brca_mutations.json

# Download file by UUID
curl "https://api.gdc.cancer.gov/data/<file_uuid>" -o filename.maf.gz
```

### Method 4: Python API Client

```python
# Install GDC client
# pip install requests

import requests
import json

# Search for TCGA-BRCA mutation files
files_endpoint = "https://api.gdc.cancer.gov/files"

filters = {
    "op": "and",
    "content": [
        {
            "op": "=",
            "content": {
                "field": "cases.project.project_id",
                "value": "TCGA-BRCA"
            }
        },
        {
            "op": "=",
            "content": {
                "field": "data_type",
                "value": "Masked Somatic Mutation"
            }
        }
    ]
}

params = {
    "filters": json.dumps(filters),
    "fields": "file_id,file_name,data_type,file_size",
    "format": "json",
    "size": 100
}

response = requests.get(files_endpoint, params=params)
data = response.json()

for file in data['data']['hits']:
    print(f"{file['file_id']}\t{file['file_name']}\t{file['file_size']}")
```

### Method 5: TCGAbiolinks (R)

```r
# Install TCGAbiolinks
# BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

# Query TCGA-BRCA
query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

# Download
GDCdownload(query, directory = "./GDC_data")

# Prepare data
data <- GDCprepare(query)
```

### Method 6: cBioPortal Alternative

```bash
# TCGA data via cBioPortal API
curl "https://www.cbioportal.org/api/studies?keyword=tcga" -o cbioportal_tcga_studies.json

# Get mutations for specific study
curl "https://www.cbioportal.org/api/mutations/fetch?studyId=brca_tcga&entrezGeneId=672" \
  -H "Content-Type: application/json" \
  -o brca_mutations_brca1.json
```

## File Inventory

### Open Access Data

| Data Type | Files | Size |
|-----------|-------|------|
| Masked Somatic Mutations | ~400 per project | ~50 MB each |
| Gene Expression (FPKM) | ~1,000 per project | ~5 MB each |
| Copy Number Segments | ~1,000 per project | ~1 MB each |
| Clinical Data | 1 per project | ~5 MB |

### TCGA Projects

| Project | Cancer Type | Samples |
|---------|-------------|---------|
| TCGA-BRCA | Breast | 1,098 |
| TCGA-LUAD | Lung Adenocarcinoma | 585 |
| TCGA-LUSC | Lung Squamous | 504 |
| TCGA-COAD | Colon | 478 |
| TCGA-PRAD | Prostate | 500 |
| TCGA-KIRC | Kidney Clear Cell | 537 |
| TCGA-HNSC | Head and Neck | 528 |
| TCGA-LIHC | Liver | 424 |
| TCGA-GBM | Glioblastoma | 617 |
| TCGA-OV | Ovarian | 608 |

## Post-Download Processing

```bash
# Parse MAF file
python3 << 'EOF'
import pandas as pd

# Read MAF file (skip comment lines)
maf = pd.read_csv('mutations.maf.gz', sep='\t', comment='#', low_memory=False)

# Key columns
cols = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'Variant_Classification', 'Variant_Type', 'Reference_Allele',
        'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']

print(maf[cols].head())

# Count mutations per gene
gene_counts = maf['Hugo_Symbol'].value_counts()
print(f"\nTop mutated genes:\n{gene_counts.head(20)}")

# Count mutation types
print(f"\nVariant classifications:\n{maf['Variant_Classification'].value_counts()}")
EOF

# Process clinical data
python3 << 'EOF'
import pandas as pd

# Read clinical data
clinical = pd.read_csv('clinical.tsv', sep='\t')

print(f"Total cases: {len(clinical)}")
print(f"\nSurvival status:\n{clinical['vital_status'].value_counts()}")
print(f"\nTumor stage:\n{clinical['ajcc_pathologic_stage'].value_counts()}")
EOF

# Merge multi-file downloads
cat *.maf.gz > combined_mutations.maf.gz

# Convert MAF to VCF
python3 << 'EOF'
import pandas as pd

maf = pd.read_csv('mutations.maf.gz', sep='\t', comment='#', low_memory=False)

with open('mutations.vcf', 'w') as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for _, row in maf.iterrows():
        chrom = row['Chromosome']
        pos = row['Start_Position']
        ref = row['Reference_Allele']
        alt = row['Tumor_Seq_Allele2']
        gene = row['Hugo_Symbol']
        f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tGENE={gene}\n")
EOF
```

## Verification

```bash
# Check downloaded files
ls -la ./downloads/

# Verify MAF format
zcat mutations.maf.gz | head -50

# Check clinical data
head clinical.tsv

# Verify file integrity (from manifest)
md5sum -c md5sums.txt
```

## Update Schedule

| Data Type | Frequency |
|-----------|-----------|
| New GDC data | Continuous |
| TCGA (legacy) | Stable |
| Pipeline updates | Periodic |

## Common Issues

- **Controlled access**: Requires dbGaP approval (~2-4 weeks)
- **Large files**: BAM/FASTQ files are very large (50-100 GB each)
- **Sample barcodes**: Use TCGA barcode format for sample identification
- **Multiple aliquots**: Some cases have multiple aliquots/samples
- **Data versions**: GDC periodically reprocesses data; track versions

## TCGA Barcode Structure

| Position | Example | Meaning |
|----------|---------|---------|
| 1-4 | TCGA | Project |
| 5-6 | A1 | Tissue source site |
| 7-10 | A0SK | Participant |
| 11-12 | 01 | Sample type |
| 13-14 | A | Vial |
| 15-16 | 11 | Portion |
| 17 | R | Analyte |
| 18-20 | A16 | Plate |
| 21-22 | 07 | Center |

### Sample Type Codes

| Code | Type |
|------|------|
| 01 | Primary Tumor |
| 02 | Recurrent Tumor |
| 06 | Metastatic |
| 10 | Blood Normal |
| 11 | Solid Tissue Normal |

## Related Resources

- [COSMIC](../../../../01.genetics.genomics/1.6.cancer.genomics/cosmic/) - Cancer mutations database
- [cBioPortal](../../../../01.genetics.genomics/1.6.cancer.genomics/cbioportal/) - Cancer genomics visualization
- [ICGC](../../../../01.genetics.genomics/1.6.cancer.genomics/icgc/) - International cancer genomics
