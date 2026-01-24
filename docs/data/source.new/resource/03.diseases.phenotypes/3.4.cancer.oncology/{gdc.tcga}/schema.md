---
id: schema-gdc.tcga
title: "GDC/TCGA Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: migrated
tags: [schema, database, cancer, genomics, tcga, nci, somatic-mutations]
---

# GDC/TCGA Schema Documentation

**Document ID:** SCHEMA-GDC-TCGA
**Version:** 1.0
**Source Version:** GDC Data Release 39.0

---

## TL;DR

GDC hosts 80K+ cancer cases including the landmark TCGA dataset (20K+ samples, 33 cancer types). Data model uses UUIDs for entities, TCGA barcodes for samples. Open data (MAF, CNV, expression) is freely available; aligned reads require dbGaP authorization.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Total GDC Cases | 80,000+ | All projects |
| TCGA Samples | 20,000+ | TCGA program |
| Cancer Types (TCGA) | 33 | Primary tumors |
| Somatic Mutations | 3M+ | Open MAF files |
| Projects in GDC | 70+ | Active projects |
| Data Volume | 2.5+ PB | All data types |

---

## Entity Relationship Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                         PROJECT                                  │
│  project_id (TCGA-BRCA), name, program, disease_type            │
└─────────────────────────────────────────────────────────────────┘
         │
         │ 1:N
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                          CASE                                    │
│  case_id (UUID), submitter_id (TCGA-XX-XXXX)                    │
│  demographics, diagnoses, exposures, family_histories           │
└─────────────────────────────────────────────────────────────────┘
         │
         │ 1:N
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                         SAMPLE                                   │
│  sample_id (UUID), submitter_id (TCGA-XX-XXXX-01A)              │
│  sample_type, tumor_descriptor, tissue_type                     │
└─────────────────────────────────────────────────────────────────┘
         │
         │ 1:N
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                        ALIQUOT                                   │
│  aliquot_id (UUID), submitter_id (TCGA-XX-XXXX-01A-01D)        │
│  analyte_type, concentration, amount                            │
└─────────────────────────────────────────────────────────────────┘
         │
         │ 1:N
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                         FILE                                     │
│  file_id (UUID), file_name, data_type, data_category           │
│  experimental_strategy, data_format, access (open/controlled)  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Core Tables/Entities

### Case

**Description:** Individual patient/participant

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| case_id | UUID | Yes | GDC unique identifier |
| submitter_id | string | Yes | TCGA barcode (TCGA-XX-XXXX) |
| project_id | string | Yes | Project (e.g., TCGA-BRCA) |
| primary_site | string | Yes | Anatomic site |
| disease_type | string | Yes | Cancer type |
| demographics | object | No | Age, gender, ethnicity |
| diagnoses | list | No | Diagnosis records |
| exposures | list | No | Exposure history |
| family_histories | list | No | Family cancer history |

### Sample

**Description:** Biological specimen from a case

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| sample_id | UUID | Yes | GDC unique identifier |
| submitter_id | string | Yes | TCGA barcode with sample code |
| sample_type | string | Yes | Primary Tumor, Blood Normal, etc. |
| tumor_descriptor | string | No | Metastatic, Recurrence, etc. |
| tissue_type | string | No | Tumor, Normal, etc. |
| preservation_method | string | No | FFPE, Frozen, etc. |
| composition | string | No | Solid Tissue, Blood, etc. |

### File

**Description:** Data file associated with case/sample

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| file_id | UUID | Yes | GDC unique identifier |
| file_name | string | Yes | File name |
| file_size | integer | Yes | Size in bytes |
| data_type | string | Yes | Aligned Reads, Gene Expression, etc. |
| data_category | string | Yes | Sequencing Reads, Transcriptome, etc. |
| data_format | string | Yes | BAM, TSV, MAF, etc. |
| experimental_strategy | string | Yes | WXS, RNA-Seq, etc. |
| access | enum | Yes | open or controlled |
| md5sum | string | Yes | MD5 checksum |

### Somatic Mutation (MAF)

**Description:** Somatic variant from mutation calling

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Hugo_Symbol | string | Yes | Gene symbol |
| Chromosome | string | Yes | Chromosome |
| Start_Position | integer | Yes | Start position |
| End_Position | integer | Yes | End position |
| Reference_Allele | string | Yes | Reference allele |
| Tumor_Seq_Allele2 | string | Yes | Alternate allele |
| Variant_Classification | string | Yes | Missense, Nonsense, etc. |
| Variant_Type | string | Yes | SNP, INS, DEL |
| Tumor_Sample_Barcode | string | Yes | TCGA barcode |
| HGVSc | string | No | HGVS coding notation |
| HGVSp | string | No | HGVS protein notation |
| t_depth | integer | No | Tumor depth |
| t_ref_count | integer | No | Tumor ref count |
| t_alt_count | integer | No | Tumor alt count |
| n_depth | integer | No | Normal depth |

---

## TCGA Barcode Format

### Structure

```
TCGA-XX-XXXX-01A-01D-A123-01
│    │  │    │  │  │  │    │
│    │  │    │  │  │  │    └── Data level (01=sequence)
│    │  │    │  │  │  └────── Plate ID
│    │  │    │  │  └───────── Portion/Analyte
│    │  │    │  └──────────── Vial
│    │  │    └─────────────── Sample type (01=Primary, 11=Normal)
│    │  └──────────────────── Participant ID
│    └─────────────────────── Tissue Source Site (TSS)
└──────────────────────────── Project
```

### Sample Type Codes

| Code | Description |
|------|-------------|
| 01 | Primary Solid Tumor |
| 02 | Recurrent Solid Tumor |
| 03 | Primary Blood Derived Cancer |
| 06 | Metastatic |
| 10 | Blood Derived Normal |
| 11 | Solid Tissue Normal |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /cases | GET | Query cases |
| /files | GET | Query files |
| /projects | GET | Query projects |
| /annotations | GET | Query annotations |
| /data/{uuid} | GET | Download file |
| /slicing/view/{uuid} | GET | BAM slicing |

### Query Examples

```bash
# Get TCGA-BRCA cases
curl "https://api.gdc.cancer.gov/cases?filters=%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D&size=10"

# Get open MAF files
curl "https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22data_format%22%2C%22value%22%3A%22MAF%22%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22access%22%2C%22value%22%3A%22open%22%7D%7D%5D%7D"
```

---

## Data Formats

| Format | Description |
|--------|-------------|
| MAF | Mutation Annotation Format (somatic variants) |
| VCF | Variant Call Format |
| BAM | Binary Alignment Map |
| TSV | Tab-separated values |
| Parquet | Columnar format for large files |
| Encoding | UTF-8 |

---

## Sample Records

### Case (JSON)

```json
{
  "case_id": "8a3a3f26-1c5c-4c7c-9d9a-1b2c3d4e5f6a",
  "submitter_id": "TCGA-A1-A0SK",
  "project": {
    "project_id": "TCGA-BRCA",
    "name": "Breast Invasive Carcinoma"
  },
  "primary_site": "Breast",
  "disease_type": "Ductal and Lobular Neoplasms",
  "diagnoses": [
    {
      "diagnosis_id": "uuid-here",
      "primary_diagnosis": "Infiltrating duct carcinoma, NOS",
      "age_at_diagnosis": 22645,
      "vital_status": "alive",
      "tumor_stage": "stage iia"
    }
  ],
  "demographic": {
    "gender": "female",
    "race": "white",
    "ethnicity": "not hispanic or latino"
  }
}
```

### MAF Record (TSV)

```tsv
Hugo_Symbol	Chromosome	Start_Position	End_Position	Reference_Allele	Tumor_Seq_Allele2	Variant_Classification	Variant_Type	Tumor_Sample_Barcode	HGVSc	HGVSp_Short
TP53	17	7577121	7577121	G	A	Missense_Mutation	SNP	TCGA-A1-A0SK-01A-11D-A10Y-09	c.742C>T	p.R248W
```

### Gene Expression (TSV)

```tsv
gene_id	gene_name	gene_type	unstranded	stranded_first	stranded_second	tpm_unstranded	fpkm_unstranded
ENSG00000000003	TSPAN6	protein_coding	1245	623	622	12.45	15.67
```

---

## Data Types Available

| Data Type | Format | Access | Description |
|-----------|--------|--------|-------------|
| Somatic Mutations | MAF | Open | Filtered somatic variants |
| Copy Number | SEG/TSV | Open | Gene-level CNV |
| Gene Expression | TSV | Open | RNA-seq FPKM/TPM |
| miRNA Expression | TSV | Open | miRNA counts |
| Methylation | TSV | Open | Beta values |
| Clinical | TSV/JSON | Open | Patient demographics, outcomes |
| Aligned Reads | BAM | Controlled | WGS/WXS alignments |
| Raw Sequences | FASTQ | Controlled | Unaligned reads |

---

## TCGA Cancer Types

| Code | Cancer Type | Samples |
|------|-------------|---------|
| BRCA | Breast invasive carcinoma | 1,098 |
| LUAD | Lung adenocarcinoma | 585 |
| LUSC | Lung squamous cell carcinoma | 504 |
| PRAD | Prostate adenocarcinoma | 500 |
| THCA | Thyroid carcinoma | 507 |
| COAD | Colon adenocarcinoma | 478 |
| KIRC | Kidney renal clear cell | 537 |
| HNSC | Head and neck squamous cell | 528 |
| LIHC | Liver hepatocellular | 377 |
| STAD | Stomach adenocarcinoma | 443 |

---

## Glossary

| Term | Definition |
|------|------------|
| TCGA | The Cancer Genome Atlas |
| GDC | Genomic Data Commons |
| MAF | Mutation Annotation Format |
| UUID | Universally Unique Identifier |
| TSS | Tissue Source Site |
| WXS | Whole Exome Sequencing |
| WGS | Whole Genome Sequencing |
| CNV | Copy Number Variation |
| dbGaP | Database of Genotypes and Phenotypes |
| Controlled Access | Requires IRB approval via dbGaP |

---

## References

1. https://portal.gdc.cancer.gov/
2. https://docs.gdc.cancer.gov/
3. https://api.gdc.cancer.gov/
4. Cancer Genome Atlas Research Network publications

---

## See Also

- [GDC/TCGA Overview](./_index.md)
- [GDC/TCGA Download Instructions](./download.md)
- [cBioPortal](../../cbioportal/_index.md) - Cancer genomics visualization
- [COSMIC](../../../01.genetics.genomics/1.6.cancer.genomics/cosmic/_index.md) - Cancer mutation database
