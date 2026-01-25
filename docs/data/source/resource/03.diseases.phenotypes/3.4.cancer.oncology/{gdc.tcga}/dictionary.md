# GDC/TCGA - Data Dictionary

## Overview

This data dictionary documents the schema for GDC/TCGA (Genomic Data Commons / The Cancer Genome Atlas).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gdc.tcga |
| **Name** | GDC/TCGA |
| **Parent** | 3.4.cancer.oncology |
| **Total Fields** | 50+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Case

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| case_id | UUID | 1:1 | Yes | GDC unique identifier | 8a3a3f26-1c5c-... |
| submitter_id | string | 1:1 | Yes | TCGA barcode | TCGA-A1-A0SK |
| project_id | string | 1:1 | Yes | Project identifier | TCGA-BRCA |
| primary_site | string | 1:1 | Yes | Anatomic site | Breast |
| disease_type | string | 1:1 | Yes | Cancer type | Ductal and Lobular Neoplasms |
| demographics | object | 1:1 | No | Age, gender, ethnicity | {gender: female} |
| diagnoses | array | 1:N | No | Diagnosis records | [{primary_diagnosis: ...}] |
| exposures | array | 1:N | No | Exposure history | [{alcohol_history: ...}] |
| family_histories | array | 1:N | No | Family cancer history | [{relative_with_cancer: ...}] |

### Sample

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| sample_id | UUID | 1:1 | Yes | GDC unique identifier | uuid-here |
| submitter_id | string | 1:1 | Yes | TCGA sample barcode | TCGA-XX-XXXX-01A |
| sample_type | string | 1:1 | Yes | Sample type | Primary Tumor |
| tumor_descriptor | string | 1:1 | No | Tumor description | Metastatic |
| tissue_type | string | 1:1 | No | Tissue type | Tumor |
| preservation_method | string | 1:1 | No | Preservation | FFPE |
| composition | string | 1:1 | No | Sample composition | Solid Tissue |

### File

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| file_id | UUID | 1:1 | Yes | GDC unique identifier | uuid-here |
| file_name | string | 1:1 | Yes | File name | sample.bam |
| file_size | integer | 1:1 | Yes | Size in bytes | 5368709120 |
| data_type | string | 1:1 | Yes | Data type | Aligned Reads |
| data_category | string | 1:1 | Yes | Data category | Sequencing Reads |
| data_format | string | 1:1 | Yes | File format | BAM |
| experimental_strategy | string | 1:1 | Yes | Experiment type | WXS |
| access | enum | 1:1 | Yes | Access level | open, controlled |
| md5sum | string | 1:1 | Yes | MD5 checksum | abc123... |

### Somatic Mutation (MAF)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Hugo_Symbol | string | 1:1 | Yes | Gene symbol | TP53 |
| Chromosome | string | 1:1 | Yes | Chromosome | 17 |
| Start_Position | integer | 1:1 | Yes | Start position | 7577121 |
| End_Position | integer | 1:1 | Yes | End position | 7577121 |
| Reference_Allele | string | 1:1 | Yes | Reference allele | G |
| Tumor_Seq_Allele2 | string | 1:1 | Yes | Alternate allele | A |
| Variant_Classification | string | 1:1 | Yes | Variant type | Missense_Mutation |
| Variant_Type | string | 1:1 | Yes | Mutation type | SNP |
| Tumor_Sample_Barcode | string | 1:1 | Yes | TCGA barcode | TCGA-A1-A0SK-01A |
| HGVSc | string | 1:1 | No | HGVS coding | c.742C>T |
| HGVSp | string | 1:1 | No | HGVS protein | p.R248W |
| t_depth | integer | 1:1 | No | Tumor depth | 150 |
| t_alt_count | integer | 1:1 | No | Tumor alt count | 45 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| UUID | 8-4-4-4-12 hex | 8a3a3f26-1c5c-... | GDC unique ID |
| TCGA Barcode | TCGA-XX-XXXX | TCGA-A1-A0SK | Participant ID |
| Sample Barcode | TCGA-XX-XXXX-##A | TCGA-A1-A0SK-01A | Sample ID |
| Project ID | TCGA-XXXX | TCGA-BRCA | Project ID |
| dbSNP | rs[0-9]+ | rs121912651 | Variant ID |
| COSMIC | COSM[0-9]+ | COSM10656 | Cancer variant |

---

## Enumerations

### Sample Type Codes

| Code | Description |
|------|-------------|
| 01 | Primary Solid Tumor |
| 02 | Recurrent Solid Tumor |
| 03 | Primary Blood Derived Cancer |
| 06 | Metastatic |
| 10 | Blood Derived Normal |
| 11 | Solid Tissue Normal |
| 14 | Bone Marrow Normal |
| 20 | Cell Line Derived |

### Variant Classifications

| Classification | Description |
|----------------|-------------|
| Missense_Mutation | Amino acid change |
| Nonsense_Mutation | Premature stop |
| Silent | Synonymous |
| Frame_Shift_Del | Frameshift deletion |
| Frame_Shift_Ins | Frameshift insertion |
| In_Frame_Del | In-frame deletion |
| In_Frame_Ins | In-frame insertion |
| Splice_Site | Splice site mutation |
| Nonstop_Mutation | Stop codon lost |

### Data Types

| Type | Description | Access |
|------|-------------|--------|
| Aligned Reads | BAM files | Controlled |
| Raw Sequencing | FASTQ | Controlled |
| Gene Expression | TPM/FPKM | Open |
| Somatic Mutations | MAF | Open |
| Copy Number | SEG/TSV | Open |
| Methylation | Beta values | Open |
| Clinical | Demographics | Open |

### Experimental Strategies

| Strategy | Description |
|----------|-------------|
| WGS | Whole Genome Sequencing |
| WXS | Whole Exome Sequencing |
| RNA-Seq | RNA Sequencing |
| miRNA-Seq | MicroRNA Sequencing |
| Methylation Array | DNA Methylation |
| Genotyping Array | SNP Genotyping |

### TCGA Cancer Types

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

### Access Levels

| Level | Description | Requirements |
|-------|-------------|--------------|
| open | Publicly available | None |
| controlled | Restricted access | dbGaP authorization |

---

## Entity Relationships

### Project to Cases
- **Cardinality:** 1:N
- **Description:** Projects contain multiple cases
- **Key Fields:** project_id, case_id

### Case to Samples
- **Cardinality:** 1:N
- **Description:** Cases have multiple samples
- **Key Fields:** case_id, sample_id

### Sample to Files
- **Cardinality:** 1:N
- **Description:** Samples have multiple data files
- **Key Fields:** sample_id, file_id

### Case to Diagnoses
- **Cardinality:** 1:N
- **Description:** Cases can have multiple diagnoses
- **Key Fields:** case_id, diagnosis_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| TCGA | The Cancer Genome Atlas | NCI program |
| GDC | Genomic Data Commons | Data portal |
| MAF | Mutation Annotation Format | Variant format |
| UUID | Universally Unique Identifier | GDC ID format |
| TSS | Tissue Source Site | Barcode component |
| WXS | Whole Exome Sequencing | Experiment type |
| WGS | Whole Genome Sequencing | Experiment type |
| CNV | Copy Number Variation | Data type |
| FFPE | Formalin-Fixed Paraffin-Embedded | Preservation |
| dbGaP | Database of Genotypes and Phenotypes | Access control |
| TPM | Transcripts Per Million | Expression unit |
| FPKM | Fragments Per Kilobase Million | Expression unit |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant annotation |
| COSMIC | COSM ID | Cancer variants |
| Ensembl | Gene ID | Gene annotation |
| HGNC | Gene symbol | Gene nomenclature |
| dbGaP | Study ID | Access control |
| PubMed | PMID | Literature |
| cBioPortal | Study ID | Visualization |

---

## Data Quality Notes

1. **GDC Cases:** 80,000+ cancer cases
2. **TCGA Samples:** 20,000+ samples across 33 cancer types
3. **Somatic Mutations:** 3M+ open-access mutations
4. **Data Volume:** 2.5+ petabytes total
5. **Open Access:** MAF, CNV, expression, clinical
6. **Controlled Access:** Aligned reads via dbGaP
7. **REST API:** Full programmatic access
8. **Regular Updates:** Data Release 39.0+

