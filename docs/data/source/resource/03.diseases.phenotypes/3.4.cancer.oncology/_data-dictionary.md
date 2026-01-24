# 3.4 Cancer/Oncology & Drug-Disease Targets - Data Dictionary

## Overview

This subcategory contains cancer genomics data including somatic mutations from tumor samples and drug-gene interaction data for oncology applications.

**Data Sources:** GDC/TCGA, DGIdb

---

## Unified Fields

These fields are harmonized across multiple data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `gene_symbol` | string | Required (1:1) | HGNC gene symbol | GDC/TCGA, DGIdb | `EGFR`, `TP53`, `BRAF` |
| `sample_id` | string | Optional (1:1) | Sample or patient identifier | GDC/TCGA | `TCGA-A1-A0SB-01A` |
| `cancer_type` | string | Optional (1:1) | Cancer type or project code | GDC/TCGA | `TCGA-BRCA`, `LUAD`, `COAD` |
| `drug_name` | string | Optional (1:1) | Drug or compound name | DGIdb | `Erlotinib`, `Vemurafenib` |
| `interaction_types` | array[string] | Optional (1:N) | Types of drug-gene interaction | DGIdb | `["inhibitor", "antagonist"]` |

---

## Source-Specific Fields

### GDC/TCGA

#### Case and Sample Information

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `case_id` | string (UUID) | Optional | GDC unique case identifier | UUID | `a1b2c3d4-e5f6-7890-abcd-ef1234567890` |
| `submitter_id` | string | Optional | TCGA barcode for patient/sample | `TCGA-[A-Z]{2}-[A-Z0-9]{4}` | `TCGA-A1-A0SB` |
| `project_id` | string | Optional | Project identifier | - | `TCGA-BRCA`, `TCGA-LUAD` |
| `primary_site` | string | Optional | Anatomic site of tumor | - | `Breast`, `Lung`, `Colon` |
| `sample_type` | enum | Optional | Type of biological specimen | `Primary Tumor`, `Recurrent Tumor`, `Metastatic`, `Blood Derived Normal`, `Solid Tissue Normal` | `Primary Tumor` |

#### Mutation Information

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `chromosome` | string | Optional | Chromosome of mutation | `chr1`, `chr17`, `chrX` |
| `start_position` | integer | Optional | Mutation start position (GRCh38) | `43044295`, `7577120` |
| `end_position` | integer | Optional | Mutation end position (GRCh38) | `43044295`, `7577121` |
| `reference_allele` | string | Optional | Reference allele | `A`, `G`, `CT` |
| `tumor_allele` | string | Optional | Alternate/tumor allele | `T`, `C`, `-` |
| `variant_classification` | enum | Optional | Functional classification of variant | `Missense_Mutation`, `Nonsense_Mutation`, `Frame_Shift_Del`, `Frame_Shift_Ins`, `Splice_Site`, `Silent` | `Missense_Mutation` |
| `variant_type` | enum | Optional | Type of sequence variant | `SNP`, `INS`, `DEL`, `DNP`, `TNP` | `SNP` |
| `hgvsc` | string | Optional | HGVS coding notation | `c.2573T>G` |
| `hgvsp` | string | Optional | HGVS protein notation | `p.Leu858Arg` |

#### Sequencing Metrics

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `t_depth` | integer | Optional | Tumor sequencing depth | `150`, `245` |
| `t_alt_count` | integer | Optional | Tumor alternate allele count | `45`, `78` |

#### Data File Information

| Field Name | Data Type | Cardinality | Description | Allowed Values | Example Values |
|------------|-----------|-------------|-------------|----------------|----------------|
| `data_type` | string | Optional | Type of data file | - | `Aligned Reads`, `Gene Expression Quantification`, `Masked Somatic Mutation` |
| `experimental_strategy` | enum | Optional | Experimental methodology | `WXS`, `WGS`, `RNA-Seq`, `miRNA-Seq`, `Methylation Array` | `WXS` |
| `access` | enum | Optional | Data access level | `open`, `controlled` | `open` |

---

### DGIdb

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `concept_id` | string | Optional | Normalized identifier for genes or drugs | `HGNC:3236`, `ChEMBL:CHEMBL553` |
| `interaction_score` | float | Optional | Strength of drug-gene interaction | `0.85`, `0.42` |
| `approved` | boolean | Optional | FDA approval status of drug | `true`, `false` |
| `immunotherapy` | boolean | Optional | Whether drug is immunotherapy | `true`, `false` |
| `antineoplastic` | boolean | Optional | Anti-cancer activity flag | `true`, `false` |
| `gene_categories` | array[string] | Optional | Gene categories | `["KINASE", "TUMOR SUPPRESSOR", "DRUGGABLE GENOME"]` |
| `dgidb_sources` | array[string] | Optional | Data sources supporting interaction | `["ChEMBL", "DrugBank", "PharmGKB"]` |
| `publications` | array[integer] | Optional | Supporting publication PMIDs | `[12345678, 23456789]` |

---

## Source Field Mappings

### GDC/TCGA Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `case_id` | `case_id` |
| `submitter_id` | `submitter_id` |
| `project_id` | `project_id` |
| `primary_site` | `primary_site` |
| `sample_type` | `sample_type` |
| `Hugo_Symbol` | `gene_symbol` |
| `Chromosome` | `chromosome` |
| `Start_Position` | `start_position` |
| `End_Position` | `end_position` |
| `Reference_Allele` | `reference_allele` |
| `Tumor_Seq_Allele2` | `tumor_allele` |
| `Variant_Classification` | `variant_classification` |
| `Variant_Type` | `variant_type` |
| `HGVSc` | `hgvsc` |
| `HGVSp` | `hgvsp` |
| `t_depth` | `t_depth` |
| `t_alt_count` | `t_alt_count` |
| `data_type` | `data_type` |
| `experimental_strategy` | `experimental_strategy` |
| `access` | `access` |

### DGIdb Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `gene_name` | `gene_symbol` |
| `drug_name` | `drug_name` |
| `conceptId` | `concept_id` |
| `interactionTypes` | `interaction_types` |
| `interactionScore` | `interaction_score` |
| `approved` | `approved` |
| `immunotherapy` | `immunotherapy` |
| `antineoplastic` | `antineoplastic` |
| `geneCategories` | `gene_categories` |
| `sources` | `dgidb_sources` |
| `publications` | `publications` |

---

## TCGA Project Codes

| Code | Cancer Type |
|------|-------------|
| `TCGA-BRCA` | Breast Invasive Carcinoma |
| `TCGA-LUAD` | Lung Adenocarcinoma |
| `TCGA-LUSC` | Lung Squamous Cell Carcinoma |
| `TCGA-COAD` | Colon Adenocarcinoma |
| `TCGA-READ` | Rectum Adenocarcinoma |
| `TCGA-PRAD` | Prostate Adenocarcinoma |
| `TCGA-OV` | Ovarian Serous Cystadenocarcinoma |
| `TCGA-GBM` | Glioblastoma Multiforme |
| `TCGA-SKCM` | Skin Cutaneous Melanoma |
| `TCGA-KIRC` | Kidney Renal Clear Cell Carcinoma |

---

## DGIdb Gene Categories

| Category | Description |
|----------|-------------|
| `KINASE` | Protein kinase |
| `TUMOR SUPPRESSOR` | Tumor suppressor gene |
| `ONCOGENE` | Oncogene |
| `DRUGGABLE GENOME` | Part of the druggable genome |
| `DRUG RESISTANCE` | Associated with drug resistance |
| `TRANSCRIPTION FACTOR` | Transcription factor |
| `DNA REPAIR` | Involved in DNA repair |
| `CELL SURFACE` | Cell surface protein |

---

## Variant Classification Definitions

| Classification | Description |
|----------------|-------------|
| `Missense_Mutation` | Single nucleotide substitution causing amino acid change |
| `Nonsense_Mutation` | Substitution creating premature stop codon |
| `Frame_Shift_Del` | Deletion causing reading frame shift |
| `Frame_Shift_Ins` | Insertion causing reading frame shift |
| `Splice_Site` | Mutation affecting splice site |
| `Silent` | Synonymous substitution (no amino acid change) |

---

## Metadata Fields

| Field Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `_source.database` | string | Name of the source database | `GDC_TCGA`, `DGIdb` |
| `_source.version` | string | Version of the source data | `v36.0`, `4.2.0` |
| `_source.access_date` | date | Date the data was accessed | `2026-01-24` |
| `_source.original_id` | string | Original identifier in source | `TCGA-A1-A0SB-01A` |
