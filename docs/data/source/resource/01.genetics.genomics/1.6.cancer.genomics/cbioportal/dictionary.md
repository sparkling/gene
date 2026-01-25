# cBioPortal - Data Dictionary

## Overview

This data dictionary documents the schema for cBioPortal cancer genomics portal.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | cbioportal |
| **Name** | cBioPortal |
| **Parent** | 1.6.cancer.genomics |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Mutation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| entrezGeneId | integer | 1:1 | Yes | NCBI Gene ID | 7157 |
| gene | object | 1:1 | Yes | Gene object | {hugoGeneSymbol: "TP53"} |
| sampleId | string | 1:1 | Yes | Sample identifier | TCGA-AA-3518-01 |
| patientId | string | 1:1 | Yes | Patient identifier | TCGA-AA-3518 |
| studyId | string | 1:1 | Yes | Study identifier | coadread_tcga |
| mutationType | string | 1:1 | Yes | Mutation consequence | Missense_Mutation |
| proteinChange | string | 1:1 | No | Protein change | V600E |
| mutationStatus | string | 1:1 | No | Somatic/germline | Somatic |
| validationStatus | string | 1:1 | No | Validation status | Validated |
| chromosome | string | 1:1 | Yes | Chromosome | 7 |
| startPosition | integer | 1:1 | Yes | Start position | 140753336 |
| endPosition | integer | 1:1 | Yes | End position | 140753336 |
| referenceAllele | string | 1:1 | Yes | Reference allele | T |
| variantAllele | string | 1:1 | Yes | Variant allele | A |

### Copy Number Alteration

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| entrezGeneId | integer | 1:1 | Yes | NCBI Gene ID | 2064 |
| gene | object | 1:1 | Yes | Gene object | {hugoGeneSymbol: "ERBB2"} |
| sampleId | string | 1:1 | Yes | Sample identifier | TCGA-AA-3518-01 |
| alteration | integer | 1:1 | Yes | CNA value | 2 |
| cnaType | string | 1:1 | Yes | CNA category | AMP |

### Structural Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| site1EntrezGeneId | integer | 1:1 | Yes | Gene 1 ID | 238 |
| site2EntrezGeneId | integer | 1:1 | Yes | Gene 2 ID | 25 |
| site1HugoSymbol | string | 1:1 | Yes | Gene 1 symbol | ALK |
| site2HugoSymbol | string | 1:1 | Yes | Gene 2 symbol | EML4 |
| sampleId | string | 1:1 | Yes | Sample identifier | TCGA-AA-3518-01 |
| variantClass | string | 1:1 | Yes | Variant class | FUSION |

### Clinical Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| patientId | string | 1:1 | Yes | Patient identifier | TCGA-AA-3518 |
| sampleId | string | 1:1 | Yes | Sample identifier | TCGA-AA-3518-01 |
| studyId | string | 1:1 | Yes | Study identifier | coadread_tcga |
| clinicalAttributeId | string | 1:1 | Yes | Attribute name | OS_STATUS |
| value | string | 1:1 | Yes | Attribute value | DECEASED |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Study ID | lowercase_underscore | coadread_tcga | Study identifier |
| Patient ID | Study-specific | TCGA-AA-3518 | Patient identifier |
| Sample ID | Patient + suffix | TCGA-AA-3518-01 | Sample identifier |
| Entrez Gene ID | Integer | 7157 | NCBI Gene ID |
| Hugo Symbol | HGNC symbol | TP53 | Gene symbol |

---

## Enumerations

### Mutation Type

| Type | Description |
|------|-------------|
| Missense_Mutation | Amino acid substitution |
| Nonsense_Mutation | Stop codon gain |
| Frame_Shift_Del | Frameshift deletion |
| Frame_Shift_Ins | Frameshift insertion |
| In_Frame_Del | In-frame deletion |
| In_Frame_Ins | In-frame insertion |
| Splice_Site | Splice junction variant |
| Nonstop_Mutation | Stop codon loss |
| Translation_Start_Site | Start codon variant |
| Silent | Synonymous variant |

### CNA Alteration Values

| Value | Type | Meaning |
|-------|------|---------|
| -2 | HOMDEL | Homozygous deletion |
| -1 | HETLOSS | Heterozygous loss |
| 0 | DIPLOID | No change |
| 1 | GAIN | Low-level gain |
| 2 | AMP | High-level amplification |

### Mutation Status

| Status | Description |
|--------|-------------|
| Somatic | Somatic mutation |
| Germline | Germline variant |
| LOH | Loss of heterozygosity |
| None | Not determined |

### Validation Status

| Status | Description |
|--------|-------------|
| Validated | Experimentally validated |
| Untested | Not validated |
| Unknown | Validation unknown |

### Clinical Attributes

| Attribute | Description |
|-----------|-------------|
| OS_STATUS | Overall survival status |
| OS_MONTHS | Overall survival time |
| DFS_STATUS | Disease-free status |
| DFS_MONTHS | Disease-free time |
| AGE | Patient age |
| SEX | Patient sex |
| CANCER_TYPE | Cancer type |
| SAMPLE_TYPE | Tumor/normal |

---

## Entity Relationships

### Patient to Sample
- **Cardinality:** 1:N
- **Description:** One patient has multiple samples
- **Key Fields:** patientId, sampleId

### Sample to Mutation
- **Cardinality:** 1:N
- **Description:** One sample has multiple mutations
- **Key Fields:** sampleId, mutation

### Study to Sample
- **Cardinality:** 1:N
- **Description:** One study contains multiple samples
- **Key Fields:** studyId, sampleId

### Gene to Mutation
- **Cardinality:** 1:N
- **Description:** One gene has multiple mutations
- **Key Fields:** entrezGeneId, mutation

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| cBioPortal | Cancer Bioinformatics Portal | Platform name |
| TCGA | The Cancer Genome Atlas | Data source |
| CNA | Copy Number Alteration | Amplification/deletion |
| MAF | Mutation Annotation Format | File format |
| OS | Overall Survival | Clinical endpoint |
| DFS | Disease-Free Survival | Clinical endpoint |
| PFS | Progression-Free Survival | Clinical endpoint |
| TMB | Tumor Mutational Burden | Biomarker |
| MSI | Microsatellite Instability | Biomarker |
| GISTIC | Genomic Identification of Significant Targets | CNA method |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| TCGA | Sample ID | Data source |
| ICGC | Donor ID | Data source |
| dbSNP | rsID | Variant identifier |
| Ensembl | ENSG | Gene annotation |
| ClinVar | VCV | Clinical significance |
| OncoKB | Gene/Variant | Actionability |
| COSMIC | COSM | Mutation reference |

---

## API Endpoints

| Endpoint | Description |
|----------|-------------|
| /studies | List all studies |
| /genes | Gene information |
| /mutations | Mutation data |
| /molecular-profiles | Data profiles |
| /clinical-data | Clinical data |
| /copy-number-alterations | CNA data |
| /structural-variants | Fusion data |

---

## Data Quality Notes

1. **Cardinality:** One entry per sample-gene-variant observation
2. **Sources:** TCGA, ICGC, published studies, institution data
3. **Normalization:** Standard MAF format for mutations
4. **Public Studies:** 300+ public studies available
5. **API Access:** REST API freely available
6. **Visualization:** OncoPrint, mutation plots, survival analysis
