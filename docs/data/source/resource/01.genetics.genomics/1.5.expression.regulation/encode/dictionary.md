# ENCODE - Data Dictionary

## Overview

This data dictionary documents the schema for ENCODE (Encyclopedia of DNA Elements) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | encode |
| **Name** | ENCODE |
| **Parent** | 1.5.expression.regulation |
| **Total Fields** | 40+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Experiment

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| accession | string | 1:1 | Yes | ENCODE accession | ENCSR000AAA |
| assay_title | string | 1:1 | Yes | Assay type | ChIP-seq |
| target | string | 1:1 | No | Target protein/histone | CTCF |
| biosample_ontology | string | 1:1 | Yes | Cell/tissue type | HepG2 |
| biosample_type | string | 1:1 | Yes | Sample category | cell line |
| status | string | 1:1 | Yes | Release status | released |
| date_released | date | 1:1 | No | Release date | 2020-01-15 |

### File

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| accession | string | 1:1 | Yes | File accession | ENCFF000AAA |
| file_format | string | 1:1 | Yes | File format | bigWig, bed |
| file_type | string | 1:1 | Yes | Data type | signal, peaks |
| output_type | string | 1:1 | Yes | Output category | fold change over control |
| assembly | string | 1:1 | Yes | Genome assembly | GRCh38 |
| file_size | integer | 1:1 | Yes | Size in bytes | 1048576 |
| md5sum | string | 1:1 | Yes | File checksum | a1b2c3d4e5... |

### Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| accession | string | 1:1 | Yes | Annotation accession | ENCSR000AAA |
| annotation_type | string | 1:1 | Yes | Annotation category | candidate cis-regulatory elements |
| collection_type | string | 1:1 | No | Collection | ENCODE4 |
| software_used | array | 1:N | No | Analysis software | MACS2, STAR |

### Biosample

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| accession | string | 1:1 | Yes | Biosample accession | ENCBS000AAA |
| biosample_ontology_id | string | 1:1 | Yes | Ontology ID | EFO:0002067 |
| organism | string | 1:1 | Yes | Species | Homo sapiens |
| sex | string | 1:1 | No | Biological sex | female |
| age | string | 1:1 | No | Age | 53 year |
| donor | string | 1:1 | No | Donor accession | ENCDO000AAA |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Experiment | ENCSR + 6 alphanumeric | ENCSR000AAA | Experiment accession |
| File | ENCFF + 6 alphanumeric | ENCFF000AAA | File accession |
| Biosample | ENCBS + 6 alphanumeric | ENCBS000AAA | Sample accession |
| Donor | ENCDO + 6 alphanumeric | ENCDO000AAA | Donor accession |
| Antibody | ENCAB + 6 alphanumeric | ENCAB000AAA | Antibody lot |
| cCRE | EH38E + 7 digits | EH38E1234567 | Regulatory element |

---

## Enumerations

### Assay Types

| Assay | Description | Target |
|-------|-------------|--------|
| ChIP-seq | Chromatin immunoprecipitation | Transcription factors, histones |
| ATAC-seq | Open chromatin | Accessible regions |
| RNA-seq | Transcriptome | Gene expression |
| DNase-seq | DNase hypersensitivity | Open chromatin |
| WGBS | Bisulfite sequencing | DNA methylation |
| Hi-C | Chromatin conformation | 3D structure |
| CAGE | Cap analysis | TSS mapping |
| RRBS | Reduced representation bisulfite | DNA methylation |

### Biosample Types

| Type | Description |
|------|-------------|
| cell line | Immortalized cell line |
| primary cell | Primary cells |
| tissue | Tissue sample |
| organoid | 3D organoid culture |
| in vitro differentiated cells | Differentiated cells |

### File Formats

| Format | Description | Use |
|--------|-------------|-----|
| bigWig | Binary signal track | Visualization |
| bigBed | Binary BED | Peaks/annotations |
| bed | Text BED | Peak calls |
| bam | Aligned reads | Alignments |
| fastq | Raw reads | Raw data |
| tsv | Tab-separated | Quantifications |

### Output Types

| Output | Description |
|--------|-------------|
| alignments | Mapped reads |
| signal | Coverage signal |
| peaks | Peak calls |
| IDR thresholded peaks | Reproducible peaks |
| fold change over control | Normalized signal |

### cCRE Classification

| Category | Description |
|----------|-------------|
| PLS | Promoter-like signature |
| pELS | Proximal enhancer-like |
| dELS | Distal enhancer-like |
| CTCF-only | CTCF-bound only |
| DNase-H3K4me3 | DNase + H3K4me3 |

---

## Entity Relationships

### Experiment to File
- **Cardinality:** 1:N
- **Description:** Experiments generate multiple files
- **Key Fields:** experiment_accession, file_accession

### Experiment to Biosample
- **Cardinality:** N:N
- **Description:** Experiments use multiple biosamples
- **Key Fields:** experiment_accession, biosample_accession

### Biosample to Donor
- **Cardinality:** N:1
- **Description:** Multiple samples from one donor
- **Key Fields:** biosample_accession, donor_accession

### File to Annotation
- **Cardinality:** N:N
- **Description:** Files contribute to annotations
- **Key Fields:** file_accession, annotation_accession

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ENCODE | Encyclopedia of DNA Elements | Project name |
| ChIP | Chromatin Immunoprecipitation | Assay type |
| ATAC | Assay for Transposase-Accessible Chromatin | Open chromatin |
| DNase | Deoxyribonuclease | Enzyme |
| WGBS | Whole Genome Bisulfite Sequencing | Methylation |
| cCRE | Candidate Cis-Regulatory Element | Registry element |
| IDR | Irreproducibility Discovery Rate | Reproducibility metric |
| TF | Transcription Factor | Protein type |
| TSS | Transcription Start Site | Gene feature |
| RPKM | Reads Per Kilobase Million | Expression unit |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| GEO | GSE/GSM | Data archive |
| NCBI SRA | SRR | Raw reads |
| Ensembl | ENSG | Gene annotation |
| EFO | EFO ID | Cell type ontology |
| UBERON | UBERON ID | Anatomy ontology |
| Roadmap | E ID | Epigenome reference |

---

## Data Quality Notes

1. **Cardinality:** One accession per experiment/file/sample
2. **Replicates:** Experiments require biological replicates
3. **Quality Metrics:** Extensive QC metrics per file
4. **Pipelines:** Standardized analysis pipelines
5. **Versioning:** Multiple phases (ENCODE2, ENCODE3, ENCODE4)
6. **Access:** All released data freely available
