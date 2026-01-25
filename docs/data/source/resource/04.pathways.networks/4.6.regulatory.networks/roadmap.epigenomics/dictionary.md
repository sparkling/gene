# Roadmap Epigenomics - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | roadmap.epigenomics |
| **Name** | Roadmap Epigenomics |
| **Total Fields** | 30 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| eid | String | Yes | Epigenome identifier | E003 |
| standardName | String | Yes | Standard epigenome name | H1 Cells |
| group | Enum | No | Epigenome group classification | ESC |
| anatomy | String | No | Anatomical source | embryonic stem cells |
| type | String | No | Sample type | Cell Line |
| age | String | No | Donor age | fetal |
| sex | Enum | No | Donor sex | Male, Female |

---

## Histone Modification Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| histone_marks.H3K4me1 | Enum | Enhancer mark availability |
| histone_marks.H3K4me3 | Enum | Promoter mark availability |
| histone_marks.H3K27ac | Enum | Active enhancer/promoter mark |
| histone_marks.H3K27me3 | Enum | Polycomb repression mark |
| histone_marks.H3K36me3 | Enum | Transcription elongation mark |
| histone_marks.H3K9me3 | Enum | Heterochromatin mark |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Epigenome ID | E### | Roadmap sample ID | E003 |
| Chromosome | chr## | Chromosome name | chr1 |
| Coordinates | Numeric | Genomic positions | 12345-67890 |

---

## Enumerations

### Epigenome Groups

| Value | Description |
|-------|-------------|
| ESC | Embryonic stem cells |
| iPSC | Induced pluripotent stem cells |
| ES-deriv | ES-derived cells |
| Blood & T-cell | Blood and T cells |
| HSC & B-cell | Hematopoietic stem cells and B cells |
| Mesench | Mesenchymal cells |
| Epithelial | Epithelial cells |
| Myosat | Muscle satellite cells |
| Brain | Brain tissues |
| Adipose | Fat tissue |
| Heart | Heart tissue |
| Sm. Muscle | Smooth muscle |
| Digestive | Digestive system |
| Thymus | Thymus |
| Other | Other cell types |

### Sample Types

| Value | Description |
|-------|-------------|
| Cell Line | Cultured cell line |
| Primary Cell | Primary isolated cells |
| Primary Tissue | Primary tissue sample |

### 15-State Chromatin Model

| State | Name | Description |
|-------|------|-------------|
| 1_TssA | Active TSS | Active transcription start site |
| 2_TssAFlnk | Flanking Active TSS | Flanking active promoter |
| 3_TxFlnk | Transcription Flanking | Transcription at gene 5'/3' |
| 4_Tx | Strong Transcription | Active transcription |
| 5_TxWk | Weak Transcription | Weak transcription |
| 6_EnhG | Genic Enhancers | Enhancers within genes |
| 7_Enh | Enhancers | Active enhancers |
| 8_ZNF/Rpts | ZNF/Repeats | Zinc finger genes and repeats |
| 9_Het | Heterochromatin | Constitutive heterochromatin |
| 10_TssBiv | Bivalent TSS | Bivalent/poised promoter |
| 11_BivFlnk | Flanking Bivalent | Flanking bivalent TSS |
| 12_EnhBiv | Bivalent Enhancer | Bivalent enhancer |
| 13_ReprPC | Repressed Polycomb | Polycomb-repressed |
| 14_ReprPCWk | Weak Repressed Polycomb | Weak repression |
| 15_Quies | Quiescent | Quiescent/low signal |

### DNA Methylation Types

| Value | Description |
|-------|-------------|
| WGBS | Whole-genome bisulfite sequencing |
| RRBS | Reduced representation bisulfite sequencing |
| not_available | No methylation data |

---

## File Formats

| Format | Description |
|--------|-------------|
| bigWig | Signal track format |
| narrowPeak | Peak calls (BED6+4) |
| broadPeak | Broad peak calls |
| BED | Chromatin states |
| bedGraph | Methylation values |

---

## Entity Relationships

### Epigenome to Marks
- **Cardinality:** 1:N
- **Description:** Each epigenome has multiple data types
- **Key Fields:** eid, histone_marks

### Region to State
- **Cardinality:** 1:1
- **Description:** Genomic regions have chromatin state
- **Key Fields:** chrom, chromStart, chromEnd, state

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| EID | Epigenome Identifier | Sample ID format |
| TSS | Transcription Start Site | Promoter region |
| ChIP-seq | Chromatin Immunoprecipitation Sequencing | Method |
| WGBS | Whole-Genome Bisulfite Sequencing | Methylation method |
| RRBS | Reduced Representation Bisulfite Sequencing | Methylation method |
| H3K4me1 | Histone H3 Lysine 4 Monomethylation | Enhancer mark |
| H3K4me3 | Histone H3 Lysine 4 Trimethylation | Promoter mark |
| H3K27ac | Histone H3 Lysine 27 Acetylation | Active mark |
| H3K27me3 | Histone H3 Lysine 27 Trimethylation | Repressive mark |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
