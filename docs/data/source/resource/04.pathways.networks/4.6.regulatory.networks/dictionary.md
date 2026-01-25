# 4.6 Regulatory Networks - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| Subcategory ID | 4.6 |
| Subcategory Name | Regulatory Networks |
| Data Sources | JASPAR, Roadmap Epigenomics |
| Schema ID | `https://gene.ai/schemas/4.6-regulatory-networks.json` |

## Unified Fields

These fields are harmonized across all data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `element_id` | string | Required (1:1) | Primary identifier for regulatory element | JASPAR, Roadmap Epigenomics | JASPAR: `MA0106.3`, Roadmap: `E003` |
| `name` | string | Optional (1:1) | Name of the element or sample | JASPAR, Roadmap Epigenomics | JASPAR: `TP53`, Roadmap: `H1 Cells` |
| `species` | string | Optional (1:N) | Organism(s) | JASPAR, Roadmap Epigenomics | `Homo sapiens` |

## JASPAR-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `matrix_id` | string | Required (1:1) | JASPAR matrix identifier with version | `MA0106.3` |
| `version` | integer | Optional (1:1) | Matrix version number | `3` |
| `collection` | enum | Optional (1:1) | JASPAR collection | `CORE`, `CNE`, `PHYLOFACTS`, `SPLICE`, `POLII`, `FAM` |
| `tax_group` | string | Optional (1:1) | Taxonomic group | `vertebrates` |
| `class` | string | Optional (1:1) | TF structural class | `Zinc-coordinating` |
| `family` | string | Optional (1:1) | TF family | `C2H2 zinc finger factors` |
| `uniprot_ids` | array&lt;string&gt; | Optional (1:N) | UniProt accessions for the TF | `["P04637"]` |
| `data_type` | string | Optional (1:1) | Experimental data source type | `ChIP-seq` |
| `source` | string | Optional (1:1) | Data source project | `ENCODE` |
| `pfm` | object | Optional (1:1) | Position Frequency Matrix (nucleotide counts) | `{"A": [286, 684], "C": [148, 58], "G": [283, 187], "T": [137, 37]}` |
| `pwm` | object | Optional (1:1) | Position Weight Matrix (log-odds scores) | `{"A": [0.14, 0.85], "C": [-0.35, -0.95], "G": [0.12, -0.10], "T": [-0.50, -1.70]}` |
| `sequence_logo` | string | Optional (1:1) | URL to sequence logo image | `http://jaspar.genereg.net/static/logos/all/MA0106.3.png` |
| `ic_values` | array&lt;float&gt; | Optional (1:N) | Information content per position (bits) | `[0.15, 0.89, 0.82]` |
| `total_ic` | float | Optional (1:1) | Total information content | `9.74` |

## Roadmap Epigenomics-Specific Fields

### Sample/Epigenome Metadata

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `eid` | string | Required (1:1) | Epigenome ID | `E003` |
| `standardName` | string | Required (1:1) | Standard epigenome name | `H1 Cells` |
| `group` | string | Optional (1:1) | Epigenome group | `ES/iPS` |
| `anatomy` | string | Optional (1:1) | Anatomical source | `ESC` |
| `type` | string | Optional (1:1) | Sample type | `Cell Line` |

### Chromatin State Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `chromatin_state` | string | Optional (1:1) | Chromatin state annotation (15-state model) | `1_TssA` |
| `state` | enum | Optional (1:1) | 15-state chromatin state code | See Chromatin States table below |
| `histone_marks` | object | Optional (1:1) | Available histone modification tracks | `{"H3K4me1": "available", "H3K4me3": "available", "H3K27ac": "available"}` |

### Genomic Coordinates (BED format)

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `chrom` | string | Optional (1:1) | Chromosome (BED format) | `chr1` |
| `chromStart` | integer | Optional (1:1) | Start position (0-based, BED format) | `10000` |
| `chromEnd` | integer | Optional (1:1) | End position (BED format) | `10600` |

### DNA Methylation Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `dna_methylation` | enum | Optional (1:1) | DNA methylation assay type | `WGBS`, `RRBS` |
| `methylation_fraction` | float | Optional (1:1) | DNA methylation level (0-1) | `0.83` |

### Peak Calling Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `signal_type` | enum | Optional (1:1) | Signal track type | `pval`, `fc`, `tagAlign` |
| `peak_score` | integer | Optional (1:1) | Peak score (0-1000, narrowPeak format) | `1000` |
| `signalValue` | float | Optional (1:1) | Enrichment signal value | `9.12` |
| `pValue` | float | Optional (1:1) | -log10(p-value) for peak | `16.23` |
| `qValue` | float | Optional (1:1) | -log10(q-value/FDR) for peak | `14.56` |
| `peak` | integer | Optional (1:1) | Summit position relative to start | `426` |

## 15-State Chromatin Model

| State | Name | Color | Key Marks |
|-------|------|-------|-----------|
| 1_TssA | Active TSS | Red | H3K4me3, H3K27ac |
| 2_TssAFlnk | Flanking Active TSS | Orange Red | H3K4me3, H3K4me1 |
| 3_TxFlnk | Transcription 5'/3' | Lime Green | H3K4me1, H3K36me3 |
| 4_Tx | Strong Transcription | Green | H3K36me3 |
| 5_TxWk | Weak Transcription | Dark Green | Weak H3K36me3 |
| 6_EnhG | Genic Enhancers | Green Yellow | H3K4me1, H3K27ac, H3K36me3 |
| 7_Enh | Enhancers | Yellow | H3K4me1, H3K27ac |
| 8_ZNF/Rpts | ZNF genes & repeats | Medium Aquamarine | H3K36me3, H3K9me3 |
| 9_Het | Heterochromatin | Pale Turquoise | H3K9me3 |
| 10_TssBiv | Bivalent/Poised TSS | Indian Red | H3K4me3, H3K27me3 |
| 11_BivFlnk | Flanking Bivalent TSS/Enh | Dark Salmon | H3K4me1, H3K27me3 |
| 12_EnhBiv | Bivalent Enhancer | Dark Khaki | H3K4me1, H3K27me3 |
| 13_ReprPC | Repressed PolyComb | Silver | H3K27me3 |
| 14_ReprPCWk | Weak Repressed PolyComb | Gainsboro | Weak H3K27me3 |
| 15_Quies | Quiescent/Low | White | No marks |

## JASPAR TF Classes

| Class | Families |
|-------|----------|
| Zinc-coordinating | C2H2 zinc finger factors, Nuclear receptors with C4 zinc fingers, C4 zinc finger-type integrases |
| Helix-turn-helix | Fork head/winged helix factors, Homeodomain factors, POU domain factors |
| Basic domains | bHLH factors, bZIP factors |

## Source Metadata

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `_source.database` | string | Required | Source database name | `JASPAR`, `Roadmap Epigenomics` |
| `_source.version` | string | Optional | Database version | `2024`, `Release 9` |
| `_source.access_date` | date | Optional | Date data was retrieved | `2026-01-24` |
| `_source.original_id` | string | Optional | Original identifier in source | `MA0106.3` |

## Field Mappings by Source

### JASPAR Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `matrix_id` | `matrix_id` |
| `name` | `name` (also `jaspar_tf_name`) |
| `collection` | `collection` |
| `tax_group` | `tax_group` |
| `class` | `class` |
| `family` | `family` |
| `uniprot_ids` | `uniprot_ids` |
| `data_type` | `data_type` |
| `source` | `source` |
| `pfm` | `pfm` |
| `pwm` | `pwm` |
| `sequence_logo` | `sequence_logo` |
| `ic_values` | `ic_values` |
| `total_ic` | `total_ic` |

### Roadmap Epigenomics Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `EID` | `eid` |
| `STANDARDIZED_NAME` | `standardName` |
| `GROUP` | `group` |
| `ANATOMY` | `anatomy` |
| `TYPE` | `type` |
| `chrom` | `chrom` |
| `chromStart` | `chromStart` |
| `chromEnd` | `chromEnd` |
| `state` | `state` |
| `signal_type` | `signal_type` |
| `methylation_fraction` | `methylation_fraction` |
| `peak_score` | `peak_score` |
| `signalValue` | `signalValue` |
| `pValue` | `pValue` |
| `qValue` | `qValue` |
| `peak` | `peak` |

## Data File Formats

### JASPAR Matrix Formats

| Format | Description | Extension |
|--------|-------------|-----------|
| JASPAR | Native JASPAR format | `.jaspar` |
| TRANSFAC | TRANSFAC format | `.transfac` |
| MEME | MEME suite format | `.meme` |
| PFM | Position frequency matrix | `.pfm` |
| PWM | Position weight matrix | `.pwm` |

### Roadmap Epigenomics Formats

| Format | Description | Extension |
|--------|-------------|-----------|
| BED | Browser Extensible Data | `.bed` |
| narrowPeak | ENCODE narrowPeak | `.narrowPeak` |
| broadPeak | ENCODE broadPeak | `.broadPeak` |
| bigWig | Compressed signal | `.bigWig` |
| bigBed | Compressed BED | `.bigBed` |
