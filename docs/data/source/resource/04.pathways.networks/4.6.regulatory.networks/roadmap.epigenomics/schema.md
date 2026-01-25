---
id: schema-roadmap-epigenomics
title: "Roadmap Epigenomics Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, database, epigenomics, chromatin]
---

# Roadmap Epigenomics - Schema Documentation

## TL;DR

The Roadmap Epigenomics project provides genome-wide epigenomic maps including histone modifications, DNA methylation, and chromatin accessibility across 127 human cell types. Data is stored in standard genomic formats (BED, bigWig, BAM).

## Chromatin State BED Format

### 15-State Model

```
chr1	10000	10600	15_Quies
chr1	10600	11137	14_ReprPCWk
chr1	11137	11737	13_ReprPC
chr1	11737	12137	7_Enh
chr1	12137	14537	1_TssA
```

### BED Columns

| Column | Description | Example |
|--------|-------------|---------|
| chrom | Chromosome | chr1 |
| chromStart | Start position (0-based) | 10000 |
| chromEnd | End position | 10600 |
| state | Chromatin state | 15_Quies |

## Chromatin States (15-State Model)

| State | Name | Color | Description | Key Marks |
|-------|------|-------|-------------|-----------|
| 1 | TssA | Red | Active TSS | H3K4me3, H3K27ac |
| 2 | TssAFlnk | Orange Red | Flanking Active TSS | H3K4me3, H3K4me1 |
| 3 | TxFlnk | Lime Green | Transcription 5'/3' | H3K4me1, H3K36me3 |
| 4 | Tx | Green | Strong Transcription | H3K36me3 |
| 5 | TxWk | Dark Green | Weak Transcription | Weak H3K36me3 |
| 6 | EnhG | Green Yellow | Genic Enhancers | H3K4me1, H3K27ac, H3K36me3 |
| 7 | Enh | Yellow | Enhancers | H3K4me1, H3K27ac |
| 8 | ZNF/Rpts | Medium Aquamarine | ZNF genes & repeats | H3K36me3, H3K9me3 |
| 9 | Het | Pale Turquoise | Heterochromatin | H3K9me3 |
| 10 | TssBiv | Indian Red | Bivalent/Poised TSS | H3K4me3, H3K27me3 |
| 11 | BivFlnk | Dark Salmon | Flanking Bivalent TSS/Enh | H3K4me1, H3K27me3 |
| 12 | EnhBiv | Dark Khaki | Bivalent Enhancer | H3K4me1, H3K27me3 |
| 13 | ReprPC | Silver | Repressed PolyComb | H3K27me3 |
| 14 | ReprPCWk | Gainsboro | Weak Repressed PolyComb | Weak H3K27me3 |
| 15 | Quies | White | Quiescent/Low | No marks |

## Epigenome Metadata

```json
{
  "eid": "E003",
  "standardName": "H1 Cells",
  "group": "ES/iPS",
  "anatomy": "ESC",
  "type": "Cell Line",
  "age": "Unknown",
  "sex": "Male",
  "description": "H1 human embryonic stem cells",
  "histone_marks": {
    "H3K4me1": "available",
    "H3K4me3": "available",
    "H3K27ac": "available",
    "H3K27me3": "available",
    "H3K36me3": "available",
    "H3K9me3": "available"
  },
  "dna_methylation": "WGBS",
  "rna_seq": "available"
}
```

## Epigenome Groups

| Group | EIDs | Description |
|-------|------|-------------|
| ESC | E001-E003 | Embryonic Stem Cells |
| iPSC | E004-E024 | Induced Pluripotent Stem Cells |
| ES-deriv | E004-E020 | ES-derived cells |
| Blood & T-cell | E029-E050 | Hematopoietic |
| HSC & B-cell | E029-E046 | Stem cells, B cells |
| Mesench | E049-E058 | Mesenchymal |
| Epithelial | E057-E061 | Epithelial |
| Myosat | E052-E054 | Muscle Satellite |
| Brain | E053-E082 | Brain regions |
| Adipose | E063 | Fat tissue |
| Heart | E095-E105 | Cardiac |
| Sm. Muscle | E076-E078 | Smooth Muscle |
| Digestive | E084-E094 | GI Tract |
| Thymus | E112-E113 | Thymus |
| Other | E106-E127 | Miscellaneous |

## Signal Track Format (bigWig)

### File Naming Convention

```
{EID}-{mark}.{type}.{analysis}.bigwig

Example: E003-H3K4me3.pval.signal.bigwig
```

### Signal Types

| Type | Description |
|------|-------------|
| pval | -log10(p-value) |
| fc | Fold change over control |
| tagAlign | Tag density |

## DNA Methylation Format

### WGBS bedGraph

```
chr1	10468	10470	0.83
chr1	10470	10472	0.91
chr1	10483	10485	0.87
```

Columns: chrom, start, end, methylation_fraction (0-1)

### RRBS Coverage

```
chr1	10468	10470	0.83	12	10
```

Columns: chrom, start, end, methylation_fraction, total_reads, methylated_reads

## Peak Files (BED narrowPeak)

```
chr1	713875	714727	E003_H3K4me3_narrowPeak.1	1000	.	9.12	16.23	14.56	426
```

### narrowPeak Columns

| Column | Description |
|--------|-------------|
| 1-3 | chrom, start, end |
| 4 | name |
| 5 | score (0-1000) |
| 6 | strand |
| 7 | signalValue |
| 8 | pValue (-log10) |
| 9 | qValue (-log10) |
| 10 | peak (summit position relative to start) |

## RNA-seq Expression

### RPKM Format

```
gene_id	E003	E004	E005
ENSG00000141510	15.32	12.45	18.76
ENSG00000012048	8.45	9.12	7.89
```

### TPM Format

Similar structure with TPM values.

## File Types Summary

| Data Type | Format | Extension |
|-----------|--------|-----------|
| Chromatin States | BED | .bed.gz |
| Signal Tracks | bigWig | .bigwig |
| Peak Calls | narrowPeak | .narrowPeak.gz |
| DNA Methylation | bedGraph | .bedGraph.gz |
| Alignments | BAM | .bam |
| Expression | TSV | .RPKM.pc.gz |

## See Also

- [Download Documentation](./download.md)
- [Roadmap Data Portal](https://egg2.wustl.edu/roadmap/)
