---
id: roadmap-epigenomics
title: "Roadmap Epigenomics Project"
type: data-source
category: pathways
subcategory: regulatory-networks
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [epigenomics, chromatin, histone-modifications, dna-methylation, regulatory]
---

# Roadmap Epigenomics Project

**Category:** [Pathways & Networks](../../_index.md) > [Regulatory Networks](../_index.md)

## Overview

The NIH Roadmap Epigenomics Mapping Consortium produced comprehensive epigenome reference maps for primary human cells and tissues. The project generated genome-wide maps of histone modifications, DNA methylation, chromatin accessibility, and RNA expression across 127 human cell types and tissues.

This resource provides the most comprehensive view of epigenomic regulation across diverse human cell types, enabling researchers to understand how epigenetic marks control gene expression in different biological contexts. The chromatin state annotations segment the genome into functional regions (promoters, enhancers, transcribed regions, etc.) based on combinatorial patterns of histone modifications.

Roadmap Epigenomics data is foundational for understanding regulatory variation, disease mechanisms, and tissue-specific gene regulation.

## Key Statistics

| Metric | Value |
|--------|-------|
| Reference Epigenomes | 127 |
| Cell/Tissue Types | 127 |
| Histone ChIP-seq Datasets | 2,800+ |
| DNA Methylation Datasets | 150+ |
| ATAC-seq/DNase-seq Datasets | 100+ |
| RNA-seq Datasets | 100+ |
| Total Data Volume | 150+ TB |

## Data Types

| Assay | Description | Resolution |
|-------|-------------|------------|
| Histone ChIP-seq | Histone modifications | 200 bp bins |
| WGBS | Whole genome bisulfite sequencing | Single CpG |
| RRBS | Reduced representation bisulfite | CpG islands |
| DNase-seq | Chromatin accessibility | 50 bp |
| ATAC-seq | Open chromatin | 50 bp |
| RNA-seq | Gene expression | Transcript-level |
| ChIA-PET | Chromatin interactions | Loop anchors |

## Histone Modifications

### Core Marks

| Mark | Association | Function |
|------|-------------|----------|
| H3K4me1 | Enhancers | Enhancer priming |
| H3K4me3 | Promoters | Active promoters |
| H3K27ac | Active enhancers | Active regulation |
| H3K27me3 | Polycomb repression | Gene silencing |
| H3K36me3 | Transcribed regions | Elongation |
| H3K9me3 | Heterochromatin | Constitutive silencing |

### Extended Marks

| Mark | Association |
|------|-------------|
| H3K9ac | Active chromatin |
| H3K4me2 | Promoters/enhancers |
| H2A.Z | Promoters |
| H3K79me2 | Transcription |

## Chromatin States (15-State Model)

| State | Color | Description | Marks |
|-------|-------|-------------|-------|
| 1_TssA | Red | Active TSS | H3K4me3, H3K27ac |
| 2_TssAFlnk | Orange | Flanking active TSS | H3K4me3, H3K4me1 |
| 3_TxFlnk | Yellow | Transcription at gene 5'/3' | H3K4me1, H3K36me3 |
| 4_Tx | Green | Strong transcription | H3K36me3 |
| 5_TxWk | Dark Green | Weak transcription | Weak H3K36me3 |
| 6_EnhG | Light Green | Genic enhancers | H3K4me1, H3K27ac |
| 7_Enh | Yellow | Enhancers | H3K4me1, H3K27ac |
| 8_ZNF/Rpts | Aqua | ZNF genes & repeats | H3K36me3, H3K9me3 |
| 9_Het | Light Blue | Heterochromatin | H3K9me3 |
| 10_TssBiv | Purple | Bivalent/poised TSS | H3K4me3, H3K27me3 |
| 11_BivFlnk | Dark Purple | Flanking bivalent TSS | H3K4me1, H3K27me3 |
| 12_EnhBiv | Pink | Bivalent enhancer | H3K4me1, H3K27me3 |
| 13_ReprPC | Gray | Repressed Polycomb | H3K27me3 |
| 14_ReprPCWk | Light Gray | Weak repressed Polycomb | Weak H3K27me3 |
| 15_Quies | White | Quiescent/low | No marks |

## Primary Use Cases

1. **GWAS variant interpretation** - Annotate variants with tissue-specific states
2. **Enhancer identification** - Find cell-type-specific enhancers
3. **Promoter activity** - Assess promoter activation across tissues
4. **Disease mechanism** - Link disease variants to regulatory elements
5. **Cell type characterization** - Define epigenomic signatures

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Epigenome ID | `E{###}` | E003 |
| Sample ID | ERS/EGS IDs | ERS192736 |
| BioSample | SAMN IDs | SAMN02632555 |
| Chromatin State | 1-15 | 7_Enh |

## Reference Epigenome Groups

| Group | Epigenomes | Description |
|-------|------------|-------------|
| ES/iPS | E001-E024 | Embryonic/induced pluripotent stem cells |
| Blood | E029-E050 | Hematopoietic cells |
| Brain | E053-E082 | Brain regions |
| Heart | E095-E105 | Cardiac tissues |
| GI Tract | E084-E094 | Gastrointestinal |
| Muscle | E107-E121 | Skeletal and smooth muscle |
| Other | E106, E122-E127 | Various tissues |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Data Portal | https://egg2.wustl.edu/roadmap/ | Browse/download |
| UCSC Browser | https://genome.ucsc.edu | Visualization |
| WashU EpiGenome Browser | https://epigenomegateway.wustl.edu | Interactive |
| AWS S3 | s3://encode-public/ | Bulk download |
| FTP | ftp://ftp.ncbi.nlm.nih.gov/pub/geo/ | GEO data |

### Download Examples

```bash
# Chromatin state files
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E003_15_coreMarks_mnemonics.bed.gz

# Signal tracks (bigWig)
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K4me3.pval.signal.bigwig

# DNA methylation
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E003_WGBS_FractionalMethylation.bigwig
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| BED | Chromatin states | Region annotation |
| bigWig | Signal tracks | Genome browser |
| BAM | Aligned reads | Re-analysis |
| bedGraph | Signal values | Processing |
| WIG | Wiggle format | Legacy tools |

## File Types

| Type | Description |
|------|-------------|
| chromhmmSegmentations | Chromatin state calls |
| signal/consolidated | Processed signal tracks |
| peaks/consolidated | Peak calls |
| dnamethylation | Methylation data |
| rnaseq | Expression data |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access |
| Commercial Use | Yes |
| Attribution | Required (cite consortium) |
| Data Sharing | GEO/SRA public |

## Cross-References

| Database | Relationship |
|----------|--------------|
| GEO | Raw data deposition |
| SRA | Sequencing reads |
| ENCODE | Complementary data |
| dbGaP | Protected data |
| Ensembl | Genome annotations |
| UCSC | Browser tracks |

## Related Projects

| Project | Description |
|---------|-------------|
| ENCODE | Encyclopedia of DNA Elements |
| Blueprint | European epigenome project |
| IHEC | International Human Epigenome Consortium |
| GTeX | Gene expression by tissue |

## Analysis Tools

| Tool | Description |
|------|-------------|
| ChromHMM | Chromatin state learning |
| GREAT | Region-gene association |
| LOLA | Region overlap analysis |
| deepTools | Signal processing |
| bedtools | Interval operations |

## Citation

Roadmap Epigenomics Consortium. (2015). Integrative analysis of 111 reference human epigenomes. Nature, 518(7539), 317-330.

## See Also

- [JASPAR](../jaspar/_index.md) - TF binding motifs
- [ENCODE](../../../01.genomic.references/1.4.functional.annotations/encode/_index.md) - Related functional data
- [Gene Ontology](../../4.5.gene.function.ontology/gene.ontology/_index.md) - Functional annotations
