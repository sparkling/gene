---
id: jaspar
title: "JASPAR - Transcription Factor Binding Profiles"
type: source
parent: ../README.md
tier: 1
status: active
category: pathways.networks
subcategory: regulatory.networks
tags:
  - transcription-factors
  - binding-motifs
  - pwm
  - regulatory
  - open-access
---

# JASPAR - Transcription Factor Binding Profiles

**Category:** [Pathways & Networks](../../README.md) > [Regulatory Networks](../README.md)

## Overview

JASPAR is the largest open-access database of curated, non-redundant transcription factor (TF) binding profiles. It provides position weight matrices (PWMs) and position frequency matrices (PFMs) derived from experimentally defined binding sites, enabling computational prediction of TF binding sites in DNA sequences.

The database covers transcription factors from multiple taxonomic groups including vertebrates, plants, insects, nematodes, fungi, and urochordates. JASPAR profiles are generated from experimentally verified binding sites obtained through techniques like SELEX, ChIP-seq, and protein binding microarrays.

JASPAR is the gold standard for TF binding motif analysis and is completely open-access under CC BY 4.0, making it the primary choice for regulatory network analysis.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Profiles | 2,000+ |
| CORE Collection | 800+ |
| Vertebrate Profiles | 750+ |
| Plant Profiles | 500+ |
| TF Species | 800+ |
| Taxonomic Groups | 6 |

## Profile Collections

| Collection | Description | Profiles |
|------------|-------------|----------|
| CORE | Non-redundant, curated | 800+ |
| CNE | Conserved non-coding elements | 200+ |
| PHYLOFACTS | Phylogenetically inferred | 150+ |
| SPLICE | Splice factor motifs | 50+ |
| POLII | RNA Pol II binding | 50+ |
| FAM | TF family motifs | 100+ |

## Primary Use Cases

1. **Binding site prediction** - Scan sequences for TF binding sites
2. **Promoter analysis** - Identify regulatory elements in promoters
3. **ChIP-seq analysis** - Motif enrichment in peak regions
4. **Regulatory network inference** - Build TF-target networks
5. **Comparative genomics** - Cross-species motif conservation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| JASPAR ID | `MA{#####}.{#}` | MA0106.3 |
| Matrix Name | Text | TP53 |
| UniProt | Protein accession | P04637 |
| TF Class | Classification | Zinc finger |

## Matrix Formats

### Position Frequency Matrix (PFM)

Raw nucleotide counts at each position:

```
>MA0106.3 TP53
A  [ 286  684  660    2  704  724   70  118   18  667  295  ]
C  [ 148   58   56  739   61    0  601  606  127   34  262  ]
G  [ 283  187   78    5    3    2   31   11   15  202  183  ]
T  [ 137   37   32   93   70   64   84   49  664   57  115  ]
```

### Position Weight Matrix (PWM)

Log-odds scores:

```
>MA0106.3 TP53
A  [ 0.14  0.85  0.78 -4.30  0.90  0.96 -0.90 -0.50 -2.10  0.80  0.20 ]
C  [-0.35 -0.95 -1.00  1.15 -0.85 -4.50  0.75  0.76 -0.58 -1.20  0.08 ]
G  [ 0.12 -0.10 -0.65 -3.70 -3.50 -4.20 -1.40 -2.20 -2.50 -0.05 -0.15 ]
T  [-0.50 -1.70 -1.80 -0.55 -0.90 -1.00 -0.70 -1.40  0.88 -1.00 -0.45 ]
```

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://jaspar.genereg.net | Browse/search |
| REST API | https://jaspar.genereg.net/api/ | JSON responses |
| Downloads | https://jaspar.genereg.net/download/ | Matrix files |
| MEME Suite | https://meme-suite.org | Analysis tools |

### API Examples

```bash
# Get matrix by ID
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.3/"

# Search by TF name
curl "https://jaspar.genereg.net/api/v1/matrix/?name=TP53"

# Get all vertebrate matrices
curl "https://jaspar.genereg.net/api/v1/matrix/?tax_group=vertebrates"

# Download PFM format
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.3/?format=pfm"

# Search by TF class
curl "https://jaspar.genereg.net/api/v1/matrix/?class=Zinc-finger"

# Get binding sites
curl "https://jaspar.genereg.net/api/v1/sites/MA0106.3/"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| PFM | Position Frequency Matrix | Raw counts |
| PWM | Position Weight Matrix | Log-odds scores |
| MEME | MEME motif format | MEME Suite tools |
| TRANSFAC | TRANSFAC format | Legacy tools |
| JSON | API responses | Programmatic access |

## TF Classification

### Structural Classes

| Class | Description | Examples |
|-------|-------------|----------|
| Zinc finger | Zinc-coordinating domains | SP1, TP53, KLF4 |
| Helix-turn-helix | HTH DNA binding | FOXA1, FOXP2 |
| Basic helix-loop-helix | bHLH dimers | MYC, MAX |
| Basic leucine zipper | bZIP dimers | JUN, FOS |
| Homeodomain | 60 aa helix-turn-helix | HOX genes |
| Nuclear receptor | Steroid receptors | ESR1, AR |

### Taxonomic Groups

| Group | Organisms | Profiles |
|-------|-----------|----------|
| Vertebrates | Human, mouse, etc. | 750+ |
| Plants | Arabidopsis, rice, etc. | 500+ |
| Insects | Drosophila, etc. | 150+ |
| Nematodes | C. elegans | 50+ |
| Fungi | Yeast, etc. | 200+ |
| Urochordates | Ciona | 50+ |

## Matrix Quality Metrics

| Metric | Description |
|--------|-------------|
| Information Content | Bits of information per position |
| Motif Length | Number of positions |
| Site Count | Number of aligned sites |
| Validation | Experimental validation level |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required |
| Redistribution | Allowed |

## Cross-References

| Database | Relationship |
|----------|--------------|
| UniProt | TF protein IDs |
| Ensembl | Gene identifiers |
| HGNC | Gene symbols |
| TFClass | TF classification |
| ENCODE | ChIP-seq validation |
| HOCOMOCO | Alternative motif DB |

## Scanning Tools

| Tool | Description |
|------|-------------|
| FIMO | Find Individual Motif Occurrences |
| MAST | Motif Alignment & Search Tool |
| MEME | Motif discovery |
| TFBSTools | R/Bioconductor package |
| RSAT | Regulatory Sequence Analysis |

## Limitations

- Motif scanning produces many false positive predictions
- Position weight matrices assume position independence
- Coverage varies across taxonomic groups
- In silico predictions require experimental validation

## Download Files

| File | Description |
|------|-------------|
| JASPAR2024_CORE_vertebrates_non-redundant_pfms.txt | Vertebrate PFMs |
| JASPAR2024_CORE_plants_non-redundant_pfms.txt | Plant PFMs |
| JASPAR2024_CORE_all_pfms.txt | All PFMs |
| JASPAR2024_CORE_all_pfms.meme | MEME format |

## See Also

- [Roadmap Epigenomics](../roadmap.epigenomics/README.md) - Chromatin states
- [ENCODE](../../../01.genomic.references/1.4.functional.annotations/encode/README.md) - TF ChIP-seq data
- [Gene Ontology](../../4.5.gene.function.ontology/gene.ontology/README.md) - TF annotations
