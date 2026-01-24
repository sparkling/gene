---
id: schema-jaspar
title: "JASPAR Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: final
tags: [schema, database, transcription-factors, motifs]
---

# JASPAR - Schema Documentation

## TL;DR

JASPAR provides transcription factor binding profiles as position matrices (PFM/PWM) in multiple formats. The REST API returns JSON, while download files use JASPAR, TRANSFAC, or MEME formats.

## Matrix Profile Object (API)

```json
{
  "matrix_id": "MA0106.3",
  "name": "TP53",
  "version": 3,
  "collection": "CORE",
  "tax_group": "vertebrates",
  "species": [
    {
      "tax_id": "9606",
      "species": "Homo sapiens"
    }
  ],
  "class": "Zinc-coordinating",
  "family": "C2H2 zinc finger factors",
  "uniprot_ids": ["P04637"],
  "data_type": "ChIP-seq",
  "source": "ENCODE",
  "pubmed_ids": ["25378396"],
  "pfm": {
    "A": [286, 684, 660, 2, 704, 724, 70, 118, 18, 667, 295],
    "C": [148, 58, 56, 739, 61, 0, 601, 606, 127, 34, 262],
    "G": [283, 187, 78, 5, 3, 2, 31, 11, 15, 202, 183],
    "T": [137, 37, 32, 93, 70, 64, 84, 49, 664, 57, 115]
  },
  "sequence_logo": "http://jaspar.genereg.net/static/logos/all/MA0106.3.png"
}
```

## Position Frequency Matrix (PFM)

Raw nucleotide counts at each position:

```
>MA0106.3 TP53
A  [ 286  684  660    2  704  724   70  118   18  667  295 ]
C  [ 148   58   56  739   61    0  601  606  127   34  262 ]
G  [ 283  187   78    5    3    2   31   11   15  202  183 ]
T  [ 137   37   32   93   70   64   84   49  664   57  115 ]
```

### Matrix Interpretation

- Each column represents a position in the binding site
- Each row represents a nucleotide (A, C, G, T)
- Values are raw counts from aligned binding sites
- Higher count = more frequently observed

## Position Weight Matrix (PWM)

Log-odds scores relative to background:

```
>MA0106.3 TP53
A  [ 0.14  0.85  0.78 -4.30  0.90  0.96 -0.90 -0.50 -2.10  0.80  0.20 ]
C  [-0.35 -0.95 -1.00  1.15 -0.85 -4.50  0.75  0.76 -0.58 -1.20  0.08 ]
G  [ 0.12 -0.10 -0.65 -3.70 -3.50 -4.20 -1.40 -2.20 -2.50 -0.05 -0.15 ]
T  [-0.50 -1.70 -1.80 -0.55 -0.90 -1.00 -0.70 -1.40  0.88 -1.00 -0.45 ]
```

### PWM Calculation

```
PWM[i,j] = log2(PFM[i,j] / (N * background[i]))
```

Where:
- N = total count at position j
- background[i] = background frequency (often 0.25)

## Collections

| Collection | Description | Content |
|------------|-------------|---------|
| CORE | Non-redundant curated | High-quality binding profiles |
| CNE | Conserved non-coding elements | Regulatory elements |
| PHYLOFACTS | Phylogenetically inferred | Predicted from evolution |
| SPLICE | Splice factor motifs | Splicing regulation |
| POLII | RNA Pol II binding | Promoter elements |
| FAM | TF family averages | Family consensus motifs |

## TF Classifications

### Structural Classes

```json
{
  "classes": [
    {
      "name": "Zinc-coordinating",
      "families": [
        "C2H2 zinc finger factors",
        "Nuclear receptors with C4 zinc fingers",
        "C4 zinc finger-type integrases"
      ]
    },
    {
      "name": "Helix-turn-helix",
      "families": [
        "Fork head/winged helix factors",
        "Homeodomain factors",
        "POU domain factors"
      ]
    },
    {
      "name": "Basic domains",
      "families": [
        "bHLH factors",
        "bZIP factors"
      ]
    }
  ]
}
```

## MEME Format

```
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0106.3 TP53
letter-probability matrix: alnum= 4 w= 11 nsites= 854 E= 0
 0.335  0.173  0.331  0.160
 0.801  0.068  0.219  0.043
 0.773  0.066  0.091  0.037
 0.002  0.865  0.006  0.109
 0.824  0.071  0.004  0.082
 0.848  0.000  0.002  0.075
 0.082  0.704  0.036  0.098
 0.138  0.710  0.013  0.057
 0.021  0.149  0.018  0.778
 0.781  0.040  0.237  0.067
 0.345  0.307  0.214  0.135
```

## TRANSFAC Format

```
AC  MA0106.3
XX
ID  TP53
XX
DE  Homo sapiens TP53
XX
PO      A      C      G      T
01    286    148    283    137
02    684     58    187     37
03    660     56     78     32
04      2    739      5     93
05    704     61      3     70
06    724      0      2     64
07     70    601     31     84
08    118    606     11     49
09     18    127     15    664
10    667     34    202     57
11    295    262    183    115
XX
CC  program: JASPAR
CC  collection: CORE
CC  tax_group: vertebrates
XX
//
```

## Binding Sites

```json
{
  "sites": [
    {
      "matrix_id": "MA0106.3",
      "seq_id": "ENSG00000141510_1",
      "start": 100,
      "end": 110,
      "strand": "+",
      "sequence": "AGACATGTCC",
      "score": 0.92,
      "source": "ChIP-seq"
    }
  ]
}
```

## Information Content

```json
{
  "matrix_id": "MA0106.3",
  "ic_values": [0.15, 0.89, 0.82, 1.75, 0.96, 1.12, 0.78, 0.82, 1.45, 0.72, 0.28],
  "total_ic": 9.74,
  "max_possible_ic": 22.0
}
```

- Information Content (IC) measured in bits
- Maximum IC per position = 2 bits
- Higher IC = more conserved position

## See Also

- [Download Documentation](./download.md)
- [JASPAR Matrix Format](https://jaspar.genereg.net/docs/)
