# CADD - Data Dictionary

## Overview

This data dictionary documents the schema for CADD (Combined Annotation Dependent Depletion) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | cadd |
| **Name** | CADD |
| **Parent** | 1.2.functional.prediction |
| **Total Fields** | 15+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Chrom | string | 1:1 | Yes | Chromosome | 1, X, MT |
| Pos | integer | 1:1 | Yes | Genomic position (1-based) | 11856378 |
| Ref | string | 1:1 | Yes | Reference allele | A |
| Alt | string | 1:1 | Yes | Alternate allele | G |

### CADD Scores

| Field Name | Data Type | Cardinality | Required | Description | Range |
|------------|-----------|-------------|----------|-------------|-------|
| RawScore | float | 1:1 | Yes | Raw CADD score | -7 to 35 |
| PHRED | float | 1:1 | Yes | Phred-scaled score | 0-99 |

### Annotation Features

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Type | string | 1:1 | No | Variant type | SNV, INS, DEL |
| Length | integer | 1:1 | No | Variant length | 1, 5, 100 |
| AnnoType | string | 1:1 | No | Annotation type | CodingTranscript, Intergenic |
| Consequence | string | 1:1 | No | VEP consequence | missense_variant |
| ConsScore | integer | 1:1 | No | Consequence severity | 0-10 |
| ConsDetail | string | 1:1 | No | Detailed consequence | stop_gained |
| GeneID | string | 1:1 | No | Ensembl gene ID | ENSG00000141510 |
| GeneName | string | 1:1 | No | Gene symbol | TP53 |
| FeatureID | string | 1:1 | No | Transcript ID | ENST00000269305 |

### Conservation Annotations

| Field Name | Data Type | Cardinality | Required | Description | Range |
|------------|-----------|-------------|----------|-------------|-------|
| GerpRS | float | 1:1 | No | GERP rejected substitutions | -12 to 6 |
| GerpN | float | 1:1 | No | GERP neutral rate | 0-10 |
| GerpS | float | 1:1 | No | GERP S score | -12 to 6 |
| phyloP | float | 1:1 | No | phyloP conservation | -20 to 10 |
| phastCons | float | 1:1 | No | phastCons probability | 0-1 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Variant Key | chr-pos-ref-alt | 17-7577121-G-A | Genomic coordinates |
| Ensembl Gene | ENSG + 11 digits | ENSG00000141510 | Gene ID |
| Ensembl Transcript | ENST + 11 digits | ENST00000269305 | Transcript ID |

---

## Enumerations

### PHRED Score Interpretation

| Score Range | Percentile | Interpretation |
|-------------|------------|----------------|
| >= 30 | Top 0.1% | Very high deleteriousness |
| >= 20 | Top 1% | High deleteriousness |
| >= 15 | Top 3% | Moderate-high deleteriousness |
| >= 10 | Top 10% | Moderate deleteriousness |
| < 10 | Bottom 90% | Lower deleteriousness |

### Consequence Types

| Value | Meaning | ConsScore |
|-------|---------|-----------|
| synonymous_variant | Silent change | 1 |
| missense_variant | Amino acid change | 5 |
| stop_gained | Premature stop codon | 8 |
| frameshift_variant | Reading frame shift | 9 |
| splice_donor_variant | Splice site disruption | 8 |
| splice_acceptor_variant | Splice site disruption | 8 |

### Annotation Types

| Value | Meaning |
|-------|---------|
| CodingTranscript | In protein-coding gene |
| NonCodingTranscript | In non-coding RNA |
| RegulatoryFeature | In regulatory region |
| Intergenic | Between genes |

---

## Entity Relationships

### Variant to Gene
- **Cardinality:** N:M
- **Description:** Variants can affect multiple genes; genes have many variants
- **Key Fields:** GeneID, GeneName

### Variant to Transcript
- **Cardinality:** 1:N
- **Description:** One variant may affect multiple transcripts
- **Key Fields:** FeatureID

### Variant to Score
- **Cardinality:** 1:1
- **Description:** One CADD score per variant
- **Key Fields:** Chrom, Pos, Ref, Alt

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CADD | Combined Annotation Dependent Depletion | Database/algorithm name |
| PHRED | Phil's Read Editor (scaled score) | -10*log10(rank/total) |
| SVM | Support Vector Machine | Classification model |
| VEP | Variant Effect Predictor | Annotation tool |
| GERP | Genomic Evolutionary Rate Profiling | Conservation metric |
| phyloP | Phylogenetic P-values | Conservation metric |
| phastCons | Phylogenetic Analysis with Space/Time | Conservation metric |
| GRCh38 | Genome Reference Consortium human 38 | Current assembly |
| GRCh37 | Genome Reference Consortium human 37 | Previous assembly |
| SNV | Single Nucleotide Variant | Variant type |
| INS | Insertion | Variant type |
| DEL | Deletion | Variant type |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| Ensembl | ENSG/ENST | Gene/transcript |
| dbSNP | rsID | Variant identifier |
| gnomAD | Variant ID | Population frequency |
| ClinVar | VCV | Clinical validation |
| dbNSFP | Variant | Prediction comparison |

---

## Score Components

| Category | Features Used |
|----------|---------------|
| Conservation | GERP, phyloP, phastCons, SiPhy |
| Regulatory | ENCODE, chromatin states |
| Protein | Grantham score, PolyPhen |
| Splice | ESE/ESS, splice site scores |
| Sequence | GC content, CpG |

---

## Data Quality Notes

1. **Cardinality:** One score per variant (not per transcript)
2. **Pre-computed:** All possible SNVs are pre-scored
3. **Indels:** Separately scored; use CADD-Indel
4. **Genome Build:** Available for both GRCh37 and GRCh38
5. **Score Meaning:** Higher PHRED = more likely deleterious
6. **Training Data:** Based on evolutionary/functional constraint
