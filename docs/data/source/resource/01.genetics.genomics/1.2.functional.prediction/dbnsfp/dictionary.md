# dbNSFP - Data Dictionary

## Overview

This data dictionary documents the schema for dbNSFP functional prediction database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dbnsfp |
| **Name** | dbNSFP |
| **Parent** | 1.2.functional.prediction |
| **Total Fields** | 40+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| chr | string | 1:1 | Yes | Chromosome | 1, X |
| pos(1-based) | integer | 1:1 | Yes | Genomic position | 11856378 |
| ref | string | 1:1 | Yes | Reference allele | A |
| alt | string | 1:1 | Yes | Alternate allele | G |
| aaref | string | 1:1 | Yes | Reference amino acid | R |
| aaalt | string | 1:1 | Yes | Alternate amino acid | Q |
| aapos | integer | 1:1 | Yes | Amino acid position | 248 |

### Gene Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| genename | string | 1:1 | Yes | Gene symbol | TP53 |
| Ensembl_geneid | string | 1:1 | No | Ensembl gene ID | ENSG00000141510 |
| Ensembl_transcriptid | string | 1:1 | No | Ensembl transcript | ENST00000269305 |
| Uniprot_acc | string | 1:1 | No | UniProt accession | P04637 |

### Prediction Scores

| Field Name | Data Type | Cardinality | Required | Description | Range |
|------------|-----------|-------------|----------|-------------|-------|
| SIFT_score | float | 1:1 | No | SIFT prediction score | 0-1 (lower=damaging) |
| SIFT_pred | string | 1:1 | No | SIFT prediction | D, T |
| Polyphen2_HDIV_score | float | 1:1 | No | PolyPhen-2 HDIV score | 0-1 (higher=damaging) |
| Polyphen2_HDIV_pred | string | 1:1 | No | PolyPhen-2 prediction | D, P, B |
| CADD_phred | float | 1:1 | No | CADD phred-scaled score | 0-99 |
| REVEL_score | float | 1:1 | No | REVEL ensemble score | 0-1 |
| MetaSVM_score | float | 1:1 | No | MetaSVM score | -2 to 3 |
| MetaLR_score | float | 1:1 | No | MetaLR score | 0-1 |
| MutationTaster_pred | string | 1:1 | No | MutationTaster prediction | A, D, N, P |

### Conservation Scores

| Field Name | Data Type | Cardinality | Required | Description | Range |
|------------|-----------|-------------|----------|-------------|-------|
| phyloP100way_vertebrate | float | 1:1 | No | phyloP conservation | -20 to 10 |
| phastCons100way_vertebrate | float | 1:1 | No | phastCons score | 0-1 |
| GERP++_RS | float | 1:1 | No | GERP rejected substitutions | -12 to 6 |

### Population Frequencies

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gnomAD_exomes_AF | float | 1:1 | No | gnomAD exome frequency | 0.0001 |
| gnomAD_genomes_AF | float | 1:1 | No | gnomAD genome frequency | 0.00005 |
| 1000Gp3_AF | float | 1:1 | No | 1000 Genomes frequency | 0.001 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Ensembl Gene | ENSG + 11 digits | ENSG00000141510 | Gene ID |
| Ensembl Transcript | ENST + 11 digits | ENST00000269305 | Transcript ID |
| UniProt | 6-10 alphanumeric | P04637 | Protein accession |
| rsID | rs + digits | rs28934576 | dbSNP identifier |

---

## Enumerations

### SIFT Prediction

| Value | Meaning | Score Range |
|-------|---------|-------------|
| D | Damaging | <0.05 |
| T | Tolerated | >=0.05 |

### PolyPhen-2 Prediction

| Value | Meaning | Score Range |
|-------|---------|-------------|
| D | Probably damaging | >0.957 |
| P | Possibly damaging | 0.453-0.956 |
| B | Benign | <0.453 |

### MutationTaster Prediction

| Value | Meaning |
|-------|---------|
| A | Disease-causing automatic |
| D | Disease-causing |
| N | Polymorphism |
| P | Polymorphism automatic |

### Clinical Significance (ClinVar)

| Value | Meaning |
|-------|---------|
| Pathogenic | Disease-causing |
| Likely_pathogenic | Probably pathogenic |
| Uncertain_significance | VUS |
| Likely_benign | Probably benign |
| Benign | Not pathogenic |

---

## Entity Relationships

### Variant to Transcript
- **Cardinality:** 1:N
- **Description:** One variant may affect multiple transcripts
- **Key Fields:** chr, pos, Ensembl_transcriptid

### Variant to Gene
- **Cardinality:** N:1
- **Description:** Multiple variants per gene
- **Key Fields:** genename

### Variant to Prediction Scores
- **Cardinality:** 1:1
- **Description:** One set of predictions per variant
- **Key Fields:** chr, pos, ref, alt

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| dbNSFP | Database of Non-Synonymous Functional Predictions | Database name |
| SIFT | Sorting Intolerant From Tolerant | Prediction algorithm |
| PolyPhen | Polymorphism Phenotyping | Prediction algorithm |
| CADD | Combined Annotation Dependent Depletion | Meta-predictor |
| REVEL | Rare Exome Variant Ensemble Learner | Meta-predictor |
| GERP | Genomic Evolutionary Rate Profiling | Conservation score |
| phyloP | Phylogenetic P-values | Conservation score |
| phastCons | Phylogenetic Analysis with Space/Time | Conservation score |
| gnomAD | Genome Aggregation Database | Population frequency |
| MAF | Minor Allele Frequency | Population metric |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| Ensembl | ENSG/ENST | Gene/transcript |
| UniProt | Accession | Protein |
| dbSNP | rsID | Variant identifier |
| ClinVar | VCV | Clinical data |
| gnomAD | Variant ID | Population frequency |
| OMIM | MIM | Disease associations |

---

## Meta-Predictor Guidance

| Predictor | Threshold | Interpretation |
|-----------|-----------|----------------|
| REVEL | >= 0.5 | Likely pathogenic |
| CADD phred | >= 20 | Top 1% deleterious |
| MetaSVM | > 0 | Damaging |
| MetaLR | > 0.5 | Damaging |

---

## Data Quality Notes

1. **Cardinality:** One row per variant-transcript combination
2. **Score Coverage:** Not all scores available for all variants
3. **Version Dependency:** Prediction scores may differ across database versions
4. **Missing Values:** Represented as "." in TSV format
5. **Ensemble Approach:** Combine multiple predictors for better accuracy
