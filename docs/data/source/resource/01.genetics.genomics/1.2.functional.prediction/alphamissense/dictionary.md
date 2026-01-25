# AlphaMissense - Data Dictionary

## Overview

This data dictionary documents the schema for AlphaMissense pathogenicity prediction database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | alphamissense |
| **Name** | AlphaMissense |
| **Parent** | 1.2.functional.prediction |
| **Total Fields** | 10 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant Prediction

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| CHROM | string | 1:1 | Yes | Chromosome | chr1, chrX |
| POS | integer | 1:1 | Yes | Genomic position (GRCh38) | 11856378 |
| REF | string | 1:1 | Yes | Reference allele | A |
| ALT | string | 1:1 | Yes | Alternate allele | G |
| uniprot_id | string | 1:1 | Yes | UniProt accession | P04637 |
| transcript_id | string | 1:1 | Yes | Ensembl transcript | ENST00000269305 |
| protein_variant | string | 1:1 | Yes | Amino acid change | R248Q |
| am_pathogenicity | float | 1:1 | Yes | Pathogenicity score (0-1) | 0.9876 |
| am_class | string | 1:1 | Yes | Classification category | likely_pathogenic |

### Gene Summary

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_name | string | 1:1 | Yes | Gene symbol | TP53 |
| mean_am_pathogenicity | float | 1:1 | No | Average gene score | 0.45 |
| pathogenic_fraction | float | 1:1 | No | Fraction pathogenic variants | 0.32 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| UniProt | 6-10 alphanumeric | P04637 | Protein accession |
| Ensembl Transcript | ENST + 11 digits | ENST00000269305 | Transcript ID |
| Protein Variant | AA + pos + AA | R248Q | Amino acid change |
| Variant Key | chr-pos-ref-alt | chr17-7577121-G-A | Genomic coordinates |

---

## Enumerations

### am_class (Classification)

| Value | Score Range | Meaning |
|-------|-------------|---------|
| likely_benign | <0.34 | Predicted benign |
| ambiguous | 0.34-0.564 | Uncertain prediction |
| likely_pathogenic | >0.564 | Predicted pathogenic |

---

## Entity Relationships

### Variant to Transcript
- **Cardinality:** 1:1
- **Description:** Each prediction is for one transcript context
- **Key Fields:** transcript_id

### Variant to Protein
- **Cardinality:** N:1
- **Description:** Multiple variants map to one protein
- **Key Fields:** uniprot_id

### Variant to Gene
- **Cardinality:** N:1
- **Description:** Multiple variants per gene
- **Key Fields:** gene_name

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| AM | AlphaMissense | Prediction tool name |
| ESM | Evolutionary Scale Modeling | Language model base |
| AF | AlphaFold | Structural prediction |
| MSA | Multiple Sequence Alignment | Input feature |
| GRCh38 | Genome Reference Consortium human 38 | Coordinate system |
| VEP | Variant Effect Predictor | Annotation tool |
| TSV | Tab-Separated Values | File format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UniProt | Accession | Protein sequence |
| Ensembl | ENST | Transcript annotation |
| dbSNP | rsID | Variant identifier |
| ClinVar | VCV | Clinical validation |
| gnomAD | Variant ID | Population frequency |

---

## Score Interpretation

| Score Range | Classification | Action |
|-------------|----------------|--------|
| 0.0-0.34 | Likely Benign | Low concern |
| 0.34-0.564 | Ambiguous | Requires additional evidence |
| 0.564-1.0 | Likely Pathogenic | High concern, validate |

---

## Data Quality Notes

1. **Cardinality:** One prediction per missense variant per transcript
2. **Coverage:** 71M possible single amino acid substitutions
3. **Model Basis:** Combines AlphaFold structure with ESM embeddings
4. **Limitations:** Missense only; no indels, splice, or structural variants
5. **Threshold Choice:** 0.564 threshold balances sensitivity/specificity
