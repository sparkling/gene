# SpliceAI - Data Dictionary

## Overview

This data dictionary documents the schema for SpliceAI splice site prediction database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | spliceai |
| **Name** | SpliceAI |
| **Parent** | 1.2.functional.prediction |
| **Total Fields** | 12 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| CHROM | string | 1:1 | Yes | Chromosome | 1, X |
| POS | integer | 1:1 | Yes | Genomic position | 11856378 |
| REF | string | 1:1 | Yes | Reference allele | A |
| ALT | string | 1:1 | Yes | Alternate allele | G |
| SYMBOL | string | 1:1 | Yes | Gene symbol | BRCA1 |
| STRAND | string | 1:1 | Yes | Gene strand | +, - |

### Delta Scores

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| DS_AG | float | 1:1 | Yes | Delta score acceptor gain | 0.92 |
| DS_AL | float | 1:1 | Yes | Delta score acceptor loss | 0.05 |
| DS_DG | float | 1:1 | Yes | Delta score donor gain | 0.01 |
| DS_DL | float | 1:1 | Yes | Delta score donor loss | 0.87 |
| DS_MAX | float | 1:1 | Yes | Maximum delta score | 0.92 |

### Position Offsets

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| DP_AG | integer | 1:1 | Yes | Position offset acceptor gain | -15 |
| DP_AL | integer | 1:1 | Yes | Position offset acceptor loss | 0 |
| DP_DG | integer | 1:1 | Yes | Position offset donor gain | 42 |
| DP_DL | integer | 1:1 | Yes | Position offset donor loss | -3 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Variant Key | chr-pos-ref-alt | 17-41197732-G-A | Genomic coordinates |
| Gene Symbol | HGNC symbol | BRCA1 | Gene identifier |
| Transcript | ENST + digits | ENST00000357654 | Ensembl transcript |

---

## Enumerations

### Score Interpretation

| Score Range | Impact | Clinical Relevance |
|-------------|--------|-------------------|
| 0.0-0.2 | Low | Unlikely to affect splicing |
| 0.2-0.5 | Moderate | May affect splicing |
| 0.5-0.8 | High | Likely affects splicing |
| 0.8-1.0 | Very High | Strong splicing effect |

### Splice Effect Type

| Value | Meaning |
|-------|---------|
| AG | Acceptor Gain (cryptic 3' splice site) |
| AL | Acceptor Loss (loss of 3' splice site) |
| DG | Donor Gain (cryptic 5' splice site) |
| DL | Donor Loss (loss of 5' splice site) |

---

## Entity Relationships

### Variant to Gene
- **Cardinality:** N:1
- **Description:** Multiple variants affect one gene
- **Key Fields:** SYMBOL

### Variant to Transcript
- **Cardinality:** 1:N
- **Description:** One variant may affect multiple transcripts
- **Key Fields:** transcript_id

### Score to Position
- **Cardinality:** 1:1
- **Description:** Each delta score has a corresponding position offset
- **Key Fields:** DS_*, DP_*

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DS | Delta Score | Splice probability change |
| DP | Delta Position | Distance to splice site |
| AG | Acceptor Gain | New 3' splice site |
| AL | Acceptor Loss | Lost 3' splice site |
| DG | Donor Gain | New 5' splice site |
| DL | Donor Loss | Lost 5' splice site |
| CNN | Convolutional Neural Network | Model architecture |
| VCF | Variant Call Format | File format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| Ensembl | ENST | Transcript annotation |
| HGNC | Symbol | Gene naming |
| dbSNP | rsID | Variant identifier |
| ClinVar | VCV | Clinical validation |
| gnomAD | Variant ID | Population frequency |

---

## Score Calculation

| Component | Description |
|-----------|-------------|
| Delta Score | Difference in splice probability (variant vs reference) |
| Position Offset | Distance from variant to affected splice junction |
| Max Score | Maximum of all four delta scores (DS_AG, DS_AL, DS_DG, DS_DL) |

---

## Data Quality Notes

1. **Cardinality:** One prediction per variant per gene (multi-gene variants have multiple entries)
2. **Context Window:** Considers 10kb flanking sequence
3. **Threshold Recommendation:** DS >= 0.2 for potential splice effects
4. **Position Range:** DP values range from -50 to +50 nucleotides
5. **Limitations:** Pre-mRNA based; does not model tissue-specific splicing
