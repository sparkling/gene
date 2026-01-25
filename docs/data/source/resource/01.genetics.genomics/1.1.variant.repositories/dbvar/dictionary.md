# dbVar - Data Dictionary

## Overview

This data dictionary documents the schema for dbVar structural variant database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dbvar |
| **Name** | dbVar |
| **Parent** | 1.1.variant.repositories |
| **Total Fields** | 14 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| accession | string | 1:1 | Yes | dbVar variant accession | nsv1234567, nssv5678901 |
| variant_type | string | 1:1 | Yes | SV classification | copy_number_gain, deletion |
| chromosome | string | 1:1 | Yes | Chromosome | 1, X, MT |
| start | integer | 1:1 | Yes | Start position | 1000000 |
| stop | integer | 1:1 | Yes | End position | 1500000 |
| length | integer | 1:1 | No | Variant length in bp | 500000 |
| inner_start | integer | 1:1 | No | Inner boundary start | 1050000 |
| outer_stop | integer | 1:1 | No | Outer boundary end | 1550000 |

### Variant Region

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| region_accession | string | 1:1 | Yes | Region accession | nsv1234567 |
| variant_count | integer | 1:1 | No | Variants in region | 15 |
| min_start | integer | 1:1 | Yes | Region minimum start | 1000000 |
| max_stop | integer | 1:1 | Yes | Region maximum end | 1600000 |

### Sample

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| sample_id | string | 1:1 | Yes | Sample identifier | SAMN12345678 |
| study_accession | string | 1:1 | Yes | Study accession | nstd123 |
| phenotype | string | 1:N | No | Associated phenotype | intellectual disability |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| nsv | nsv + 7 digits | nsv1234567 | Variant region accession |
| nssv | nssv + 7 digits | nssv5678901 | Submitted SV accession |
| esv | esv + 7 digits | esv1234567 | DGVa variant (external) |
| essv | essv + 7 digits | essv5678901 | DGVa submitted variant |
| nstd | nstd + 3 digits | nstd123 | Study accession |

---

## Enumerations

### Variant Type

| Value | Meaning |
|-------|---------|
| copy_number_gain | Increased copy number (duplication) |
| copy_number_loss | Decreased copy number (deletion) |
| copy_number_variation | Unspecified CNV |
| deletion | Sequence deletion |
| duplication | Sequence duplication |
| insertion | Sequence insertion |
| inversion | Sequence inversion |
| translocation | Chromosomal translocation |
| complex | Complex rearrangement |
| mobile_element_insertion | Transposable element |
| tandem_duplication | Adjacent duplication |

### Clinical Significance

| Value | Meaning |
|-------|---------|
| Pathogenic | Disease-causing |
| Likely pathogenic | Probably disease-causing |
| Uncertain significance | VUS |
| Likely benign | Probably benign |
| Benign | Not disease-causing |

### Validation Status

| Value | Meaning |
|-------|---------|
| validated | Experimentally confirmed |
| not validated | Not experimentally confirmed |
| unknown | Validation status unknown |

---

## Entity Relationships

### nsv to nssv
- **Cardinality:** 1:N
- **Description:** One variant region contains multiple submitted variants
- **Key Fields:** region_accession, submitted_accession

### Variant to Study
- **Cardinality:** N:1
- **Description:** Many variants belong to one study
- **Key Fields:** accession, study_accession

### Variant to Sample
- **Cardinality:** N:1
- **Description:** Many variants from one sample
- **Key Fields:** accession, sample_id

### Variant to Gene
- **Cardinality:** N:M
- **Description:** SVs can span multiple genes
- **Key Fields:** accession, gene_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| SV | Structural Variant | Variants >50bp |
| CNV | Copy Number Variant | Gain or loss |
| nsv | NCBI Structural Variant | Region accession |
| nssv | NCBI Submitted Structural Variant | Submission accession |
| nstd | NCBI Structural Study | Study accession |
| DGVa | Database of Genomic Variants archive | Partner database |
| VCF | Variant Call Format | File format |
| BED | Browser Extensible Data | Coordinate format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| ClinVar | VCV | Clinical significance |
| DGVa | esv/essv | Partner SVs |
| dbSNP | rsID | Small variants |
| BioSample | SAMN | Sample metadata |
| BioProject | PRJNA | Study metadata |

---

## Data Quality Notes

1. **Cardinality:** nsv→nssv aggregates submitted variants into regions
2. **Coordinate Uncertainty:** Inner/outer boundaries indicate breakpoint uncertainty
3. **Size Range:** Focuses on variants >50bp (structural variants)
4. **Clinical Data:** Subset of variants have clinical assertions from ClinVar
