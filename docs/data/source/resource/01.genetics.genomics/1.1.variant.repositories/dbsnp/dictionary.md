# dbSNP - Data Dictionary

## Overview

This data dictionary documents the schema for dbSNP reference SNP database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dbsnp |
| **Name** | dbSNP |
| **Parent** | 1.1.variant.repositories |
| **Total Fields** | 15 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### RefSNP

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| refsnp_id | string | 1:1 | Yes | Reference SNP identifier | rs12345, rs699 |
| create_date | date | 1:1 | Yes | First submission date | 2001-03-15 |
| last_update_date | date | 1:1 | No | Most recent update | 2024-06-20 |
| gene_symbol | string | 1:N | No | Associated gene symbol | TP53, BRCA1 |
| alleles | array | 1:N | Yes | Reference and alternate alleles | A/G, C/T/G |
| ancestral_allele | string | 1:1 | No | Inferred ancestral allele | A |

### Placement

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| seq_id | string | 1:1 | Yes | Sequence accession | NC_000001.11 |
| position | integer | 1:1 | Yes | Genomic position (0-based) | 11856378 |
| allele_string | string | 1:1 | Yes | Ref/alt alleles | G/A |
| is_aln_opposite_orientation | boolean | 1:1 | No | Alignment orientation | false |

### Population Frequency

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| study_name | string | 1:1 | Yes | Contributing study | 1000Genomes, gnomAD |
| allele_count | integer | 1:1 | Yes | Observed allele count | 1500 |
| total_count | integer | 1:1 | Yes | Total chromosomes | 5000 |
| freq | float | 1:1 | Yes | Allele frequency | 0.30 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| rsID | rs + digits | rs12345 | Reference SNP ID |
| ssID | ss + digits | ss123456789 | Submitted SNP ID |
| Seq Accession | NC_ + digits | NC_000001.11 | Chromosome sequence |
| HGVS | transcript:change | NM_000546.6:c.743G>A | Variant notation |

---

## Enumerations

### Variant Class

| Value | Meaning |
|-------|---------|
| snv | Single nucleotide variant |
| mnv | Multi-nucleotide variant |
| ins | Insertion |
| del | Deletion |
| delins | Deletion-insertion |
| identity | Reference match |

### Validation Status

| Value | Meaning |
|-------|---------|
| validated | Experimentally validated |
| by-frequency | Frequency data available |
| by-cluster | Multiple submissions |
| by-1000genomes | 1000 Genomes validated |
| by-hapmap | HapMap validated |

### Clinical Significance

| Value | Meaning |
|-------|---------|
| pathogenic | Disease-causing |
| likely-pathogenic | Probably disease-causing |
| uncertain-significance | VUS |
| likely-benign | Probably benign |
| benign | Not disease-causing |

---

## Entity Relationships

### rsID to ssID
- **Cardinality:** 1:N
- **Description:** One reference SNP aggregates multiple submitted SNPs
- **Key Fields:** refsnp_id, subsnp_id

### rsID to Gene
- **Cardinality:** N:M
- **Description:** SNPs can affect multiple genes; genes have many SNPs
- **Key Fields:** refsnp_id, gene_id

### rsID to Placement
- **Cardinality:** 1:N
- **Description:** One rsID maps to multiple genome assemblies
- **Key Fields:** refsnp_id, seq_id, assembly

### rsID to Frequency
- **Cardinality:** 1:N
- **Description:** One rsID has frequencies from multiple studies
- **Key Fields:** refsnp_id, study_name

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| SNP | Single Nucleotide Polymorphism | Common variant type |
| SNV | Single Nucleotide Variant | Broader term than SNP |
| rsID | Reference SNP Identifier | Primary dbSNP ID |
| ssID | Submitted SNP Identifier | Per-submission ID |
| MAF | Minor Allele Frequency | Population frequency |
| HGVS | Human Genome Variation Society | Nomenclature standard |
| VCF | Variant Call Format | File format |
| GRCh38 | Genome Reference Consortium human 38 | Current assembly |
| GRCh37 | Genome Reference Consortium human 37 | Previous assembly |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| ClinVar | VCV/RCV | Clinical significance |
| gnomAD | Variant ID | Population frequency |
| 1000 Genomes | Sample ID | Population data |
| Ensembl | ENST | Transcript annotation |

---

## Data Quality Notes

1. **Cardinality:** rsID→ssID forms hierarchical aggregation
2. **Assembly Mapping:** Same rsID has different coordinates per assembly
3. **Merging:** rsIDs can merge when proven identical (retired IDs redirect)
4. **Frequency Sources:** Multiple studies provide independent frequency estimates
