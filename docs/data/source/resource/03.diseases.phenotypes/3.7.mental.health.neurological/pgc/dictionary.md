# PGC - Data Dictionary

## Overview

This data dictionary documents the schema for PGC (Psychiatric Genomics Consortium).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pgc |
| **Name** | PGC |
| **Parent** | 3.7.mental.health.neurological |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### GWAS Summary Statistics

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| snp_id | string | 1:1 | Yes | dbSNP identifier | rs1006737 |
| chromosome | string | 1:1 | Yes | Chromosome | 12 |
| position | integer | 1:1 | Yes | Base position | 2345678 |
| effect_allele | string | 1:1 | Yes | Effect allele | A |
| other_allele | string | 1:1 | Yes | Reference allele | G |
| effect_allele_freq | float | 1:1 | No | Effect allele frequency | 0.35 |
| beta | float | 1:1 | Yes | Effect size (log OR) | 0.045 |
| se | float | 1:1 | Yes | Standard error | 0.008 |
| p_value | float | 1:1 | Yes | Association p-value | 2.5e-12 |
| n_cases | integer | 1:1 | No | Number of cases | 50000 |
| n_controls | integer | 1:1 | No | Number of controls | 150000 |

### Study Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| study_id | string | 1:1 | Yes | PGC study ID | PGC_SCZ3 |
| disorder | string | 1:1 | Yes | Psychiatric disorder | Schizophrenia |
| total_samples | integer | 1:1 | Yes | Total sample size | 320000 |
| n_cases | integer | 1:1 | Yes | Case count | 76755 |
| n_controls | integer | 1:1 | Yes | Control count | 243649 |
| ancestry | string | 1:1 | No | Primary ancestry | European |
| publication_pmid | integer | 1:1 | No | PubMed ID | 35396580 |
| publication_year | integer | 1:1 | No | Publication year | 2022 |
| genome_build | string | 1:1 | Yes | Reference build | GRCh37 |

### Locus Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| locus_id | integer | 1:1 | Yes | Locus number | 45 |
| study_id | string | 1:1 | Yes | Study reference | PGC_SCZ3 |
| lead_snp | string | 1:1 | Yes | Lead variant | rs1006737 |
| chromosome | string | 1:1 | Yes | Chromosome | 12 |
| start | integer | 1:1 | Yes | Locus start | 2250000 |
| end | integer | 1:1 | Yes | Locus end | 2450000 |
| genes | array | 1:N | No | Genes in locus | [CACNA1C] |
| p_value | float | 1:1 | Yes | Lead SNP p-value | 4.3e-18 |
| novel | boolean | 1:1 | No | Novel discovery | true |

### Genetic Correlation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| trait1 | string | 1:1 | Yes | First trait | Schizophrenia |
| trait2 | string | 1:1 | Yes | Second trait | Bipolar Disorder |
| rg | float | 1:1 | Yes | Genetic correlation | 0.68 |
| se | float | 1:1 | Yes | Standard error | 0.04 |
| p_value | float | 1:1 | Yes | P-value | 1.2e-45 |
| method | string | 1:1 | Yes | Method used | LDSC |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| dbSNP ID | rs[0-9]+ | rs1006737 | Variant identifier |
| Study ID | PGC_[A-Z]+[0-9]* | PGC_SCZ3 | PGC study code |
| PMID | Integer | 35396580 | PubMed reference |
| Gene Symbol | Text | CACNA1C | HGNC symbol |
| Disorder Code | Text | SCZ | Disorder abbreviation |

---

## Enumerations

### Psychiatric Disorders

| Code | Name | Latest GWAS Loci |
|------|------|------------------|
| SCZ | Schizophrenia | 287 loci |
| BIP | Bipolar Disorder | 64 loci |
| MDD | Major Depressive Disorder | 243 loci |
| ADHD | Attention Deficit Hyperactivity Disorder | 27 loci |
| ASD | Autism Spectrum Disorder | 5 loci |
| AN | Anorexia Nervosa | 8 loci |
| OCD | Obsessive Compulsive Disorder | 1 locus |
| TS | Tourette Syndrome | 1 locus |
| PTSD | Post-Traumatic Stress Disorder | 3 loci |
| ANX | Anxiety Disorders | Multiple |
| SUD | Substance Use Disorders | Multiple |

### Study Versions

| Study | Version | Sample Size | Reference |
|-------|---------|-------------|-----------|
| SCZ3 | Wave 3 | 320,000 | Trubetskoy 2022 |
| BIP3 | Wave 3 | 413,000 | Mullins 2021 |
| MDD3 | Wave 3 | 1,200,000 | Howard 2019 |
| ADHD2 | Wave 2 | 225,000 | Demontis 2023 |
| ASD | Wave 1 | 46,000 | Grove 2019 |

### Ancestry Categories

| Ancestry | Description |
|----------|-------------|
| EUR | European |
| EAS | East Asian |
| AFR | African |
| AMR | Admixed American |
| SAS | South Asian |
| MULTI | Multi-ancestry |

### Genome Builds

| Build | Description |
|-------|-------------|
| GRCh37 | hg19 |
| GRCh38 | hg38 |

### Effect Measures

| Measure | Description |
|---------|-------------|
| OR | Odds Ratio |
| Beta | Log odds ratio |
| Z | Z-score |

### Significance Thresholds

| Threshold | Value | Description |
|-----------|-------|-------------|
| Genome-wide | 5e-8 | Standard GWAS threshold |
| Suggestive | 1e-5 | Suggestive significance |
| Bonferroni | Varies | Multiple testing correction |

---

## Entity Relationships

### Study to Variants
- **Cardinality:** 1:N
- **Description:** Studies report many variant associations
- **Key Fields:** study_id, snp_id

### Study to Loci
- **Cardinality:** 1:N
- **Description:** Studies identify multiple loci
- **Key Fields:** study_id, locus_id

### Locus to Genes
- **Cardinality:** 1:N
- **Description:** Loci contain multiple genes
- **Key Fields:** locus_id, genes

### Disorder Pairs to Correlations
- **Cardinality:** 1:1
- **Description:** Genetic correlations between disorders
- **Key Fields:** trait1, trait2

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PGC | Psychiatric Genomics Consortium | Organization |
| GWAS | Genome-Wide Association Study | Study type |
| SCZ | Schizophrenia | Disorder |
| BIP | Bipolar Disorder | Disorder |
| MDD | Major Depressive Disorder | Disorder |
| ADHD | Attention Deficit Hyperactivity Disorder | Disorder |
| ASD | Autism Spectrum Disorder | Disorder |
| OR | Odds Ratio | Effect size |
| SE | Standard Error | Uncertainty |
| MAF | Minor Allele Frequency | Population frequency |
| LDSC | LD Score Regression | Analysis method |
| SNP | Single Nucleotide Polymorphism | Variant type |
| LD | Linkage Disequilibrium | Genetic correlation |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant data |
| GWAS Catalog | Study ID | Study metadata |
| PubMed | PMID | Publications |
| Ensembl | Gene ID | Gene annotation |
| HGNC | Gene symbol | Gene nomenclature |
| LD Hub | Trait ID | Genetic correlations |
| Open Targets | Disease ID | Target-disease |

---

## Data Quality Notes

1. **Participating Groups:** 800+ research groups
2. **Total Samples:** 4M+ across all studies
3. **Published GWAS:** 50+ major publications
4. **Genetic Loci:** 1,000+ identified loci
5. **Cross-Disorder:** Extensive genetic overlap analysis
6. **Multi-Ancestry:** Growing non-European samples
7. **Summary Stats:** Available for download
8. **Open Science:** Data sharing policies

