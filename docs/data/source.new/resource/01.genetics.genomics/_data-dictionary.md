# Category 01: Genetics & Genomics - Data Dictionary

## Overview

| Attribute | Value |
|-----------|-------|
| Category ID | 01 |
| Category Name | Genetics & Genomics |
| Total Data Sources | 25 |
| Total Subcategories | 6 |
| Total Unified Fields | 34 |
| Total Source-Specific Fields | 287 |
| Schema Version | 1.0 |

## Subcategories

| ID | Name | Data Sources | Unified Schema |
|----|------|--------------|----------------|
| 1.1 | Variant Repositories | ClinVar, dbSNP, dbVar | [1.1.variant.repositories/_data-dictionary.md](1.1.variant.repositories/_data-dictionary.md) |
| 1.2 | Functional Prediction | AlphaMissense, CADD, dbNSFP, SpliceAI | [1.2.functional.prediction/_data-dictionary.md](1.2.functional.prediction/_data-dictionary.md) |
| 1.3 | Population Genetics | 1000 Genomes, gnomAD, TOPMed, UK Biobank | [1.3.population.genetics/_data-dictionary.md](1.3.population.genetics/_data-dictionary.md) |
| 1.4 | Pharmacogenomics | CPIC, DPWG, PharmGKB, PharmVar | [1.4.pharmacogenomics/_data-dictionary.md](1.4.pharmacogenomics/_data-dictionary.md) |
| 1.5 | Expression & Regulation | ENCODE, eQTLGen, GTEx, GWAS Catalog | [1.5.expression.regulation/_data-dictionary.md](1.5.expression.regulation/_data-dictionary.md) |
| 1.6 | Cancer Genomics | BRCA Exchange, Cancer Gene Census, cBioPortal, CIViC, COSMIC, OncoKB | [1.6.cancer.genomics/_data-dictionary.md](1.6.cancer.genomics/_data-dictionary.md) |

---

## Cross-Subcategory Common Fields

### Genomic Coordinates

Fields shared across variant, functional prediction, population genetics, expression, and cancer databases.

| Field | Type | Description | Present In |
|-------|------|-------------|------------|
| chromosome | string | Chromosome identifier (1-22, X, Y, MT) | 1.1, 1.2, 1.3, 1.5, 1.6 |
| position | integer | 1-based genomic coordinate on reference assembly | 1.1, 1.2, 1.3, 1.5, 1.6 |
| reference_allele | string | Reference genome allele (IUPAC nucleotides) | 1.1, 1.2, 1.3, 1.5, 1.6 |
| alternate_allele | string | Observed variant allele (IUPAC nucleotides) | 1.1, 1.2, 1.3, 1.5, 1.6 |

### Gene Identifiers

Gene identification fields shared across all genomics databases.

| Field | Type | Description | Present In |
|-------|------|-------------|------------|
| gene_symbol | string | HGNC official gene symbol | 1.1, 1.2, 1.3, 1.4, 1.5, 1.6 |
| entrez_gene_id | integer | NCBI Entrez Gene identifier | 1.1, 1.2, 1.3, 1.4, 1.5, 1.6 |
| ensembl_gene_id | string | Ensembl gene identifier | 1.1, 1.2, 1.3, 1.4, 1.5, 1.6 |

### Clinical Significance

Variant pathogenicity classifications from clinical and cancer databases.

| Field | Type | Description | Present In |
|-------|------|-------------|------------|
| clinical_significance | string | ACMG/AMP pathogenicity classification | 1.1, 1.2, 1.6 |
| pathogenicity | string | Variant pathogenicity assessment | 1.1, 1.2, 1.6 |
| oncogenicity | string | Cancer oncogenicity classification | 1.1, 1.6 |

### Population Frequency

Population genetics metrics shared between databases.

| Field | Type | Description | Present In |
|-------|------|-------------|------------|
| allele_frequency | float | Proportion of alternate allele in population (0-1) | 1.2, 1.3 |
| allele_count | integer | Number of observed alternate alleles | 1.2, 1.3 |
| allele_number | integer | Total alleles in called genotypes | 1.2, 1.3 |

---

## Field Type Reference

| Type | Description | Examples |
|------|-------------|----------|
| string | Text value | "BRCA1", "rs334", "Pathogenic" |
| integer | Whole number | 5227002, 42, 1000 |
| float | Decimal number | 0.425, 1e-12, 0.95 |
| boolean | True/false | true, false |
| array | List of values | ["BRCA1", "BRCA2"], [123, 456] |
| object | Nested structure | {"name": "value", "id": 123} |
| date | ISO 8601 date | "2024-01-15" |
| datetime | ISO 8601 datetime | "2024-01-15T10:30:00Z" |

## Cardinality Reference

| Cardinality | Description |
|-------------|-------------|
| 1:1 | One value per record (required) |
| 1:N | One or more values per record |
| 0:1 | Zero or one value (optional) |
| 0:N | Zero or more values (optional array) |

---

## Data Source Coverage Summary

| Subcategory | Sources |
|-------------|---------|
| 1.1 Variant Repositories | 3 |
| 1.2 Functional Prediction | 4 |
| 1.3 Population Genetics | 4 |
| 1.4 Pharmacogenomics | 4 |
| 1.5 Expression & Regulation | 4 |
| 1.6 Cancer Genomics | 6 |

---

*Generated: 2026-01-24*
