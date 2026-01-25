# GWAS Catalog - Data Dictionary

## Overview

This data dictionary documents the schema for NHGRI-EBI GWAS Catalog.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gwas.catalog |
| **Name** | GWAS Catalog |
| **Parent** | 1.5.expression.regulation |
| **Total Fields** | 35+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ASSOCIATION_ID | integer | 1:1 | Yes | Unique association ID | 12345678 |
| SNPS | string | 1:N | Yes | Variant rsID(s) | rs12345, rs12345; rs67890 |
| STRONGEST_SNP_RISK_ALLELE | string | 1:1 | No | Risk allele | rs12345-A |
| P-VALUE | float | 1:1 | Yes | Association p-value | 1.0e-10 |
| P-VALUE_MLOG | float | 1:1 | Yes | -log10(p-value) | 10 |
| OR_or_BETA | float | 1:1 | No | Effect size | 1.25 |
| CI_95 | string | 1:1 | No | 95% confidence interval | [1.15-1.35] |
| RISK_ALLELE_FREQUENCY | float | 1:1 | No | Risk allele frequency | 0.35 |

### Study

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| STUDY_ACCESSION | string | 1:1 | Yes | Study identifier | GCST000001 |
| PUBMEDID | integer | 1:1 | Yes | PubMed ID | 12345678 |
| FIRST_AUTHOR | string | 1:1 | Yes | First author | Smith AB |
| DATE | date | 1:1 | Yes | Publication date | 2020-01-15 |
| JOURNAL | string | 1:1 | Yes | Journal name | Nature Genetics |
| INITIAL_SAMPLE_SIZE | string | 1:1 | Yes | Discovery sample | 50,000 European |
| REPLICATION_SAMPLE_SIZE | string | 1:1 | No | Replication sample | 25,000 European |

### Trait

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| DISEASE_TRAIT | string | 1:1 | Yes | Reported trait | Type 2 diabetes |
| MAPPED_TRAIT | string | 1:N | No | EFO mapped trait | type II diabetes mellitus |
| MAPPED_TRAIT_URI | string | 1:N | No | EFO URI | http://www.ebi.ac.uk/efo/EFO_0001360 |
| EFO_TRAIT_ID | string | 1:N | No | EFO identifier | EFO_0001360 |

### Gene Mapping

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| REPORTED_GENE | string | 1:N | No | Author-reported genes | TP53, BRCA1 |
| MAPPED_GENE | string | 1:N | No | Catalog-mapped genes | TP53 |
| UPSTREAM_GENE | string | 1:1 | No | Nearest upstream gene | GENE1 |
| DOWNSTREAM_GENE | string | 1:1 | No | Nearest downstream gene | GENE2 |
| SNP_GENE_IDS | string | 1:N | No | Ensembl gene IDs | ENSG00000141510 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Study Accession | GCST + digits | GCST000001 | Study identifier |
| Association ID | Integer | 12345678 | Association identifier |
| rsID | rs + digits | rs12345 | Variant identifier |
| EFO ID | EFO_ + digits | EFO_0001360 | Trait ontology |
| PubMed ID | Integer | 12345678 | Publication reference |

---

## Enumerations

### Study Types

| Type | Description |
|------|-------------|
| GWAS | Standard genome-wide association |
| Meta-analysis | Multiple study meta-analysis |
| Targeted array | Custom/focused genotyping |
| Sequencing | WGS/WES based |
| Imputation | Imputed variants |

### Association Context

| Context | Description |
|---------|-------------|
| Initial | Discovery sample |
| Replication | Replication sample |
| Combined | Discovery + replication |

### Risk Allele Direction

| Direction | Meaning |
|-----------|---------|
| Increase | Increases trait/risk |
| Decrease | Decreases trait/risk |
| Unknown | Direction not reported |

### Ancestry Categories

| Category | Description |
|----------|-------------|
| European | European ancestry |
| East Asian | East Asian ancestry |
| African | African ancestry |
| South Asian | South Asian ancestry |
| Hispanic/Latino | Hispanic/Latino ancestry |
| Mixed | Multiple ancestries |

---

## Entity Relationships

### Study to Association
- **Cardinality:** 1:N
- **Description:** One study has multiple associations
- **Key Fields:** STUDY_ACCESSION, ASSOCIATION_ID

### Association to SNP
- **Cardinality:** 1:N
- **Description:** Association may involve multiple SNPs (haplotype)
- **Key Fields:** ASSOCIATION_ID, SNPS

### Association to Trait
- **Cardinality:** N:M
- **Description:** Associations map to multiple EFO traits
- **Key Fields:** ASSOCIATION_ID, EFO_TRAIT_ID

### SNP to Gene
- **Cardinality:** N:M
- **Description:** SNPs mapped to multiple nearby genes
- **Key Fields:** SNPS, MAPPED_GENE

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GWAS | Genome-Wide Association Study | Study type |
| EFO | Experimental Factor Ontology | Trait ontology |
| GCST | GWAS Catalog Study | Study accession prefix |
| OR | Odds Ratio | Effect size for binary traits |
| CI | Confidence Interval | Uncertainty range |
| MAF | Minor Allele Frequency | Population metric |
| LD | Linkage Disequilibrium | Haplotype structure |
| SNP | Single Nucleotide Polymorphism | Variant type |
| NHGRI | National Human Genome Research Institute | Funder |
| EBI | European Bioinformatics Institute | Host |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubMed | PMID | Publication |
| dbSNP | rsID | Variant identifier |
| EFO | EFO ID | Trait ontology |
| Ensembl | ENSG | Gene mapping |
| gnomAD | Variant ID | Population frequency |
| Open Targets | Association | Drug targets |

---

## Trait Categories (EFO)

| Category | Example Traits |
|----------|---------------|
| Disease | Type 2 diabetes, Alzheimer's |
| Measurement | Height, BMI, blood pressure |
| Biomarker | LDL cholesterol, HbA1c |
| Response | Drug response, treatment outcome |
| Cancer | Breast cancer, prostate cancer |
| Cardiovascular | Coronary artery disease |
| Psychiatric | Schizophrenia, depression |

---

## Data Quality Notes

1. **Cardinality:** One entry per study-trait-variant association
2. **Curation:** Expert-curated from published literature
3. **P-value Threshold:** Typically P < 5e-8 for inclusion
4. **Updates:** Weekly data releases
5. **Ancestry Bias:** Historically European-centric
6. **Summary Stats:** Full summary statistics for subset of studies
