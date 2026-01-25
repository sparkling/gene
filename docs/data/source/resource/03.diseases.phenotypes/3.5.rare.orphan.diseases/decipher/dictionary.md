# DECIPHER - Data Dictionary

## Overview

This data dictionary documents the schema for DECIPHER (Database of Chromosomal Imbalance and Phenotype in Humans using Ensembl Resources).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | decipher |
| **Name** | DECIPHER |
| **Parent** | 3.5.rare.orphan.diseases |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Patient Record

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| patient_id | integer | 1:1 | Yes | DECIPHER patient ID | 250000 |
| sex | enum | 1:1 | No | Patient sex | Male, Female, Unknown |
| phenotypes | array | 1:N | No | HPO phenotype terms | [HP:0001249, HP:0001250] |
| inheritance | string | 1:1 | No | Inheritance pattern | de novo, inherited |
| age_of_onset | string | 1:1 | No | Onset age category | Childhood |
| consent_level | enum | 1:1 | Yes | Data sharing consent | open, registered |

### CNV (Copy Number Variant)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| variant_id | integer | 1:1 | Yes | Variant identifier | 350000 |
| patient_id | integer | 1:1 | Yes | Patient reference | 250000 |
| chromosome | string | 1:1 | Yes | Chromosome | 22 |
| start | integer | 1:1 | Yes | Start position (GRCh38) | 19000000 |
| end | integer | 1:1 | Yes | End position | 21000000 |
| variant_type | enum | 1:1 | Yes | CNV type | deletion, duplication |
| size | integer | 1:1 | Yes | Size in base pairs | 2000000 |
| pathogenicity | enum | 1:1 | No | Classification | pathogenic, likely_pathogenic |
| genes | array | 1:N | No | Affected genes | [SHANK3, ACR] |
| inheritance | string | 1:1 | No | Inheritance | de novo |

### SNV (Single Nucleotide Variant)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| variant_id | integer | 1:1 | Yes | Variant identifier | 450000 |
| patient_id | integer | 1:1 | Yes | Patient reference | 250000 |
| chromosome | string | 1:1 | Yes | Chromosome | 7 |
| position | integer | 1:1 | Yes | Position (GRCh38) | 117559593 |
| reference | string | 1:1 | Yes | Reference allele | C |
| alternate | string | 1:1 | Yes | Alternate allele | T |
| gene | string | 1:1 | No | Affected gene | CFTR |
| consequence | string | 1:1 | No | Variant consequence | missense_variant |
| pathogenicity | enum | 1:1 | No | Classification | pathogenic |
| zygosity | enum | 1:1 | No | Zygosity | heterozygous |

### Gene Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_symbol | string | 1:1 | Yes | HGNC symbol | SHANK3 |
| ensembl_id | string | 1:1 | Yes | Ensembl gene ID | ENSG00000251322 |
| hi_score | float | 1:1 | No | Haploinsufficiency score | 0.95 |
| pli | float | 1:1 | No | pLI score | 1.0 |
| loeuf | float | 1:1 | No | LOEUF score | 0.15 |
| disease_genes | array | 1:N | No | Associated diseases | [OMIM:606232] |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| DECIPHER ID | [0-9]+ | 250000 | Patient/variant ID |
| Patient ID | DECIPHER:[0-9]+ | DECIPHER:250001 | CURIE format |
| HPO Term | HP:[0-9]{7} | HP:0001249 | Phenotype term |
| Ensembl Gene | ENSG[0-9]{11} | ENSG00000251322 | Gene identifier |
| Ensembl Transcript | ENST[0-9]{11} | ENST00000262795 | Transcript ID |
| OMIM | OMIM:[0-9]{6} | OMIM:606232 | Disease reference |

---

## Enumerations

### Variant Types (CNV)

| Type | Description |
|------|-------------|
| deletion | Copy number loss |
| duplication | Copy number gain |
| insertion | Insertion |
| inversion | Inversion |
| translocation | Translocation |

### Pathogenicity Classification

| Class | Description |
|-------|-------------|
| pathogenic | Disease-causing |
| likely_pathogenic | Probably disease-causing |
| uncertain_significance | VUS |
| likely_benign | Probably benign |
| benign | Not disease-causing |

### Inheritance Patterns

| Pattern | Description |
|---------|-------------|
| de novo | New mutation in patient |
| inherited_from_father | Paternal inheritance |
| inherited_from_mother | Maternal inheritance |
| inherited | Unknown parent |
| unknown | Not determined |

### Consent Levels

| Level | Description |
|-------|-------------|
| open | Publicly visible |
| registered | Requires DECIPHER login |
| private | Not shared |

### Sex Categories

| Value | Description |
|-------|-------------|
| Male | Male patient |
| Female | Female patient |
| Unknown | Not specified |

### Variant Consequences

| Consequence | Description |
|-------------|-------------|
| frameshift_variant | Frameshift mutation |
| stop_gained | Nonsense mutation |
| missense_variant | Amino acid change |
| splice_donor_variant | Splice site (5') |
| splice_acceptor_variant | Splice site (3') |
| start_lost | Start codon lost |
| stop_lost | Stop codon lost |

### Constraint Scores

| Score | Range | Description |
|-------|-------|-------------|
| HI Score | 0-1 | Haploinsufficiency prediction |
| pLI | 0-1 | Probability of LoF intolerance |
| LOEUF | 0-2 | LoF observed/expected upper bound |
| Z-score | -5 to +5 | Missense constraint |

---

## Entity Relationships

### Patient to Variants
- **Cardinality:** 1:N
- **Description:** Patients have multiple variants
- **Key Fields:** patient_id, variant_id

### Patient to Phenotypes
- **Cardinality:** 1:N
- **Description:** Patients have multiple phenotypes
- **Key Fields:** patient_id, phenotypes

### Variant to Genes
- **Cardinality:** N:M
- **Description:** CNVs affect multiple genes
- **Key Fields:** variant_id, gene_symbol

### Gene to Constraint Scores
- **Cardinality:** 1:1
- **Description:** Each gene has constraint metrics
- **Key Fields:** gene_symbol, hi_score

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DECIPHER | Database of Chromosomal Imbalance and Phenotype in Humans using Ensembl Resources | Full name |
| CNV | Copy Number Variant | Structural variant |
| SNV | Single Nucleotide Variant | Point mutation |
| HPO | Human Phenotype Ontology | Phenotype terms |
| VUS | Variant of Uncertain Significance | Classification |
| pLI | Probability of LoF Intolerance | Constraint metric |
| LOEUF | LoF Observed/Expected Upper Bound Fraction | Constraint metric |
| HI | Haploinsufficiency | Gene dosage |
| LoF | Loss of Function | Mutation type |
| GRCh38 | Genome Reference Consortium Human Build 38 | Reference genome |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| HPO | HP ID | Phenotype terms |
| OMIM | MIM number | Disease associations |
| Ensembl | Gene/Transcript ID | Gene annotation |
| HGNC | Gene symbol | Gene nomenclature |
| ClinVar | Variation ID | Variant overlap |
| gnomAD | Frequency data | Population frequency |
| DGV | CNV ID | Population CNVs |

---

## Data Quality Notes

1. **Patient Records:** 50,000+ anonymized patients
2. **CNV Records:** 40,000+ chromosomal variants
3. **SNV Records:** 60,000+ sequence variants
4. **Contributing Centers:** 300+ clinical genetics centers
5. **Countries:** 35+ contributing countries
6. **HPO Annotation:** Standardized phenotyping
7. **Expert Curation:** Clinical geneticist review
8. **Open Access:** Registered access required for detailed data

