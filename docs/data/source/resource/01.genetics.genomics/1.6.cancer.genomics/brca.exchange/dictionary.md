# BRCA Exchange - Data Dictionary

## Overview

This data dictionary documents the schema for BRCA Exchange variant database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | brca.exchange |
| **Name** | BRCA Exchange |
| **Parent** | 1.6.cancer.genomics |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | Internal variant ID | 12345 |
| Gene_Symbol | string | 1:1 | Yes | Gene (BRCA1/BRCA2) | BRCA1 |
| HGVS_cDNA | string | 1:1 | Yes | cDNA notation | NM_007294.4:c.5266dupC |
| HGVS_Protein | string | 1:1 | No | Protein notation | NP_009225.1:p.Gln1756ProfsTer74 |
| Chr | string | 1:1 | Yes | Chromosome | 17 |
| Pos | integer | 1:1 | Yes | GRCh38 position | 43057051 |
| Ref | string | 1:1 | Yes | Reference allele | C |
| Alt | string | 1:1 | Yes | Alternate allele | CC |

### Clinical Classification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Pathogenicity_expert | string | 1:1 | Yes | Expert classification | Pathogenic |
| Clinical_significance_ENIGMA | string | 1:1 | No | ENIGMA classification | Pathogenic |
| Clinical_Significance_ClinVar | string | 1:1 | No | ClinVar classification | Pathogenic |
| Date_last_evaluated_ENIGMA | date | 1:1 | No | ENIGMA evaluation date | 2023-06-15 |
| Assertion_method_ENIGMA | string | 1:1 | No | Classification method | ENIGMA BRCA1/2 Classification Criteria |
| IARC_class_ENIGMA | string | 1:1 | No | IARC class | 5 |

### Population Frequencies

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Allele_frequency_gnomAD | float | 1:1 | No | gnomAD frequency | 0.00001 |
| Allele_frequency_ExAC | float | 1:1 | No | ExAC frequency | 0.00002 |
| Allele_count_gnomAD | integer | 1:1 | No | gnomAD count | 2 |

### Functional Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Functional_analysis_result_LOVD | string | 1:1 | No | LOVD functional result | affects function |
| Functional_analysis_technique_LOVD | string | 1:1 | No | Assay method | HDR assay |
| BIC_Designation | string | 1:1 | No | BIC nomenclature | 5382insC |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| BRCA Exchange ID | Integer | 12345 | Internal variant ID |
| HGVS cDNA | NM_:c. notation | NM_007294.4:c.5266dupC | Coding change |
| HGVS Protein | NP_:p. notation | NP_009225.1:p.Gln1756fs | Protein change |
| ClinVar ID | VCV + digits | VCV000017661 | ClinVar accession |
| BIC | Legacy notation | 5382insC | BIC designation |
| rsID | rs + digits | rs80357906 | dbSNP identifier |

---

## Enumerations

### Pathogenicity (Expert)

| Classification | Description |
|----------------|-------------|
| Pathogenic | Disease-causing |
| Likely Pathogenic | >90% certainty pathogenic |
| Uncertain Significance | VUS |
| Likely Benign | >90% certainty benign |
| Benign | Not disease-causing |
| Not Yet Classified | Pending review |

### IARC Classification

| Class | Probability | Meaning |
|-------|-------------|---------|
| 5 | >0.99 | Definitely pathogenic |
| 4 | 0.95-0.99 | Likely pathogenic |
| 3 | 0.05-0.949 | Uncertain |
| 2 | 0.001-0.049 | Likely not pathogenic |
| 1 | <0.001 | Not pathogenic |

### Data Sources

| Source | Description |
|--------|-------------|
| ENIGMA | Expert panel consortium |
| ClinVar | NCBI clinical database |
| LOVD | Leiden Open Variation Database |
| BIC | Breast Cancer Information Core |
| gnomAD | Population frequencies |
| ExAC | Legacy population data |
| 1000 Genomes | Population reference |

### Variant Type

| Type | Description |
|------|-------------|
| substitution | Single nucleotide change |
| deletion | Nucleotide deletion |
| insertion | Nucleotide insertion |
| duplication | Nucleotide duplication |
| delins | Deletion-insertion |

---

## Entity Relationships

### Variant to Gene
- **Cardinality:** N:1
- **Description:** Variants belong to BRCA1 or BRCA2
- **Key Fields:** Gene_Symbol

### Variant to Classification
- **Cardinality:** 1:N
- **Description:** Multiple source classifications per variant
- **Key Fields:** id, source

### Variant to Population
- **Cardinality:** 1:N
- **Description:** Frequencies from multiple populations
- **Key Fields:** id, population

### Variant to Evidence
- **Cardinality:** 1:N
- **Description:** Multiple evidence sources per variant
- **Key Fields:** id, source

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| BRCA | Breast Cancer gene | Gene family |
| BRCA1 | Breast Cancer 1 | Gene on chr17 |
| BRCA2 | Breast Cancer 2 | Gene on chr13 |
| ENIGMA | Evidence-based Network for Interpretation of Germline Mutant Alleles | Expert consortium |
| IARC | International Agency for Research on Cancer | Classification system |
| VUS | Variant of Uncertain Significance | Classification category |
| HGVS | Human Genome Variation Society | Nomenclature |
| BIC | Breast Cancer Information Core | Legacy database |
| LOVD | Leiden Open Variation Database | Variant database |
| GA4GH | Global Alliance for Genomics and Health | Standards body |
| HBOC | Hereditary Breast and Ovarian Cancer | Syndrome |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| ClinVar | VCV | Clinical classification |
| dbSNP | rsID | Variant identifier |
| gnomAD | Variant ID | Population frequency |
| LOVD | Variant ID | Functional data |
| BIC | BIC Designation | Legacy notation |
| Ensembl | ENST | Transcript reference |
| OMIM | MIM | Disease reference |

---

## Gene Coverage

| Gene | Chromosome | RefSeq | Variants |
|------|------------|--------|----------|
| BRCA1 | 17q21.31 | NM_007294.4 | 8,000+ |
| BRCA2 | 13q13.1 | NM_000059.4 | 10,000+ |

---

## Data Quality Notes

1. **Cardinality:** One expert classification per variant
2. **Scope:** BRCA1/BRCA2 variants only
3. **Expert Curation:** ENIGMA consortium review
4. **Classification Criteria:** ACMG/AMP adapted for BRCA
5. **Data Sharing:** GA4GH Beacon-compatible
6. **Updates:** Regular sync with ClinVar/ENIGMA
7. **Historical:** Integrates legacy BIC designations
