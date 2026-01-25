# COSMIC - Data Dictionary

## Overview

This data dictionary documents the schema for COSMIC (Catalogue of Somatic Mutations in Cancer) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | cosmic |
| **Name** | COSMIC |
| **Parent** | 1.6.cancer.genomics |
| **Total Fields** | 35+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Mutation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| COSMIC_MUTATION_ID | string | 1:1 | Yes | Legacy mutation ID | COSM476 |
| GENOMIC_MUTATION_ID | string | 1:1 | Yes | Genomic mutation ID | COSV52071926 |
| LEGACY_MUTATION_ID | string | 1:1 | No | Legacy COSM ID | COSM476 |
| GENE_SYMBOL | string | 1:1 | Yes | Gene symbol | BRAF |
| ACCESSION_NUMBER | string | 1:1 | No | Transcript accession | ENST00000288602 |
| MUTATION_CDS | string | 1:1 | Yes | CDS-level change | c.1799T>A |
| MUTATION_AA | string | 1:1 | No | Protein change | p.V600E |
| MUTATION_DESCRIPTION | string | 1:1 | Yes | Consequence type | Substitution - Missense |
| FATHMM_PREDICTION | string | 1:1 | No | Pathogenicity prediction | PATHOGENIC |
| FATHMM_SCORE | float | 1:1 | No | FATHMM score | 0.97 |
| MUTATION_SOMATIC_STATUS | string | 1:1 | No | Somatic confirmation | Confirmed somatic |

### Sample

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| COSMIC_SAMPLE_ID | integer | 1:1 | Yes | Sample identifier | 1234567 |
| SAMPLE_NAME | string | 1:1 | Yes | Sample name | TCGA-AA-3518-01 |
| PRIMARY_SITE | string | 1:1 | Yes | Primary tumor site | large_intestine |
| SITE_SUBTYPE_1 | string | 1:1 | No | Anatomical subtype | colon |
| PRIMARY_HISTOLOGY | string | 1:1 | Yes | Histology | carcinoma |
| HISTOLOGY_SUBTYPE_1 | string | 1:1 | No | Histology subtype | adenocarcinoma |
| TUMOUR_ORIGIN | string | 1:1 | No | Tumor origin | primary |
| AGE | integer | 1:1 | No | Patient age | 65 |

### Gene

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| COSMIC_GENE_ID | integer | 1:1 | Yes | COSMIC gene ID | 376 |
| GENE_SYMBOL | string | 1:1 | Yes | Gene symbol | TP53 |
| GENE_NAME | string | 1:1 | Yes | Full gene name | tumor protein p53 |
| ENTREZ_GENE_ID | integer | 1:1 | No | NCBI Gene ID | 7157 |
| GENOME_LOCATION | string | 1:1 | No | Cytogenetic band | 17p13.1 |

### CNV

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| CNV_ID | integer | 1:1 | Yes | CNV identifier | 12345 |
| GENE_SYMBOL | string | 1:1 | Yes | Gene symbol | ERBB2 |
| CNV_TYPE | string | 1:1 | Yes | Gain or loss | gain |
| TOTAL_CN | integer | 1:1 | No | Total copy number | 8 |
| MINOR_CN | integer | 1:1 | No | Minor allele CN | 0 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| COSM | COSM + digits | COSM476 | Legacy mutation ID |
| COSV | COSV + digits | COSV52071926 | Genomic mutation ID |
| COSMIC Sample | Integer | 1234567 | Sample identifier |
| COSMIC Gene | Integer | 376 | Gene identifier |
| Study ID | Integer | 500 | Study identifier |

---

## Enumerations

### Mutation Description

| Type | Description |
|------|-------------|
| Substitution - Missense | Amino acid change |
| Substitution - Nonsense | Stop codon |
| Insertion - Frameshift | Frameshift insertion |
| Deletion - Frameshift | Frameshift deletion |
| Insertion - In frame | In-frame insertion |
| Deletion - In frame | In-frame deletion |
| Complex | Complex change |
| Unknown | Unknown effect |

### Somatic Status

| Status | Description |
|--------|-------------|
| Confirmed somatic variant | Validated somatic |
| Reported in another sample | Seen elsewhere |
| Variant of unknown origin | Not confirmed |

### FATHMM Prediction

| Prediction | Description |
|------------|-------------|
| PATHOGENIC | Likely pathogenic |
| NEUTRAL | Likely neutral |
| - | Not scored |

### Primary Site Categories

| Site | Examples |
|------|----------|
| breast | Breast cancer |
| lung | Lung cancer |
| large_intestine | Colorectal cancer |
| haematopoietic_and_lymphoid | Blood cancers |
| skin | Melanoma, SCC |
| CNS | Brain tumors |
| ovary | Ovarian cancer |
| prostate | Prostate cancer |

### Tumor Origin

| Origin | Description |
|--------|-------------|
| primary | Original tumor |
| metastasis | Metastatic site |
| recurrence | Recurrent tumor |
| secondary | Secondary primary |

---

## Entity Relationships

### Mutation to Sample
- **Cardinality:** N:M
- **Description:** Mutations observed in multiple samples
- **Key Fields:** COSMIC_MUTATION_ID, COSMIC_SAMPLE_ID

### Mutation to Gene
- **Cardinality:** N:1
- **Description:** Multiple mutations per gene
- **Key Fields:** COSMIC_MUTATION_ID, COSMIC_GENE_ID

### Sample to Study
- **Cardinality:** N:1
- **Description:** Multiple samples per study
- **Key Fields:** COSMIC_SAMPLE_ID, STUDY_ID

### Gene to CNV
- **Cardinality:** 1:N
- **Description:** One gene has multiple CNV entries
- **Key Fields:** COSMIC_GENE_ID, CNV_ID

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| COSMIC | Catalogue of Somatic Mutations in Cancer | Database name |
| COSM | COSMIC Mutation | Legacy mutation prefix |
| COSV | COSMIC Variant | Genomic variant prefix |
| CGC | Cancer Gene Census | Curated gene list |
| CNV | Copy Number Variation | Amplification/deletion |
| FATHMM | Functional Analysis Through Hidden Markov Models | Prediction tool |
| CDS | Coding DNA Sequence | DNA notation |
| AA | Amino Acid | Protein notation |
| HGVS | Human Genome Variation Society | Nomenclature |
| TCGA | The Cancer Genome Atlas | Data source |
| ICGC | International Cancer Genome Consortium | Data source |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| ClinVar | VCV | Clinical significance |
| HGNC | Symbol | Gene naming |
| TCGA | Sample ID | Sample source |
| ICGC | Donor ID | Sample source |
| PubMed | PMID | Publication |

---

## Data Subsets

| Dataset | Description |
|---------|-------------|
| CosmicMutantExport | All coding mutations |
| CosmicNCV | Non-coding variants |
| CosmicCompleteCNA | Copy number data |
| CosmicFusionExport | Gene fusions |
| CosmicResistanceMutations | Drug resistance |
| CancerGeneCensus | Curated cancer genes |
| CosmicMutantExportCensus | CGC gene mutations |

---

## Data Quality Notes

1. **Cardinality:** One entry per mutation-sample observation
2. **Curation:** Expert-curated from literature + large studies
3. **Sources:** TCGA, ICGC, published studies, cell lines
4. **Versioning:** Regular releases (v99+)
5. **Access:** Requires license for commercial use
6. **ID Migration:** COSM→COSV ID system transition
