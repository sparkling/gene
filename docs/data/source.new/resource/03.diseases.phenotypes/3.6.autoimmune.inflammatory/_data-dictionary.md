# 3.6 Autoimmune/Inflammatory Diseases - Data Dictionary

## Overview

This subcategory contains data about autoimmune and inflammatory diseases, including GWAS associations and HLA allele information relevant to disease susceptibility.

**Data Sources:** ImmunoBase, IPD-IMGT/HLA

---

## Unified Fields

These fields are harmonized across multiple data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `gene_symbol` | string | Required (1:1) | Gene symbol (HGNC or HLA locus) | ImmunoBase, IPD-IMGT/HLA | `IL23R`, `HLA-DRB1`, `PTPN22` |
| `disease_associations` | array[string] | Optional (N:M) | Associated autoimmune/inflammatory diseases | ImmunoBase, IPD-IMGT/HLA | `["Type 1 diabetes", "Rheumatoid arthritis"]` |
| `snp_ids` | array[string] | Optional (1:N) | dbSNP variant identifiers. Pattern: `rs[0-9]+` | ImmunoBase | `["rs2476601", "rs3129889"]` |

---

## Source-Specific Fields

### ImmunoBase

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `disease_code` | string | Optional | Disease abbreviation | `T1D`, `RA`, `CEL`, `MS`, `CD`, `UC` |
| `region_id` | string | Optional | Chromosomal region identifier | `1p13.2`, `6p21.32`, `2q33.2` |
| `credible_set` | array[string] | Optional | Fine-mapped credible set of variants | `["rs2476601", "rs6679677"]` |
| `gwas_loci` | integer | Optional | Number of associated GWAS loci | `42`, `123` |

---

### IPD-IMGT/HLA

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `hla_allele` | string | Optional | Full HLA allele designation | `HLA-[A-Z]+[0-9]*\*[0-9:]+` | `HLA-A*02:01:01:01`, `HLA-DRB1*04:01` |
| `allele_group` | string | Optional | Allele group (serological equivalent) | `HLA-[A-Z]+\*[0-9]+` | `HLA-A*02`, `HLA-DRB1*04` |
| `ipd_accession` | string | Optional | IPD database accession | `HLA[0-9]+` | `HLA00001`, `HLA14859` |
| `gene_locus` | string | Optional | HLA gene locus | `HLA-[A-Z]+[0-9]*` | `HLA-A`, `HLA-DRB1`, `HLA-B` |
| `sequence_type` | string | Optional | Type of sequence data | - | `genomic`, `cDNA`, `protein` |

---

## Source Field Mappings

### ImmunoBase Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `gene_symbol` | `gene_symbol` |
| `disease` | `disease_associations` |
| `disease_code` | `disease_code` |
| `snp_id` | `snp_ids` |
| `region_id` | `region_id` |
| `credible_set` | `credible_set` |
| `gwas_loci` | `gwas_loci` |

### IPD-IMGT/HLA Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `locus` | `gene_symbol` |
| `allele` | `hla_allele` |
| `allele_group` | `allele_group` |
| `accession` | `ipd_accession` |
| `gene_locus` | `gene_locus` |
| `sequence_type` | `sequence_type` |
| `disease_association` | `disease_associations` |

---

## Disease Codes (ImmunoBase)

| Code | Full Name |
|------|-----------|
| `T1D` | Type 1 Diabetes |
| `RA` | Rheumatoid Arthritis |
| `CEL` | Celiac Disease |
| `MS` | Multiple Sclerosis |
| `CD` | Crohn's Disease |
| `UC` | Ulcerative Colitis |
| `PSO` | Psoriasis |
| `AS` | Ankylosing Spondylitis |
| `SLE` | Systemic Lupus Erythematosus |
| `JIA` | Juvenile Idiopathic Arthritis |
| `ATD` | Autoimmune Thyroid Disease |
| `VIT` | Vitiligo |
| `AA` | Alopecia Areata |
| `PBC` | Primary Biliary Cholangitis |
| `IBD` | Inflammatory Bowel Disease |

---

## HLA Nomenclature

### Allele Naming Convention

HLA alleles follow a standardized naming system:

```
HLA-[Locus]*[Allele Group]:[Protein]:[Synonymous]:[Non-coding]
```

| Component | Description | Example |
|-----------|-------------|---------|
| Locus | Gene name | `A`, `B`, `DRB1` |
| Allele Group | Major type (2 digits) | `02` |
| Protein | Protein sequence differences | `01` |
| Synonymous | Coding region synonymous differences | `01` |
| Non-coding | Intronic/UTR differences | `01` |

### Example Allele Resolution Levels

| Resolution | Example | Description |
|------------|---------|-------------|
| Low (2-digit) | `HLA-A*02` | Allele group only |
| Intermediate | `HLA-A*02:01` | Allele group + protein |
| High | `HLA-A*02:01:01` | Full coding sequence |
| Highest | `HLA-A*02:01:01:01` | Full genomic sequence |

---

## HLA Gene Classes

| Class | Loci | Function |
|-------|------|----------|
| Class I | HLA-A, HLA-B, HLA-C | Present intracellular peptides to CD8+ T cells |
| Class II | HLA-DR, HLA-DQ, HLA-DP | Present extracellular peptides to CD4+ T cells |
| Class III | Complement, TNF | Various immune functions |

---

## Common HLA-Disease Associations

| HLA Allele | Disease | Risk |
|------------|---------|------|
| `HLA-B*27` | Ankylosing spondylitis | Strong positive |
| `HLA-DRB1*04` | Rheumatoid arthritis | Positive |
| `HLA-DQ2/DQ8` | Celiac disease | Required |
| `HLA-DRB1*15:01` | Multiple sclerosis | Positive |
| `HLA-DRB1*03/04` | Type 1 diabetes | Positive |
| `HLA-B*57:01` | Abacavir hypersensitivity | Strong positive |

---

## Metadata Fields

| Field Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `_source.database` | string | Name of the source database | `ImmunoBase`, `IPD_IMGT_HLA` |
| `_source.version` | string | Version of the source data | `2024-01`, `3.55.0` |
| `_source.access_date` | date | Date the data was accessed | `2026-01-24` |
| `_source.original_id` | string | Original identifier in source | `HLA00001` |
