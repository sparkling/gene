# 3.5 Rare/Orphan Diseases - Data Dictionary

## Overview

This subcategory contains data about rare and orphan diseases, including clinical features, genetic associations, and diagnostic gene panels.

**Data Sources:** DECIPHER, Orphanet, PanelApp

---

## Unified Fields

These fields are harmonized across multiple data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `disease_id` | string | Required (1:1) | Rare disease identifier | DECIPHER, Orphanet, PanelApp | `Orphanet:558`, `DECIPHER:12345` |
| `disease_name` | string | Required (1:1) | Name of rare disease | DECIPHER, Orphanet, PanelApp | `Marfan syndrome`, `Cystic fibrosis` |
| `gene_symbols` | array[string] | Optional (N:M) | Associated HGNC gene symbols | DECIPHER, Orphanet, PanelApp | `["FBN1", "TGFBR2"]` |
| `hpo_ids` | array[string] | Optional (N:M) | HPO phenotype identifiers | DECIPHER, Orphanet | `["HP:0001166", "HP:0001519"]` |
| `inheritance_patterns` | array[string] | Optional (1:N) | Modes of inheritance | Orphanet, PanelApp | `["autosomal dominant", "autosomal recessive", "X-linked", "mitochondrial"]` |

---

## Source-Specific Fields

### DECIPHER

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `decipher_id` | integer | Optional | DECIPHER patient/record identifier | - | `12345`, `67890` |
| `cnv_type` | string | Optional | Copy number variant type | - | `gain`, `loss` |
| `hi_score` | float | Optional | Haploinsufficiency prediction score | 0-1 | `0.85`, `0.12` |
| `pli` | float | Optional | Probability of loss-of-function intolerance | 0-1 | `0.99`, `0.45` |
| `loeuf` | float | Optional | LoF observed/expected upper bound fraction | 0-2 | `0.15`, `0.89` |
| `pathogenicity` | string | Optional | Variant classification | - | `Pathogenic`, `Likely pathogenic`, `Uncertain significance` |
| `ensembl_gene` | string | Optional | Ensembl gene identifier | `ENSG[0-9]{11}` | `ENSG00000166147` |

---

### Orphanet

#### Disease Classification

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `orpha_code` | integer | Optional | Orphanet disease code | `558`, `79383` |
| `disorder_type` | object | Optional | Classification (disease, malformation, group) | `{"name": "Disease", "id": "21394"}` |

#### Epidemiology

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `average_age_of_onset` | array[string] | Optional | Age categories for disease onset | `["Infancy", "Childhood", "Adolescence"]` |
| `prevalence` | string | Optional | Disease prevalence estimate | `1-9 / 100 000`, `Unknown` |
| `prevalence_class` | string | Optional | Categorical prevalence classification | `Rare`, `Very rare`, `Ultra-rare` |
| `prevalence_geographic` | string | Optional | Geographic scope of prevalence | `Worldwide`, `Europe` |

#### Gene Associations

| Field Name | Data Type | Cardinality | Description | Allowed Values | Example Values |
|------------|-----------|-------------|-------------|----------------|----------------|
| `gene_association_type` | string | Optional | Type of gene-disease relationship | - | `Disease-causing germline mutation(s) in`, `Major susceptibility factor in` |
| `gene_association_status` | enum | Optional | Curation status of gene link | `Assessed`, `Not validated` | `Assessed` |
| `hgnc_id` | string | Optional | HGNC gene identifier | - | `HGNC:3603` |
| `gene_locus` | string | Optional | Chromosomal location | - | `15q21.1`, `7q31.2` |

#### Cross-References

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `external_references` | array[object] | Optional | Cross-references to OMIM, ICD, UMLS, MeSH | `[{"source": "OMIM", "id": "154700"}]` |

---

### PanelApp

| Field Name | Data Type | Cardinality | Description | Range/Allowed Values | Example Values |
|------------|-----------|-------------|-------------|----------------------|----------------|
| `panel_id` | integer | Optional | Gene panel identifier | - | `123`, `456` |
| `panel_name` | string | Optional | Name of gene panel | - | `Intellectual disability`, `Cardiomyopathy` |
| `evidence_level` | integer | Optional | Gene evidence level | 0-3 | `3` (high evidence) |
| `confidence_color` | enum | Optional | Visual confidence indicator | `Green`, `Amber`, `Red`, `Gray` | `Green` |
| `hgnc_id` | string | Optional | HGNC gene identifier | `HGNC:[0-9]+` | `HGNC:3603` |
| `expert_reviews` | integer | Optional | Number of expert reviews | - | `5`, `12` |

---

## Source Field Mappings

### DECIPHER Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `decipher_id` | `decipher_id` |
| `syndrome_name` | `disease_name` |
| `gene_symbol` | `gene_symbols` |
| `hpo_id` | `hpo_ids` |
| `cnv_type` | `cnv_type` |
| `hi_score` | `hi_score` |
| `pLI` | `pli` |
| `LOEUF` | `loeuf` |
| `pathogenicity` | `pathogenicity` |
| `ensembl_gene` | `ensembl_gene` |

### Orphanet Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `OrphaCode` | `orpha_code` |
| `Name` | `disease_name` |
| `Gene_Symbol` | `gene_symbols` |
| `HPO_ID` | `hpo_ids` |
| `Inheritance` | `inheritance_patterns` |
| `DisorderType` | `disorder_type` |
| `AverageAgeOfOnset` | `average_age_of_onset` |
| `Prevalence` | `prevalence` |
| `PrevalenceClass` | `prevalence_class` |
| `PrevalenceGeographic` | `prevalence_geographic` |
| `GeneAssociationType` | `gene_association_type` |
| `GeneAssociationStatus` | `gene_association_status` |
| `HGNC_ID` | `hgnc_id` |
| `GeneType` | `gene_locus` |
| `ExternalReference` | `external_references` |

### PanelApp Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `panel_id` | `panel_id` |
| `panel_name` | `panel_name` |
| `gene_symbol` | `gene_symbols` |
| `hgnc_id` | `hgnc_id` |
| `evidence_level` | `evidence_level` |
| `confidence_level` | `confidence_color` |
| `mode_of_inheritance` | `inheritance_patterns` |
| `expert_reviews` | `expert_reviews` |

---

## PanelApp Evidence Levels

| Level | Color | Description |
|-------|-------|-------------|
| 3 | Green | High evidence - diagnostic grade |
| 2 | Amber | Moderate evidence - supporting evidence |
| 1 | Red | Limited evidence - only limited evidence |
| 0 | Gray | No evidence - no current evidence |

---

## Inheritance Patterns

| Pattern | Description |
|---------|-------------|
| `autosomal dominant` | One copy of altered gene sufficient to cause condition |
| `autosomal recessive` | Two copies of altered gene required to cause condition |
| `X-linked dominant` | One copy on X chromosome sufficient; affects males and females |
| `X-linked recessive` | Usually affects males; females are carriers |
| `mitochondrial` | Inherited through maternal mitochondrial DNA |
| `de novo` | New mutation not inherited from parents |
| `digenic` | Requires variants in two different genes |

---

## Constraint Scores

| Score | Description | Interpretation |
|-------|-------------|----------------|
| `pLI` | Probability of Loss-of-function Intolerance | Higher values (>0.9) indicate gene is intolerant to LoF |
| `LOEUF` | Loss-of-function Observed/Expected Upper Fraction | Lower values (<0.35) indicate LoF intolerance |
| `hi_score` | Haploinsufficiency Score | Higher values indicate likely haploinsufficient |

---

## Metadata Fields

| Field Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `_source.database` | string | Name of the source database | `DECIPHER`, `Orphanet`, `PanelApp` |
| `_source.version` | string | Version of the source data | `v11.0`, `2024-01`, `1.28` |
| `_source.access_date` | date | Date the data was accessed | `2026-01-24` |
| `_source.original_id` | string | Original identifier in source | `Orphanet:558` |
