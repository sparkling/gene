# 3.2 Phenotype Databases - Data Dictionary

## Overview

This subcategory contains phenotype annotation data linking diseases to observable clinical characteristics, including symptom frequencies and evidence sources.

**Data Sources:** HPO Annotations, OMIM, Orphanet Phenotypes

---

## Unified Fields

These fields are harmonized across multiple data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `disease_id` | string | Required (1:1) | Disease database identifier | HPO Annotations, OMIM, Orphanet Phenotypes | `OMIM:154700`, `ORPHA:558`, `MONDO:0007947` |
| `disease_name` | string | Required (1:1) | Name of the disease | HPO Annotations, OMIM, Orphanet Phenotypes | `Marfan syndrome`, `Huntington disease` |
| `hpo_ids` | array[string] | Optional (N:M) | HPO phenotype term identifiers. Pattern: `HP:[0-9]{7}` | HPO Annotations, Orphanet Phenotypes | `["HP:0001166", "HP:0001519", "HP:0002650"]` |
| `frequency` | string | Optional (1:1) | Frequency of phenotype occurrence in affected individuals | HPO Annotations, Orphanet Phenotypes | `HP:0040281 (Very frequent)`, `12/45`, `80-99%` |
| `gene_id` | integer | Optional (N:M) | NCBI Gene identifier for associated genes | HPO Annotations, OMIM | `2200`, `3077` |
| `gene_symbol` | string | Optional (N:M) | HGNC gene symbol | HPO Annotations, OMIM | `FBN1`, `HTT`, `BRCA1` |

---

## Source-Specific Fields

### HPO Annotations

| Field Name | Data Type | Cardinality | Description | Allowed Values | Example Values |
|------------|-----------|-------------|-------------|----------------|----------------|
| `qualifier` | string | Optional (1:1) | NOT if phenotype is excluded from disease | `NOT`, null | `NOT` |
| `reference` | string | Optional (1:1) | Evidence source (PMID, OMIM number) | - | `PMID:12345678`, `OMIM:154700` |
| `evidence` | enum | Optional | Evidence code for annotation | `IEA`, `TAS`, `PCS` | `TAS` |
| `onset` | string | Optional (1:1) | Age of onset (HPO term) | - | `HP:0003577` (Congenital onset) |
| `sex` | enum | Optional | Sex specificity of phenotype | `MALE`, `FEMALE` | `MALE` |
| `modifier` | array[string] | Optional (1:N) | Clinical modifiers (HP terms) | - | `["HP:0012825", "HP:0012828"]` |
| `aspect` | enum | Optional | HPO branch | `P` (Phenotype), `I` (Inheritance), `C` (Course), `M` (Modifier), `H` (History) | `P` |
| `biocuration` | string | Optional (1:1) | Curator ID and date stamp | - | `HPO:skoehler[2020-01-15]` |

---

### OMIM

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `mim_number` | string | Optional | MIM number with prefix indicating entry type | `[*#%+]?[0-9]{6}` | `#154700`, `*134797`, `%231200` |
| `entry_type` | enum | Optional | Type of OMIM entry based on prefix | `gene (*)`, `phenotype (#)`, `combined (+)`, `unknown (%)` | `phenotype (#)` |
| `allelic_variants` | array[string] | Optional (1:N) | Known pathogenic variants | - | `["FBN1, CYS1039TYR", "FBN1, ARG2726TRP"]` |

---

### Orphanet Phenotypes

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `hpo_frequency` | object | Optional | Structured frequency classification | See frequency classification table below |
| `diagnostic` | boolean | Optional | Whether phenotype is a diagnostic criterion | `true`, `false` |

#### HPO Frequency Classification

| HPO ID | Label | Percentage |
|--------|-------|------------|
| `HP:0040280` | Obligate | 100% |
| `HP:0040281` | Very frequent | 80-99% |
| `HP:0040282` | Frequent | 30-79% |
| `HP:0040283` | Occasional | 5-29% |
| `HP:0040284` | Very rare | 1-4% |
| `HP:0040285` | Excluded | 0% |

---

## Source Field Mappings

### HPO Annotations Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `DatabaseID` | `disease_id` |
| `DiseaseName` | `disease_name` |
| `HPO_ID` | `hpo_ids` |
| `Frequency` | `frequency` |
| `Gene_ID` | `gene_id` |
| `Gene_Symbol` | `gene_symbol` |
| `Qualifier` | `qualifier` |
| `Reference` | `reference` |
| `Evidence` | `evidence` |
| `Onset` | `onset` |
| `Sex` | `sex` |
| `Modifier` | `modifier` |
| `Aspect` | `aspect` |
| `Biocuration` | `biocuration` |

### OMIM Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `MIM_Number` | `mim_number` |
| `Title` | `disease_name` |
| `MIM_Entry_Type` | `entry_type` |
| `Gene_Symbols` | `gene_symbol` |
| `Allelic_Variants` | `allelic_variants` |

### Orphanet Phenotypes Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `OrphaCode` | `disease_id` |
| `Name` | `disease_name` |
| `HPO_ID` | `hpo_ids` |
| `HPO_Frequency` | `hpo_frequency` |
| `Diagnostic` | `diagnostic` |

---

## Evidence Codes

| Code | Full Name | Description |
|------|-----------|-------------|
| `IEA` | Inferred from Electronic Annotation | Computationally assigned annotation |
| `TAS` | Traceable Author Statement | Curated from published literature |
| `PCS` | Published Clinical Study | Based on clinical study data |

---

## Metadata Fields

| Field Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `_source.database` | string | Name of the source database | `HPO_Annotations`, `OMIM`, `Orphanet_Phenotypes` |
| `_source.version` | string | Version of the source data | `2024-01-15`, `v2024.1` |
| `_source.access_date` | date | Date the data was accessed | `2026-01-24` |
| `_source.original_id` | string | Original identifier in source | `OMIM:154700` |
