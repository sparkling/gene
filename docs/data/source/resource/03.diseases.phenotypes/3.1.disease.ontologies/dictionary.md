# 3.1 Disease Ontologies - Data Dictionary

## Overview

This subcategory contains disease ontology data from structured vocabularies and classification systems used to standardize disease terminology across biomedical research.

**Data Sources:** EFO, HPO, ICD, MeSH, MONDO, Orphanet/ORDO

---

## Unified Fields

These fields are harmonized across multiple data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `id` | string | Required (1:1) | Primary unique identifier for the disease/ontology term | EFO, HPO, ICD, MeSH, MONDO, Orphanet | `EFO_0000685`, `HP:0001250`, `E11.9`, `5A11`, `D003920`, `MONDO:0005015`, `Orphanet:558` |
| `name` | string | Required (1:1) | Primary/preferred name of the disease or ontology term | EFO, HPO, ICD, MeSH, MONDO, Orphanet | `Diabetes mellitus`, `Seizure`, `Marfan syndrome` |
| `definition` | string | Optional (1:1) | Textual definition of the disease or phenotype. Often includes citation references. | HPO, MONDO, Orphanet | `A disorder characterized by abnormally high blood glucose levels...` |
| `synonyms` | array[string] | Optional (1:N) | Alternative names or aliases for the term. May include type qualifiers (exact, broad, narrow, related). | EFO, HPO, MeSH, MONDO, Orphanet | `["Type 2 diabetes", "Adult-onset diabetes", "NIDDM"]` |
| `is_a` | array[string] | Optional (N:M) | Parent term relationships in ontology hierarchy (subclass relationship). | EFO, HPO, MONDO, Orphanet/ORDO | `["HP:0012638", "MONDO:0005066"]` |
| `xref` | array[string] | Optional (1:N) | Cross-references to external databases and ontologies (UMLS, SNOMED, OMIM, etc.). | EFO, HPO, ICD, MeSH, MONDO, Orphanet | `["UMLS:C0011847", "SNOMED:73211009", "OMIM:125853"]` |

---

## Source-Specific Fields

### EFO (Experimental Factor Ontology)

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `efo_id` | string | Optional | Experimental Factor Ontology identifier | `EFO_[0-9]{7}` | `EFO_0000685` |

---

### HPO (Human Phenotype Ontology)

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `hpo_comment` | string | Optional (1:1) | Additional notes about the phenotype term | - | `This term should be used when...` |
| `created_by` | string | Optional (1:1) | Curator who created the term | - | `peter`, `sebastian` |
| `creation_date` | datetime | Optional (1:1) | Date the term was created | ISO 8601 | `2012-03-15T10:30:00Z` |
| `subset` | array[string] | Optional (1:N) | Subsets the term belongs to | - | `["hpo_slim", "hposlim_core"]` |
| `alt_id` | array[string] | Optional (1:N) | Alternative/obsolete identifiers for the term | - | `["HP:0000001", "HP:0000002"]` |

---

### ICD (International Classification of Diseases)

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `icd10_code` | string | Optional | ICD-10 classification code | `[A-Z][0-9]{2}(\.[0-9]{1,4})?` | `E11.9`, `J45.0`, `G40.309` |
| `icd11_code` | string | Optional | ICD-11 entity code | `[A-Z0-9]{4,6}` | `5A11`, `CA40` |
| `block_code` | string | Optional | ICD chapter block range | `[A-Z][0-9]{2}-[A-Z][0-9]{2}` | `E10-E14`, `J40-J47` |

---

### MeSH (Medical Subject Headings)

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `mesh_ui` | string | Optional | MeSH unique identifier | `D[0-9]{6}` | `D003920`, `D001249` |
| `tree_number` | string | Optional | MeSH tree hierarchy position | `C[0-9]{2}(\.[0-9]{3})+` | `C18.452.394.750`, `C08.381.495` |
| `entry_term` | array[string] | Optional (1:N) | Entry terms for indexing | - | `["Diabetes, Type 2", "Maturity-Onset Diabetes"]` |

---

### MONDO (Mondo Disease Ontology)

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `exact_match` | array[string] | Optional (1:N) | 1:1 equivalence axioms to other disease IDs | `["OMIM:125853", "DOID:9352"]` |
| `close_match` | array[string] | Optional (1:N) | Very similar but not exact matches | `["HP:0000819"]` |
| `narrow_match` | array[string] | Optional (1:N) | More specific concept mappings | `["Orphanet:79383"]` |
| `broad_match` | array[string] | Optional (1:N) | More general concept mappings | `["MONDO:0005015"]` |

---

### Orphanet

| Field Name | Data Type | Cardinality | Description | Allowed Values | Example Values |
|------------|-----------|-------------|-------------|----------------|----------------|
| `orpha_code` | integer | Optional | Numeric Orphanet disease code | - | `558`, `79383`, `166024` |
| `orpha_group` | enum | Optional | Classification of entity type | `Disorder`, `Group`, `Subtype` | `Disorder` |
| `disorder_type` | object | Optional | Type of disorder (disease, malformation, etc.) | - | `{"name": "Disease", "id": "21394"}` |
| `association_type` | string | Optional | Type of gene-disease relationship | `Disease-causing germline mutation(s) in`, `Disease-causing somatic mutation(s) in`, `Major susceptibility factor in`, `Modifying germline mutation in` | `Disease-causing germline mutation(s) in` |
| `association_status` | enum | Optional | Curation status of gene-disease association | `Assessed`, `Not validated` | `Assessed` |

---

## Source Field Mappings

### EFO Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `efo_id` | `id` |
| `label` | `name` |
| `definition` | `definition` |
| `synonyms` | `synonyms` |
| `is_a` | `is_a` |
| `xref` | `xref` |

### HPO Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `id` | `id` |
| `name` | `name` |
| `def` | `definition` |
| `synonym` | `synonyms` |
| `is_a` | `is_a` |
| `xref` | `xref` |
| `comment` | `hpo_comment` |
| `created_by` | `created_by` |
| `creation_date` | `creation_date` |
| `subset` | `subset` |
| `alt_id` | `alt_id` |

### ICD Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `code` | `id` |
| `title` | `name` |
| `icd10_code` | `icd10_code` |
| `icd11_code` | `icd11_code` |
| `block_code` | `block_code` |

### MeSH Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `UI` | `mesh_ui` |
| `MH` | `name` |
| `MS` | `definition` |
| `SY` | `synonyms` |
| `MN` | `tree_number` |
| `ENTRY` | `entry_term` |

### MONDO Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `id` | `id` |
| `label` | `name` |
| `definition` | `definition` |
| `synonyms` | `synonyms` |
| `is_a` | `is_a` |
| `xref` | `xref` |
| `exactMatch` | `exact_match` |
| `closeMatch` | `close_match` |
| `narrowMatch` | `narrow_match` |
| `broadMatch` | `broad_match` |

### Orphanet Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `OrphaCode` | `orpha_code` |
| `Name` | `name` |
| `Definition` | `definition` |
| `Synonym` | `synonyms` |
| `OrphaGroup` | `orpha_group` |
| `DisorderType` | `disorder_type` |
| `AssociationType` | `association_type` |
| `AssociationStatus` | `association_status` |
| `ExternalReference` | `xref` |

---

## Metadata Fields

| Field Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `_source.database` | string | Name of the source database | `EFO`, `HPO`, `MONDO` |
| `_source.version` | string | Version of the source data | `2024-01`, `v3.2.0` |
| `_source.access_date` | date | Date the data was accessed | `2026-01-24` |
| `_source.original_id` | string | Original identifier in source | `EFO_0000685` |
