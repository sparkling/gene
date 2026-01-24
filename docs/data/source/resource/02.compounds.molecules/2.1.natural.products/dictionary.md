# 2.1 Natural Products - Data Dictionary

## Overview

This data dictionary documents the unified schema for natural product compound data from five major databases: COCONUT, Dr. Duke's Phytochemical and Ethnobotanical Databases, LOTUS, NPASS, and NPAtlas.

**Subcategory ID:** 2.1
**Subcategory Name:** Natural Products
**Data Sources:** COCONUT, Dr. Duke's, LOTUS, NPASS, NPAtlas

---

## Unified Fields

### Core Identification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| compound_id | string | 1:1, Required | Primary identifier for natural product compound | `CNP0123456`, `NPC12345`, `NPA012345` | COCONUT: coconut_id, NPASS: npc_id, NPAtlas: npa_id |
| name | string | 1:1, Required | Common or trivial name of the natural product | `Artemisinin`, `Curcumin`, `Quercetin` | COCONUT, Dr. Duke's, LOTUS, NPASS, NPAtlas |

### Structure Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| canonical_smiles | string | 1:1, Optional | Simplified Molecular Input Line Entry System notation (no stereochemistry) | `CC(C)CCCC(C)C`, `CC(=O)Oc1ccccc1C(=O)O` | COCONUT, LOTUS, NPASS, NPAtlas |
| isomeric_smiles | string | 1:1, Optional | SMILES with stereochemistry preserved | `C[C@@H](O)c1ccccc1` | COCONUT, LOTUS |
| inchi | string | 1:1, Optional | IUPAC International Chemical Identifier | `InChI=1S/C6H12O6/c7-1-2-3(8)...` | COCONUT, LOTUS, NPASS, NPAtlas |
| inchi_key | string | 1:1, Optional | 27-character hashed InChI for fast lookups | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` | COCONUT, LOTUS, NPASS, NPAtlas |

### Molecular Properties

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| molecular_formula | string | 1:1, Optional | Chemical formula showing atom composition | `C21H30O2`, `C9H8O4` | COCONUT, Dr. Duke's, LOTUS, NPASS, NPAtlas |
| molecular_weight | decimal | 1:1, Optional | Molecular weight in Daltons | `282.35`, `180.16` | COCONUT, Dr. Duke's, LOTUS, NPASS, NPAtlas |

### Biological Source Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| organism_sources | object[] | 1:N, Optional | Biological organisms producing the natural product | `[{"name": "Artemisia annua", "ncbi_taxon_id": 35608}]` | COCONUT, Dr. Duke's, LOTUS, NPASS, NPAtlas |
| ncbi_taxon_id | integer | 1:1, Optional | NCBI Taxonomy identifier for source organism | `35608`, `3702` | COCONUT, LOTUS, NPASS, NPAtlas |

### Classification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| compound_classes | string[] | 1:N, Optional | Natural product classification | `Terpenoid`, `Alkaloid`, `Flavonoid`, `Polyketide` | COCONUT, NPASS, NPAtlas |

### Reference Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| references | object[] | 1:N, Optional | Literature citations for compound isolation/characterization | `[{"doi": "10.1016/j.phytochem.2018.01.001", "pmid": 29371045}]` | LOTUS, NPASS, NPAtlas |

---

## Source-Specific Fields

### COCONUT

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| coconut_id | string | 1:1, Required | COCONUT identifier (CNP prefix) | `CNP0123456` |
| sugar_free_smiles | string | 1:1, Optional | SMILES with glycosidic moieties removed | - |
| qed_drug_likeness | decimal | 1:1, Optional | Quantitative Estimate of Drug-likeness (0-1) | `0.65` |
| lipinski_violations | integer | 1:1, Optional | Number of Lipinski rule violations (0-4) | `0`, `1`, `2` |
| murko_framework | string | 1:1, Optional | Core molecular scaffold structure | - |

### Dr. Duke's Phytochemical Database

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| cas_number | string | 1:1, Optional | CAS Registry Number | `458-37-7` |
| biological_activities | string[] | 1:N, Optional | Documented biological activities (1,900+ types) | `Antimicrobial`, `Anti-inflammatory` |
| ethnobotanical_uses | object[] | 1:N, Optional | Traditional/ethnobotanical uses with culture and region | `[{"use": "fever treatment", "culture": "Chinese", "region": "Asia"}]` |
| plant_part | string | 1:1, Optional | Part of plant containing compound | `Leaf`, `Root`, `Bark` |
| concentration_range | object | 1:1, Optional | Concentration in plant (low/high values in ppm) | `{"low": 100, "high": 5000, "unit": "ppm"}` |

### LOTUS

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| wikidata_id | string | 1:1, Optional | Wikidata Q-identifier | `Q312879` |
| npl_score | decimal | 1:1, Optional | Natural product likeness score | `0.85` |
| deep_smiles | string | 1:1, Optional | DeepSMILES representation | - |

### NPASS

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| npc_id | string | 1:1, Required | NPASS compound ID (NPC prefix) | `NPC12345` |
| activity_value | decimal | 1:1, Optional | Quantitative bioactivity (IC50, Ki, etc.) | `0.5` |
| activity_type | string | 1:1, Optional | Type of bioactivity measurement | `IC50`, `Ki`, `EC50` |
| target_protein | object | 1:1, Optional | Biological target with UniProt reference | `{"name": "COX-2", "uniprot_id": "P35354"}` |

### NPAtlas

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| npa_id | string | 1:1, Required | NPAtlas identifier (NPA prefix) | `NPA012345` |
| origin_type | string | 1:1, Optional | Source organism type | `bacterial`, `fungal`, `marine` |
| cluster_type | string | 1:1, Optional | Biosynthetic gene cluster type | `NRPS`, `PKS`, `Terpene` |
| isolation_source | string | 1:1, Optional | Environment of isolation | `marine`, `soil`, `endophytic` |

---

## Metadata Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| _source.database | string | 1:1, Required | Source database name | `COCONUT`, `LOTUS`, `NPASS` |
| _source.version | string | 1:1, Optional | Database version | `2024.01` |
| _source.access_date | date | 1:1, Optional | Date data was accessed | `2026-01-24` |
| _source.original_id | string | 1:1, Optional | Original identifier from source | - |

---

## Field Mappings by Source

### COCONUT to Unified Schema

| COCONUT Field | Unified Field |
|---------------|---------------|
| coconut_id | coconut_id |
| name | name |
| canonical_smiles | canonical_smiles |
| isomeric_smiles | isomeric_smiles |
| inchi | inchi |
| inchikey | inchi_key |
| molecular_formula | molecular_formula |
| molecular_weight | molecular_weight |
| organism_source | organism_sources |
| ncbi_taxon_id | ncbi_taxon_id |
| compound_class | compound_classes |
| sugar_free_smiles | sugar_free_smiles |
| qed_drug_likeliness | qed_drug_likeness |
| lipinski_rule_of_5 | lipinski_violations |
| murko_framework | murko_framework |

### Dr. Duke's to Unified Schema

| Dr. Duke's Field | Unified Field |
|------------------|---------------|
| common_name | name |
| cas_number | cas_number |
| biological_activity | biological_activities |
| ethnobotanical_use | ethnobotanical_uses |
| plant_part | plant_part |
| concentration_range | concentration_range |

### LOTUS to Unified Schema

| LOTUS Field | Unified Field |
|-------------|---------------|
| name | name |
| canonical_smiles | canonical_smiles |
| isomeric_smiles | isomeric_smiles |
| inchi | inchi |
| inchi_key | inchi_key |
| molecular_formula | molecular_formula |
| molecular_weight | molecular_weight |
| organism_source | organism_sources |
| ncbi_taxon_id | ncbi_taxon_id |
| wikidataId | wikidata_id |
| npl_score | npl_score |
| deep_smiles | deep_smiles |

### NPASS to Unified Schema

| NPASS Field | Unified Field |
|-------------|---------------|
| npc_id | npc_id |
| name | name |
| canonical_smiles | canonical_smiles |
| inchi | inchi |
| inchi_key | inchi_key |
| molecular_formula | molecular_formula |
| molecular_weight | molecular_weight |
| organism_source | organism_sources |
| ncbi_taxon_id | ncbi_taxon_id |
| compound_class | compound_classes |
| activity_value | activity_value |
| activity_type | activity_type |
| target_protein | target_protein |
| references | references |

### NPAtlas to Unified Schema

| NPAtlas Field | Unified Field |
|---------------|---------------|
| npa_id | npa_id |
| name | name |
| canonical_smiles | canonical_smiles |
| inchi | inchi |
| inchi_key | inchi_key |
| molecular_formula | molecular_formula |
| molecular_weight | molecular_weight |
| organism_source | organism_sources |
| ncbi_taxon_id | ncbi_taxon_id |
| compound_class | compound_classes |
| origin_type | origin_type |
| cluster_type | cluster_type |
| isolation_source | isolation_source |
| references | references |

---

## Data Quality Notes

- **Required Field:** Only `name` is strictly required across all sources
- **Structure Data:** Chemical structure fields (SMILES, InChI) may be null for some records
- **Taxonomy:** NCBI taxon IDs are available from most sources except Dr. Duke's
- **Bioactivity:** Quantitative bioactivity data primarily from NPASS
- **Ethnobotanical Data:** Unique to Dr. Duke's database
