# 2.2 Pharmaceuticals - Data Dictionary

## Overview

This data dictionary documents the unified schema for pharmaceutical data from five major databases: ChEMBL, DailyMed, DrugBank, FDA Orange Book, and RxNorm.

**Subcategory ID:** 2.2
**Subcategory Name:** Pharmaceuticals
**Data Sources:** ChEMBL, DailyMed, DrugBank, Orange Book, RxNorm

---

## Unified Fields

### Core Identification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| drug_id | string | 1:1, Required | Primary identifier for drug compound | `CHEMBL25`, `DB00945`, `2670` | ChEMBL: chembl_id, DrugBank: drugbank_id, RxNorm: rxcui |
| name | string | 1:1, Required | Drug name (generic or brand) | `Aspirin`, `Acetaminophen` | ChEMBL, DailyMed, DrugBank, Orange Book, RxNorm |

### Structure Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| canonical_smiles | string | 1:1, Optional | Canonical SMILES notation | `CC(=O)Oc1ccccc1C(=O)O` | ChEMBL, DrugBank |
| inchi_key | string | 1:1, Optional | InChI Key identifier | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` | ChEMBL, DrugBank |

### Molecular Properties

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| molecular_weight | decimal | 1:1, Optional | Molecular weight in Daltons | `180.16` | ChEMBL, DrugBank |
| molecular_formula | string | 1:1, Optional | Chemical formula | `C9H8O4` | ChEMBL, DrugBank |

### Classification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| drug_type | string | 1:1, Optional | Type classification of drug | `small molecule`, `biotech`, `antibody` | ChEMBL, DrugBank |
| approval_status | string | 1:1, Optional | Regulatory approval status | `approved`, `investigational`, `Phase III` | ChEMBL, DrugBank, Orange Book |

### Target Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| targets | object[] | 1:N, Optional | Drug targets with actions | `[{"name": "Cyclooxygenase-1", "uniprot_id": "P23219", "action": "inhibitor"}]` | ChEMBL, DrugBank |
| uniprot_id | string | 1:1, Optional | UniProt accession for target protein | `P23219`, `P35354` | ChEMBL, DrugBank |

---

## Source-Specific Fields

### ChEMBL

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| chembl_id | string | 1:1, Required | ChEMBL identifier | `CHEMBL25` |
| molregno | integer | 1:1, Optional | Internal molecule registration number | `123456` |
| max_phase | integer | 1:1, Optional | Highest clinical phase (0-4) | `0`, `1`, `2`, `3`, `4` |
| pchembl_value | decimal | 1:1, Optional | Standardized -log10 activity value | `6.5` |
| natural_product | boolean | 1:1, Optional | Natural product flag | `true`, `false` |
| np_likeness_score | decimal | 1:1, Optional | Natural product likeness score | `0.75` |
| qed_weighted | decimal | 1:1, Optional | QED drug-likeness score | `0.68` |
| confidence_score | integer | 1:1, Optional | Target assignment confidence (0-9) | `9` |

### DailyMed

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| set_id | uuid | 1:1, Required | Unique SPL document identifier | `a1b2c3d4-e5f6-7890-abcd-ef1234567890` |
| ndc | string | 1:1, Optional | National Drug Code (10-11 digits) | `0591-5050-01` |
| spl_sections | object[] | 1:N, Optional | Label sections with LOINC codes | `[{"loinc_code": "34084-4", "content": "..."}]` |
| unii | string | 1:1, Optional | FDA Unique Ingredient Identifier | `R16CO5Y76E` |
| application_number | string | 1:1, Optional | NDA/ANDA/BLA number | `NDA021457` |

### DrugBank

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| drugbank_id | string | 1:1, Required | DrugBank identifier (DB00001 format) | `DB00945` |
| cas_number | string | 1:1, Optional | CAS registry number | `50-78-2` |
| mechanism_of_action | text | 1:1, Optional | Detailed MOA description | `Inhibits cyclooxygenase...` |
| drug_interactions | object[] | 1:N, Optional | Drug-drug interactions with severity | `[{"drug": "Warfarin", "description": "...", "severity": "Major"}]` |
| half_life | string | 1:1, Optional | Elimination half-life | `15-20 minutes` |
| protein_binding | string | 1:1, Optional | Protein binding percentage | `80-90%` |
| pathways | object[] | 1:N, Optional | Associated KEGG/SMPDB pathways | `[{"pathway_id": "SMP0000083", "name": "Aspirin Action Pathway", "source": "SMPDB"}]` |

### FDA Orange Book

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| appl_no | string | 1:1, Required | NDA/ANDA/BLA number | `NDA018281` |
| te_code | string | 1:1, Optional | Therapeutic equivalence code | `AA`, `AB`, `BC`, `BD` |
| patent_no | string | 1:1, Optional | Patent number | `US1234567` |
| patent_expire_date | date | 1:1, Optional | Patent expiration date | `2028-12-31` |
| exclusivity_code | string | 1:1, Optional | Market exclusivity code | `NCE`, `NP`, `ODE` |
| reference_drug | boolean | 1:1, Optional | Is Reference Listed Drug | `true`, `false` |

### RxNorm

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| rxcui | string | 1:1, Required | RxNorm Concept Unique Identifier | `2670` |
| tty | string | 1:1, Optional | Term type | `IN` (Ingredient), `SCD` (Semantic Clinical Drug), `SBD` (Semantic Branded Drug) |
| sab | string | 1:1, Optional | Source abbreviation | `RXNORM`, `MMSL` |
| rela | string | 1:1, Optional | Relationship attribute | `has_ingredient`, `tradename_of` |
| dose_form | string | 1:1, Optional | Dosage form | `Tablet`, `Capsule`, `Solution` |
| strength | string | 1:1, Optional | Drug strength | `500 MG`, `10 MG/ML` |

---

## Metadata Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| _source.database | string | 1:1, Required | Source database name | `ChEMBL`, `DrugBank`, `RxNorm` |
| _source.version | string | 1:1, Optional | Database version | `34`, `5.1.10` |
| _source.access_date | date | 1:1, Optional | Date data was accessed | `2026-01-24` |
| _source.original_id | string | 1:1, Optional | Original identifier from source | - |

---

## Field Mappings by Source

### ChEMBL to Unified Schema

| ChEMBL Field | Unified Field |
|--------------|---------------|
| chembl_id | chembl_id |
| pref_name | name |
| canonical_smiles | canonical_smiles |
| standard_inchi_key | inchi_key |
| full_mwt | molecular_weight |
| molformula | molecular_formula |
| molecule_type | drug_type |
| max_phase | max_phase |
| pchembl_value | pchembl_value |
| natural_product | natural_product |
| np_likeness_score | np_likeness_score |
| qed_weighted | qed_weighted |
| confidence_score | confidence_score |
| molregno | molregno |

### DailyMed to Unified Schema

| DailyMed Field | Unified Field |
|----------------|---------------|
| set_id | set_id |
| product_name | name |
| ndc | ndc |
| spl_section | spl_sections |
| unii | unii |
| application_number | application_number |

### DrugBank to Unified Schema

| DrugBank Field | Unified Field |
|----------------|---------------|
| drugbank_id | drugbank_id |
| name | name |
| smiles | canonical_smiles |
| inchikey | inchi_key |
| average_mass | molecular_weight |
| chemical_formula | molecular_formula |
| drug_type | drug_type |
| status | approval_status |
| cas_number | cas_number |
| mechanism_of_action | mechanism_of_action |
| drug_interactions | drug_interactions |
| half_life | half_life |
| protein_binding | protein_binding |
| targets | targets |
| pathways | pathways |

### Orange Book to Unified Schema

| Orange Book Field | Unified Field |
|-------------------|---------------|
| trade_name | name |
| appl_no | appl_no |
| te_code | te_code |
| patent_no | patent_no |
| patent_expire_date | patent_expire_date |
| exclusivity_code | exclusivity_code |
| rld | reference_drug |

### RxNorm to Unified Schema

| RxNorm Field | Unified Field |
|--------------|---------------|
| rxcui | rxcui |
| name | name |
| tty | tty |
| sab | sab |
| rela | rela |
| dose_form | dose_form |
| strength | strength |

---

## Data Quality Notes

- **Required Field:** Only `name` is strictly required across all sources
- **Structure Data:** Chemical structure fields available primarily from ChEMBL and DrugBank
- **Regulatory Data:** Approval status and patent information from Orange Book and DrugBank
- **Clinical Data:** Drug interactions and mechanism of action primarily from DrugBank
- **Terminology:** RxNorm provides standardized drug naming and classification
- **Label Information:** DailyMed provides comprehensive drug labeling with SPL sections

### Therapeutic Equivalence Codes (Orange Book)

| Code | Description |
|------|-------------|
| AA | No bioequivalence problems |
| AB | Bioequivalent to reference |
| AN | Aerosol, not bioequivalent |
| AO | Injectable oil, not bioequivalent |
| AP | Injectable aqueous, not bioequivalent |
| AT | Topical, not bioequivalent |
| BC | Extended-release, different formulation |
| BD | Different active ingredient |
| BE | Delayed-release, enteric-coated |
| BN | Nebulizer, not bioequivalent |
| BP | Potential bioequivalence problems |
| BR | Controlled release, not bioequivalent |
| BS | Standard, not bioequivalent |
| BT | Topical, not bioequivalent |
| BX | Insufficient data |

### RxNorm Term Types (TTY)

| TTY | Description |
|-----|-------------|
| IN | Ingredient |
| MIN | Multiple Ingredients |
| PIN | Precise Ingredient |
| BN | Brand Name |
| SCD | Semantic Clinical Drug |
| SBD | Semantic Branded Drug |
| SCDC | Semantic Clinical Drug Component |
| SBDC | Semantic Branded Drug Component |
| GPCK | Generic Pack |
| BPCK | Branded Pack |
