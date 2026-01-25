# NPACT - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | npact |
| **Name** | NPACT - Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target |
| **Total Fields** | 22 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| npact_id | String | Yes | NPACT compound identifier | NAP_CPD_001234 |
| compound_name | String | Yes | Chemical compound name | Curcumin |
| plant_id | String | No | Source plant identifier | NAP_ORG_045678 |
| scientific_name | String | No | Latin botanical name | Curcuma longa |
| family | String | No | Botanical family | Zingiberaceae |
| plant_part | Array | No | Plant parts used | ["rhizome"] |

---

## Activity Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| activity_id | String | Activity record ID | NAP_ACT_00567 |
| activity_type | Enum | Activity measurement type | IC50, EC50, GI50 |
| activity_value | Number | Quantitative value | 2.5 |
| activity_unit | String | Unit of measurement | uM, nM |
| cell_line_name | String | Cancer cell line tested | MCF-7 |
| cancer_type | String | Cancer type | Breast |
| assay_type | String | Assay method | MTT assay |
| pubmed_id | Integer | Literature reference | 23456789 |

---

## Chemical Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| pubchem_cid | Integer | PubChem Compound ID |
| cas_number | String | CAS Registry Number |
| smiles | String | SMILES notation |
| molecular_formula | String | Chemical formula |
| molecular_weight | Number | Molecular weight (Da) |
| bioavailability_score | Number | Predicted oral bioavailability |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Compound ID | NAP_CPD_###### | NPACT compound ID | NAP_CPD_001234 |
| Organism ID | NAP_ORG_###### | NPACT organism ID | NAP_ORG_045678 |
| Activity ID | NAP_ACT_###### | NPACT activity ID | NAP_ACT_00567 |
| PubMed ID | Numeric | Literature reference | 23456789 |

---

## Enumerations

### Activity Types

| Value | Description |
|-------|-------------|
| IC50 | Half-maximal inhibitory concentration |
| EC50 | Half-maximal effective concentration |
| GI50 | Growth inhibition 50% |
| LC50 | Lethal concentration 50% |
| CC50 | Cytotoxic concentration 50% |

### Cancer Types

| Value | Description |
|-------|-------------|
| Breast | Breast cancer |
| Lung | Lung cancer |
| Colon | Colorectal cancer |
| Leukemia | Blood cancer |
| Liver | Hepatocellular carcinoma |
| Cervical | Cervical cancer |
| Prostate | Prostate cancer |

---

## Entity Relationships

### Plant to Compound
- **Cardinality:** 1:N
- **Description:** Plants contain multiple anti-cancer compounds
- **Key Fields:** plant_id, npact_id

### Compound to Activity
- **Cardinality:** 1:N
- **Description:** Compounds have multiple activity measurements
- **Key Fields:** npact_id, activity_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NPACT | Naturally occurring Plant-based Anti-cancer Compound-Activity-Target | Full name |
| IC50 | Inhibitory Concentration 50% | Activity metric |
| GI50 | Growth Inhibition 50% | Activity metric |
| MTT | 3-(4,5-dimethylthiazol-2-yl)-2,5-diphenyl tetrazolium bromide | Assay type |
| SRB | Sulforhodamine B | Assay type |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
