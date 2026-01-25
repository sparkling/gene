# TCMBank - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | tcmbank |
| **Name** | TCMBank Database |
| **Total Fields** | 28 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| herb_id | String | Yes | TCMBank herb identifier | TB_H000001 |
| herb_name_en | String | Yes | English herb name | Ginseng |
| herb_name_cn | String | No | Chinese name | 人参 |
| herb_name_latin | String | No | Latin name | Panax ginseng |
| ingredient_id | String | No | Ingredient identifier | TB_I000001 |
| ingredient_name | String | No | Compound name | Ginsenoside Rg1 |
| target_id | String | No | Target identifier | TB_T000001 |
| target_name | String | No | Target protein name | TP53 |

---

## Chemical Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| pubchem_cid | Integer | PubChem Compound ID | 441923 |
| cas_number | String | CAS Registry Number | 22427-39-0 |
| smiles | String | SMILES notation | CC(CCC=C... |
| molecular_formula | String | Chemical formula | C42H72O14 |
| molecular_weight | Number | Molecular weight (Da) | 801.01 |
| drug_likeness | Object | Drug-likeness properties | {lipinski: true} |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Herb ID | TB_H###### | Herb identifier | TB_H000001 |
| Ingredient ID | TB_I###### | Ingredient identifier | TB_I000001 |
| Target ID | TB_T###### | Target identifier | TB_T000001 |
| Disease ID | TB_D###### | Disease identifier | TB_D000001 |

---

## Entity Relationships

### Herb to Ingredient
- **Cardinality:** 1:N
- **Description:** Herbs contain multiple ingredients
- **Key Fields:** herb_id, ingredient_id

### Ingredient to Target
- **Cardinality:** N:M
- **Description:** Ingredients interact with targets
- **Key Fields:** ingredient_id, target_id

### Target to Disease
- **Cardinality:** N:M
- **Description:** Targets linked to diseases
- **Key Fields:** target_id, disease_id

---

## Cross-References

| Database | Field | Description |
|----------|-------|-------------|
| PubChem | pubchem_cid | Compound structures |
| CAS | cas_number | Chemical registry |
| UniProt | uniprot_id | Protein data |
| OMIM | omim_id | Genetic diseases |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| TCMBank | Traditional Chinese Medicine Bank | Database name |
| TCM | Traditional Chinese Medicine | Medical system |
| TB | TCMBank | ID prefix |
| CAS | Chemical Abstracts Service | Chemical registry |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
