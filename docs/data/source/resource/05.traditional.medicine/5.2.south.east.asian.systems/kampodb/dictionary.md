# KampoDB - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | kampodb |
| **Name** | KampoDB - Kampo Medicine Database |
| **Total Fields** | 30 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| formula_id | String | Yes | Kampo formula code | KT |
| formula_name | String | No | Romanized formula name | Kakkonto |
| formula_name_jp | String | No | Japanese kanji name | 葛根湯 |
| crude_drug_id | Integer | No | Crude drug ID (0-179) | 40 |
| crude_drug_name | String | No | English crude drug name | Cinnamon Bark |
| compound_id | Integer | No | PubChem CID | 2244 |
| compound_name | String | No | Chemical name | Salicylic acid |
| protein_id | Integer | No | NCBI Gene ID | 5743 |
| gene_symbol | String | No | Gene symbol | PTGS2 |

---

## Docking Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| docking.domain_id | Integer | Protein domain ID | 1 |
| docking.affinity_kcal_mol | Number | Binding energy (kcal/mol) | -7.2 |
| docking.ligand_atoms | Integer | Ligand atom count | 12 |
| docking.protein_atoms | Integer | Protein atom count | 580 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Formula ID | String code | Formula abbreviation | KT, YKS |
| Crude Drug ID | 0-179 | Numeric crude drug ID | 40 |
| Compound ID | Numeric | PubChem CID | 969516 |
| Protein ID | Numeric | NCBI Gene ID | 2475 |
| Pathway ID | hsa##### | KEGG pathway | hsa04668 |

---

## 4-Layer Hierarchy

| Layer | Entity | Count | Description |
|-------|--------|-------|-------------|
| 1 | Formula | 298 | Kampo prescriptions |
| 2 | Crude Drug | 180 | Raw medicinal materials |
| 3 | Compound | 3,002 | Natural chemicals |
| 4 | Protein | 62,906 | Human target proteins |

---

## Enumerations

### Target Prediction Types

| Value | Description |
|-------|-------------|
| known | Experimentally validated |
| predicted | Computationally predicted |
| docking | Based on docking simulation |

---

## Entity Relationships

### Formula to Crude Drug
- **Cardinality:** N:M
- **Description:** Formulas contain multiple crude drugs
- **Key Fields:** formula_id, crude_drug_id

### Crude Drug to Compound
- **Cardinality:** 1:N
- **Description:** Crude drugs contain compounds
- **Key Fields:** crude_drug_id, compound_id

### Compound to Protein
- **Cardinality:** N:M
- **Description:** Compounds interact with proteins (via docking)
- **Key Fields:** compound_id, protein_id

---

## API Endpoints

| Endpoint | Description |
|----------|-------------|
| /api/formula/ | List all formulas |
| /api/formula/{id}/info | Formula details |
| /api/formula/{id}/compound | Formula compounds |
| /api/docking/compound/{id}/protein/{id} | Docking score |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| KampoDB | Kampo Database | Database name |
| Kampo | Japanese traditional medicine | Medical system |
| CID | Compound Identifier | PubChem ID type |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
