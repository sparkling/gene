# TCMSID - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | tcmsid |
| **Name** | TCM Structural Identification Database |
| **Total Fields** | 25 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| compound_id | String | Yes | TCMSID compound identifier | TCMSID00001 |
| compound_name | String | Yes | Compound name | Ginsenoside Rg1 |
| herb_id | String | No | Source herb identifier | TCMSID_H00001 |
| herb_name | String | No | Herb name | Ginseng |
| herb_name_cn | String | No | Chinese herb name | 人参 |

---

## Structural Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| smiles | String | SMILES notation | CC(CCC=C(C)C)... |
| inchi | String | InChI identifier | InChI=1S/C42H72O14/... |
| inchi_key | String | InChIKey hash | ABCDEFGHIJ-KLMNOPQRST-U |
| molecular_formula | String | Chemical formula | C42H72O14 |
| molecular_weight | Number | Molecular weight | 801.01 |
| exact_mass | Number | Monoisotopic mass | 800.4922 |
| 2d_structure | String | SDF/MOL file path | /structures/TCMSID00001.sdf |
| 3d_structure | String | 3D structure file path | /structures/TCMSID00001.mol2 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Compound ID | TCMSID##### | Compound identifier | TCMSID00001 |
| Herb ID | TCMSID_H##### | Herb identifier | TCMSID_H00001 |
| PubChem CID | Numeric | Cross-reference | 441923 |
| ChEMBL ID | CHEMBL#### | Cross-reference | CHEMBL12345 |

---

## Structural Properties

| Property | Description | Range |
|----------|-------------|-------|
| HBD | Hydrogen bond donors | 0-20 |
| HBA | Hydrogen bond acceptors | 0-30 |
| LogP | Partition coefficient | -5 to 10 |
| TPSA | Topological polar surface area | 0-300 |
| RotBonds | Rotatable bonds | 0-20 |
| RingCount | Number of rings | 0-10 |

---

## Entity Relationships

### Compound to Herb
- **Cardinality:** N:M
- **Description:** Compounds found in multiple herbs
- **Key Fields:** compound_id, herb_id

### Compound to Structure
- **Cardinality:** 1:1
- **Description:** Each compound has unique structure
- **Key Fields:** compound_id, inchi_key

---

## File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| SDF | .sdf | 2D structure |
| MOL2 | .mol2 | 3D structure |
| PDB | .pdb | Protein structure |
| SMILES | .smi | Linear notation |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| TCMSID | TCM Structural Identification Database | Full name |
| SMILES | Simplified Molecular Input Line Entry System | Notation |
| InChI | International Chemical Identifier | Identifier |
| SDF | Structure Data File | File format |
| TPSA | Topological Polar Surface Area | Property |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
