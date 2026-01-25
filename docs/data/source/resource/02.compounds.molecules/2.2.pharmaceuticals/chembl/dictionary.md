# ChEMBL - Data Dictionary

## Overview

This data dictionary documents the schema for ChEMBL bioactivity database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | chembl |
| **Name** | ChEMBL |
| **Parent** | 2.2.pharmaceuticals |
| **Total Fields** | 100+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Molecule Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| molregno | integer | 1:1 | Yes | Internal molecule registration number | 1234567 |
| chembl_id | string | 1:1 | Yes | Public ChEMBL identifier | CHEMBL25 |
| pref_name | string | 1:1 | No | Preferred compound name | Aspirin |
| max_phase | numeric | 1:1 | No | Highest clinical phase (0-4) | 4 |
| therapeutic_flag | smallint | 1:1 | No | Is therapeutic (0/1) | 1 |
| natural_product | smallint | 1:1 | No | Natural product flag (0/1) | 0 |
| first_approval | smallint | 1:1 | No | Year of first approval | 1950 |

### Structure Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| canonical_smiles | string | 1:1 | No | Canonical SMILES notation | CC(=O)Oc1ccccc1C(=O)O |
| standard_inchi | string | 1:1 | No | Standard InChI string | InChI=1S/C9H8O4/... |
| standard_inchi_key | string | 1:1 | No | InChI Key (27 char) | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| molfile | text | 1:1 | No | MDL MOL file format | MOL block |

### Compound Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| mw_freebase | numeric | 1:1 | No | Molecular weight (freebase) | 180.16 |
| alogp | numeric | 1:1 | No | Calculated LogP (Wildman-Crippen) | 1.31 |
| hba | smallint | 1:1 | No | Hydrogen bond acceptors | 3 |
| hbd | smallint | 1:1 | No | Hydrogen bond donors | 1 |
| psa | numeric | 1:1 | No | Polar surface area | 63.60 |
| rtb | smallint | 1:1 | No | Rotatable bonds | 2 |
| num_ro5_violations | smallint | 1:1 | No | Lipinski violations | 0 |
| qed_weighted | numeric | 1:1 | No | QED drug-likeness score | 0.55 |
| np_likeness_score | numeric | 1:1 | No | Natural product likeness | -0.86 |
| aromatic_rings | smallint | 1:1 | No | Number of aromatic rings | 1 |
| heavy_atoms | smallint | 1:1 | No | Heavy atom count | 13 |

### Activity Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| activity_id | integer | 1:1 | Yes | Primary key | 12345 |
| standard_type | string | 1:1 | No | Activity type | IC50, Ki, EC50 |
| standard_relation | string | 1:1 | No | Relation operator | =, <, > |
| standard_value | numeric | 1:1 | No | Standardized value | 300.0 |
| standard_units | string | 1:1 | No | Units | nM |
| pchembl_value | numeric | 1:1 | No | -log10 of molar activity | 6.52 |

### Target Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| tid | integer | 1:1 | Yes | Target identifier | 100126 |
| target_type | string | 1:1 | No | Target classification | SINGLE PROTEIN |
| pref_name | string | 1:1 | Yes | Preferred name | Cyclooxygenase-1 |
| organism | string | 1:1 | No | Organism name | Homo sapiens |
| tax_id | integer | 1:1 | No | NCBI taxonomy ID | 9606 |
| accession | string | 1:1 | No | UniProt accession | P23219 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| ChEMBL ID | CHEMBL + digits | CHEMBL25 | Molecule/target/assay ID |
| molregno | Integer | 1234567 | Internal molecule ID |
| InChI Key | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | Structure hash |
| UniProt | Accession | P23219 | Protein target |
| PubMed ID | Integer | 12345678 | Literature reference |

---

## Enumerations

### Max Phase Values

| Value | Description |
|-------|-------------|
| 0 | Preclinical |
| 0.5 | Early clinical |
| 1 | Phase I |
| 2 | Phase II |
| 3 | Phase III |
| 4 | Approved |

### Assay Types

| Code | Description |
|------|-------------|
| B | Binding |
| F | Functional |
| A | ADMET |
| T | Toxicity |
| P | Physicochemical |
| U | Unassigned |

### Confidence Scores (Target Assignment)

| Score | Description |
|-------|-------------|
| 9 | Direct single protein target |
| 8 | Homologous single protein |
| 7 | Direct protein complex |
| 6 | Homologous protein complex |
| 5 | Direct selectivity assay |
| 4 | Indirect protein complex |
| 3 | Cell-based |
| 2 | Tissue/whole organism |
| 1 | Uncurated |
| 0 | Target not assigned |

### Target Types

| Type | Description |
|------|-------------|
| SINGLE PROTEIN | Single protein target |
| PROTEIN COMPLEX | Multi-protein complex |
| PROTEIN FAMILY | Related protein group |
| CELL-LINE | Cellular target |
| ORGANISM | Whole organism |

### Drug Mechanism Actions

| Action | Description |
|--------|-------------|
| INHIBITOR | Blocks target activity |
| AGONIST | Activates receptor |
| ANTAGONIST | Blocks receptor |
| MODULATOR | Modifies activity |
| SUBSTRATE | Metabolized by target |

---

## Entity Relationships

### Molecule to Activities
- **Cardinality:** 1:N
- **Description:** Each molecule can have multiple bioactivity measurements
- **Key Fields:** molregno, activity_id

### Activity to Assay
- **Cardinality:** N:1
- **Description:** Multiple activities measured in one assay
- **Key Fields:** activity_id, assay_id

### Assay to Target
- **Cardinality:** N:1
- **Description:** Multiple assays for one target
- **Key Fields:** assay_id, tid

### Target to Components
- **Cardinality:** 1:N
- **Description:** Targets composed of protein components
- **Key Fields:** tid, component_id

### Molecule to Drug Mechanism
- **Cardinality:** 1:N
- **Description:** Approved drugs have documented mechanisms
- **Key Fields:** molregno, mec_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ChEMBL | Chemical database of bioactive molecules | EMBL-EBI database |
| SMILES | Simplified Molecular Input Line Entry System | Structure notation |
| InChI | International Chemical Identifier | IUPAC standard |
| IC50 | Half-maximal Inhibitory Concentration | Activity measure |
| Ki | Inhibition Constant | Binding affinity |
| EC50 | Half-maximal Effective Concentration | Potency measure |
| pChEMBL | -log10 molar activity | Standardized activity |
| HBA | Hydrogen Bond Acceptor | Lipinski parameter |
| HBD | Hydrogen Bond Donor | Lipinski parameter |
| PSA | Polar Surface Area | Absorption predictor |
| RTB | Rotatable Bonds | Flexibility measure |
| MW | Molecular Weight | Daltons |
| QED | Quantitative Estimate of Drug-likeness | 0-1 score |
| MOA | Mechanism of Action | Drug-target interaction |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity | PK/PD |
| BAO | BioAssay Ontology | Assay classification |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UniProt | Accession | Protein targets |
| PubChem | CID | Chemical data |
| DrugBank | DrugBank ID | Drug information |
| PDB | PDB ID | 3D structures |
| InterPro | IPR ID | Protein domains |
| Gene Ontology | GO ID | Functional annotation |
| Reactome | Pathway ID | Metabolic pathways |
| ChEBI | ChEBI ID | Chemical ontology |
| MeSH | MeSH ID | Disease terms |
| EFO | EFO ID | Disease classification |

---

## Data Quality Notes

1. **pChEMBL Values:** Standardized -log10 activity enables cross-assay comparison
2. **Confidence Scores:** Higher scores indicate more reliable target assignment
3. **Data Validity:** Flags indicate suspicious or potentially duplicate entries
4. **Natural Products:** Binary flag plus NP-likeness score for structural similarity
5. **License:** CC BY-SA 3.0 (attribution and share-alike required)
6. **Update Frequency:** Quarterly major releases (ChEMBL 36 = July 2025)
