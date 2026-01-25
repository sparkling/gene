# BindingDB - Data Dictionary

## Overview

This data dictionary documents the schema for BindingDB binding affinity database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | bindingdb |
| **Name** | BindingDB |
| **Parent** | 2.7.compound.target.interactions |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Ligand Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| monomer_id | integer | 1:1 | Yes | BindingDB ligand ID | 12345 |
| name | string | 1:1 | No | Compound name | Example Inhibitor |
| smiles | string | 1:1 | Yes | SMILES structure | CC(=O)Oc1ccccc1C(=O)O |
| inchi | string | 1:1 | No | InChI identifier | InChI=1S/C9H8O4/... |
| inchi_key | string | 1:1 | Yes | InChI Key | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| molecular_weight | decimal | 1:1 | No | MW in Daltons | 180.16 |
| molecular_formula | string | 1:1 | No | Chemical formula | C9H8O4 |
| alogp | decimal | 1:1 | No | Calculated LogP | 1.31 |
| psa | decimal | 1:1 | No | Polar surface area | 63.6 |
| hbd | integer | 1:1 | No | H-bond donors | 1 |
| hba | integer | 1:1 | No | H-bond acceptors | 4 |
| rotatable_bonds | integer | 1:1 | No | Rotatable bond count | 3 |

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| target_id | integer | 1:1 | Yes | BindingDB target ID | 5001 |
| name | string | 1:1 | Yes | Target name | Cyclooxygenase-2 |
| organism | string | 1:1 | Yes | Species | Homo sapiens |
| uniprot_id | string | 1:1 | No | UniProt accession | P35354 |
| gene_symbol | string | 1:1 | No | Gene name | PTGS2 |
| pdb_ids | array | 1:N | No | PDB structure IDs | [1CX2, 3LN1] |
| ec_number | string | 1:1 | No | EC classification | 1.14.99.1 |
| target_type | string | 1:1 | No | Target class | Enzyme |

### Binding Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| bindingdb_id | integer | 1:1 | Yes | Primary identifier | 50000001 |
| monomer_id | integer | 1:1 | Yes | Foreign key to ligands | 12345 |
| target_id | integer | 1:1 | Yes | Foreign key to targets | 5001 |
| affinity_type | string | 1:1 | Yes | Measurement type | IC50, Ki, Kd, EC50 |
| affinity_value | decimal | 1:1 | Yes | Numeric value | 0.6 |
| affinity_unit | string | 1:1 | Yes | Unit of measurement | nM, uM, mM |
| affinity_relation | string | 1:1 | No | Value relation | =, <, >, ~ |
| ph | decimal | 1:1 | No | Assay pH | 7.4 |
| temperature | decimal | 1:1 | No | Temperature (C) | 25 |
| assay_description | text | 1:1 | No | Assay method | - |
| reference_id | integer | 1:1 | No | Literature source | 101 |
| source_type | string | 1:1 | No | Source type | Publication, Patent |

### Reference Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| reference_id | integer | 1:1 | Yes | Primary identifier | 101 |
| pmid | integer | 1:1 | No | PubMed ID | 12345678 |
| doi | string | 1:1 | No | DOI | 10.1021/jm001234 |
| patent_number | string | 1:1 | No | Patent number | US7654321 |
| title | string | 1:1 | No | Article/patent title | - |
| authors | string | 1:1 | No | Author list | - |
| year | integer | 1:1 | No | Publication year | 2020 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Monomer ID | Integer | 12345 | Ligand identifier |
| Target ID | Integer | 5001 | Target identifier |
| BindingDB ID | Integer | 50000001 | Binding record ID |
| UniProt ID | Accession | P35354 | Protein identifier |
| PDB ID | 4 characters | 1CX2 | Structure identifier |
| InChI Key | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | Structure hash |
| PubMed ID | Integer | 12345678 | Literature reference |

---

## Enumerations

### Affinity Types

| Type | Description | Typical Range | Use Case |
|------|-------------|---------------|----------|
| IC50 | Half-maximal inhibitory concentration | nM - uM | Enzyme inhibition |
| Ki | Inhibition constant | nM - uM | True binding affinity |
| Kd | Dissociation constant | nM - uM | Binding equilibrium |
| EC50 | Half-maximal effective concentration | nM - uM | Functional response |
| Kon | Association rate | M-1s-1 | Binding kinetics |
| Koff | Dissociation rate | s-1 | Binding kinetics |

### Affinity Units

| Unit | Description | Conversion |
|------|-------------|------------|
| nM | Nanomolar | 10^-9 M |
| uM | Micromolar | 10^-6 M |
| mM | Millimolar | 10^-3 M |
| pM | Picomolar | 10^-12 M |

### Affinity Relations

| Relation | Description |
|----------|-------------|
| = | Exact value |
| < | Less than |
| > | Greater than |
| ~ | Approximate |
| <= | Less than or equal |
| >= | Greater than or equal |

### Target Types

| Type | Description |
|------|-------------|
| Enzyme | Catalytic protein |
| GPCR | G protein-coupled receptor |
| Kinase | Protein kinase |
| Ion channel | Ion channel protein |
| Transporter | Membrane transporter |
| Nuclear receptor | NHR |
| Other | Other target types |

### Source Types

| Type | Description |
|------|-------------|
| Publication | Peer-reviewed article |
| Patent | Patent document |
| Thesis | Academic thesis |
| Book | Book chapter |

---

## Entity Relationships

### Ligand to Binding Data
- **Cardinality:** 1:N
- **Description:** One ligand tested against multiple targets
- **Key Fields:** monomer_id, bindingdb_id

### Target to Binding Data
- **Cardinality:** 1:N
- **Description:** One target has multiple ligand measurements
- **Key Fields:** target_id, bindingdb_id

### Binding Data to Reference
- **Cardinality:** N:1
- **Description:** Multiple measurements from one publication
- **Key Fields:** bindingdb_id, reference_id

### Target to PDB Structures
- **Cardinality:** 1:N
- **Description:** One target may have multiple crystal structures
- **Key Fields:** target_id, pdb_ids

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IC50 | Half-maximal Inhibitory Concentration | Potency measure |
| Ki | Inhibition Constant | Binding affinity |
| Kd | Dissociation Constant | Equilibrium binding |
| EC50 | Half-maximal Effective Concentration | Functional potency |
| pIC50 | Negative log10 of IC50 | -log10(IC50 in M) |
| pKi | Negative log10 of Ki | -log10(Ki in M) |
| SAR | Structure-Activity Relationship | Drug design analysis |
| TSV | Tab-Separated Values | Data format |
| SDF | Structure-Data File | Chemical format |
| GPCR | G Protein-Coupled Receptor | Target class |
| PDB | Protein Data Bank | Structure database |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UniProt | Accession | Target proteins |
| PDB | PDB ID | 3D structures |
| PubChem | CID | Chemical data |
| ChEMBL | ChEMBL ID | Bioactivity |
| PubMed | PMID | Literature |
| DrugBank | DrugBank ID | Drug information |

---

## Data Quality Notes

1. **Quantitative Data:** 2.9M+ binding measurements
2. **Ligand Coverage:** 1.3M+ distinct compounds
3. **Target Coverage:** ~9,400 protein targets
4. **Literature Support:** 40,000+ publications, 35,000+ patents
5. **Weekly Updates:** Regular data additions
6. **Experimental Conditions:** pH and temperature recorded
7. **Multiple Affinity Types:** IC50, Ki, Kd, EC50, kinetics
8. **REST API:** Programmatic access available

