# COCONUT - Data Dictionary

## Overview

This data dictionary documents the schema for COCONUT (COlleCtion of Open Natural prodUcTs) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | coconut |
| **Name** | COCONUT |
| **Parent** | 2.1.natural.products |
| **Total Fields** | 40+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| coconut_id | string | 1:1 | Yes | Unique COCONUT identifier | CNP0123456 |
| name | string | 1:1 | No | Common or trivial name | Artemisinin |
| iupac_name | string | 1:1 | No | IUPAC systematic name | Full systematic name |
| molecular_formula | string | 1:1 | No | Chemical formula | C21H30O2 |
| molecular_weight | decimal | 1:1 | No | Molecular weight in Daltons | 314.46 |

### Structure Representations

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| canonical_smiles | string | 1:1 | No | Standardized SMILES | CC(C)CCCC(C)C |
| isomeric_smiles | string | 1:1 | No | SMILES with stereochemistry | C[C@@H](O)c1ccccc1 |
| inchi | string | 1:1 | No | International Chemical Identifier | InChI=1S/C15H24O5/... |
| inchi_key | string | 1:1 | No | 27-character InChI hash | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| sugar_free_smiles | string | 1:1 | No | SMILES without glycosidic moieties | Aglycone structure |

### Physicochemical Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| total_atom_number | integer | 1:1 | No | Total atoms in molecule | 45 |
| heavy_atom_number | integer | 1:1 | No | Non-hydrogen atoms | 25 |
| rotatable_bond_number | integer | 1:1 | No | Conformational flexibility | 5 |
| hbond_acceptor | integer | 1:1 | No | H-bond accepting groups | 5 |
| hbond_donor | integer | 1:1 | No | H-bond donating groups | 2 |
| number_of_rings | integer | 1:1 | No | Total ring count | 3 |
| aromatic_rings_count | integer | 1:1 | No | Aromatic ring count | 1 |
| qed_drug_likeliness | decimal | 1:1 | No | Drug-likeness score (0-1) | 0.85 |
| lipinski_rule_of_5 | integer | 1:1 | No | Lipinski violations (0-4) | 0 |
| topological_polar_surface_area | decimal | 1:1 | No | TPSA in Angstrom squared | 40.46 |
| alogp | decimal | 1:1 | No | Calculated partition coefficient | 2.5 |
| xlogp | decimal | 1:1 | No | Alternative logP calculation | 2.8 |
| murko_framework | string | 1:1 | No | Core scaffold structure | Benzene ring |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| COCONUT ID | CNP + 7 digits | CNP0123456 | Primary identifier |
| InChI Key | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | Structure hash |
| CAS Number | varies | 50-78-2 | CAS Registry (when available) |

---

## Enumerations

### Data Sources

| Source | Description |
|--------|-------------|
| ZINC | ZINC Database |
| ChEMBL | Bioactivity Database |
| PubChem | NIH Chemical Database |
| LOTUS | Natural Products Occurrence |
| NPASS | Natural Products Activity |
| NAPRALERT | Natural Products Database |
| TCMDB | Traditional Chinese Medicine |

### Biological Activity Types

| Activity Type | Description |
|---------------|-------------|
| IC50 | Half maximal inhibitory concentration |
| EC50 | Half maximal effective concentration |
| Ki | Inhibition constant |
| MIC | Minimum inhibitory concentration |
| LD50 | Median lethal dose |
| ED50 | Median effective dose |

### Taxonomic Kingdoms

| Kingdom | Description |
|---------|-------------|
| plant | Plant-derived natural products |
| bacteria | Bacterial metabolites |
| fungi | Fungal metabolites |
| animal | Animal-derived compounds |

---

## Entity Relationships

### Compound to Properties
- **Cardinality:** 1:1
- **Description:** Each compound has one set of physicochemical properties
- **Key Fields:** compound_id, properties.compound_id

### Compound to Sources
- **Cardinality:** N:M
- **Description:** Compounds can be found in multiple source databases
- **Key Fields:** compound_id, source_id, external_id

### Compound to Organisms
- **Cardinality:** N:M
- **Description:** Compounds can be isolated from multiple organisms
- **Key Fields:** compound_id, organism_id

### Compound to Biological Activities
- **Cardinality:** 1:N
- **Description:** Compounds may have multiple documented activities
- **Key Fields:** compound_id, activity_id

### Compound to Citations
- **Cardinality:** 1:N
- **Description:** Compounds may be referenced in multiple publications
- **Key Fields:** compound_id, doi, pmid

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| COCONUT | COlleCtion of Open Natural prodUcTs | Database name |
| NP | Natural Product | Compound type |
| SMILES | Simplified Molecular Input Line Entry System | Structure notation |
| InChI | International Chemical Identifier | IUPAC identifier |
| QED | Quantitative Estimate of Drug-likeness | Drug-likeness score |
| TPSA | Topological Polar Surface Area | Absorption predictor |
| HBD | Hydrogen Bond Donor | Lipinski parameter |
| HBA | Hydrogen Bond Acceptor | Lipinski parameter |
| MW | Molecular Weight | Daltons |
| SDF | Structure-Data File | Chemical file format |
| NCBI | National Center for Biotechnology Information | Taxonomy source |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| ChEMBL | ChEMBL ID | Bioactivity data |
| PubChem | PubChem CID | Chemical information |
| ZINC | ZINC ID | Commercial availability |
| ChEBI | ChEBI ID | Chemical ontology |
| HMDB | HMDB ID | Human metabolome |
| DrugBank | DrugBank ID | Drug information |
| FooDB | FooDB ID | Food composition |

---

## Data Quality Notes

1. **Cardinality:** Fields marked 1:1 are single-valued; 1:N are multi-valued arrays
2. **N:M Relationships:** Compound-source and compound-organism are many-to-many
3. **Cross-References:** InChI Key enables multi-database integration
4. **Update Frequency:** Quarterly releases with new compounds
5. **Validation:** Structures validated against RDKit/CDK
