# GtoPdb - Data Dictionary

## Overview

This data dictionary documents the schema for GtoPdb (Guide to Pharmacology database).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gtopdb |
| **Name** | GtoPdb |
| **Parent** | 2.7.compound.target.interactions |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| target_id | integer | 1:1 | Yes | GtoPdb target ID | 290 |
| name | string | 1:1 | Yes | Target name | Serine/threonine-protein kinase B-raf |
| abbreviation | string | 1:1 | No | Short name | BRAF |
| systematic_name | string | 1:1 | No | IUPHAR systematic name | - |
| family_id | integer | 1:1 | Yes | Target family | 538 |
| type | string | 1:1 | Yes | Target type | kinase |
| hgnc_symbol | string | 1:1 | No | HGNC gene symbol | BRAF |
| hgnc_id | integer | 1:1 | No | HGNC ID | 1097 |
| uniprot_id | string | 1:1 | No | UniProt accession | P15056 |
| species | string | 1:1 | Yes | Organism | Human |

### Ligand Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ligand_id | integer | 1:1 | Yes | GtoPdb ligand ID | 5085 |
| name | string | 1:1 | Yes | Ligand name | vemurafenib |
| type | string | 1:1 | Yes | Ligand type | synthetic |
| approved | boolean | 1:1 | No | Approved drug status | true |
| radioactive | boolean | 1:1 | No | Radioligand | false |
| pubchem_sid | integer | 1:1 | No | PubChem SID | 178100714 |
| pubchem_cid | integer | 1:1 | No | PubChem CID | 42611257 |
| chembl_id | string | 1:1 | No | ChEMBL ID | CHEMBL1667 |
| inn | string | 1:1 | No | International nonproprietary name | vemurafenib |
| smiles | string | 1:1 | No | Chemical structure | CCCS(=O)(=O)Nc1ccc... |
| inchi | string | 1:1 | No | InChI | InChI=1S/... |
| inchi_key | string | 1:1 | No | InChI Key | GPXBXXGIAQBQNI-UHFFFAOYSA-N |

### Interaction Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| interaction_id | integer | 1:1 | Yes | Primary identifier | 10001 |
| target_id | integer | 1:1 | Yes | Foreign key to targets | 290 |
| ligand_id | integer | 1:1 | Yes | Foreign key to ligands | 5085 |
| type | string | 1:1 | Yes | Interaction type | inhibitor |
| action | string | 1:1 | No | Action type | Inhibition |
| selectivity | string | 1:1 | No | Selectivity class | Selective |
| endogenous | boolean | 1:1 | No | Endogenous ligand flag | false |
| affinity_type | string | 1:1 | No | Affinity notation | pKi, pIC50, pEC50 |
| affinity_median | decimal | 1:1 | No | Median pX value | 8.3 |
| affinity_low | decimal | 1:1 | No | Low pX value | 7.8 |
| affinity_high | decimal | 1:1 | No | High pX value | 8.8 |
| species | string | 1:1 | No | Assay species | Human |
| reference_ids | array | 1:N | No | Literature references | [1234, 5678] |

### Target Family Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| family_id | integer | 1:1 | Yes | Family identifier | 538 |
| name | string | 1:1 | Yes | Family name | Kinases |
| parent_id | integer | 1:1 | No | Parent family | - |
| type | string | 1:1 | Yes | Family type | kinase |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Target ID | Integer | 290 | GtoPdb target identifier |
| Ligand ID | Integer | 5085 | GtoPdb ligand identifier |
| Family ID | Integer | 538 | Target family identifier |
| HGNC ID | Integer | 1097 | Gene nomenclature |
| UniProt ID | Accession | P15056 | Protein identifier |
| ChEMBL ID | CHEMBL + digits | CHEMBL1667 | Chemical database |
| PubChem CID | Integer | 42611257 | Compound identifier |
| InChI Key | 27 characters | GPXBXXGIAQBQNI-UHFFFAOYSA-N | Structure hash |

---

## Enumerations

### Target Family Types

| Family | Description | Examples |
|--------|-------------|----------|
| GPCRs | G protein-coupled receptors | Adrenoceptors, opioid receptors |
| Ion channels | Voltage/ligand-gated channels | Nav, Cav, nAChR |
| NHRs | Nuclear hormone receptors | Estrogen receptor, PPAR |
| Kinases | Protein kinases | EGFR, BRAF |
| Enzymes | Other enzymes | COX, PDE |
| Transporters | SLC, ABC transporters | SERT, P-gp |
| Other proteins | Miscellaneous targets | Various |

### Interaction Types

| Type | Description |
|------|-------------|
| agonist | Activates receptor |
| antagonist | Blocks receptor |
| inhibitor | Inhibits enzyme/kinase |
| activator | Activates enzyme |
| allosteric modulator | Allosteric site binding |
| partial agonist | Partial receptor activation |
| inverse agonist | Inverse receptor activity |
| channel blocker | Blocks ion channel |

### Ligand Types

| Type | Description |
|------|-------------|
| synthetic | Synthetic compound |
| metabolite | Endogenous metabolite |
| endogenous | Endogenous ligand |
| natural product | Natural product origin |
| peptide | Peptide ligand |
| antibody | Monoclonal antibody |
| inorganic | Inorganic compound |

### Action Types

| Action | Description |
|--------|-------------|
| Activation | Increases target activity |
| Inhibition | Decreases target activity |
| Modulation | Modifies target function |
| Binding | Binds without functional effect |

### Selectivity Classes

| Class | Description |
|-------|-------------|
| Selective | Selective for target |
| Non-selective | Acts on multiple targets |

### Affinity Notation (pX)

| Type | Description | Conversion |
|------|-------------|------------|
| pKi | -log10(Ki in M) | Ki (nM) = 10^(9-pKi) |
| pIC50 | -log10(IC50 in M) | IC50 (nM) = 10^(9-pIC50) |
| pEC50 | -log10(EC50 in M) | EC50 (nM) = 10^(9-pEC50) |
| pKd | -log10(Kd in M) | Kd (nM) = 10^(9-pKd) |

---

## Entity Relationships

### Target to Family
- **Cardinality:** N:1
- **Description:** Targets belong to families
- **Key Fields:** target_id, family_id

### Target to Interactions
- **Cardinality:** 1:N
- **Description:** Each target has multiple ligand interactions
- **Key Fields:** target_id, interaction_id

### Ligand to Interactions
- **Cardinality:** 1:N
- **Description:** Each ligand has multiple target interactions
- **Key Fields:** ligand_id, interaction_id

### Family Hierarchy
- **Cardinality:** N:1
- **Description:** Families can have parent families
- **Key Fields:** family_id, parent_id

### Interaction to References
- **Cardinality:** 1:N
- **Description:** Interactions linked to literature
- **Key Fields:** interaction_id, reference_ids

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IUPHAR | International Union of Basic and Clinical Pharmacology | Standards body |
| BPS | British Pharmacological Society | Database partner |
| GtoPdb | Guide to Pharmacology Database | Database name |
| pKi | Negative log10 of Ki | Affinity notation |
| pIC50 | Negative log10 of IC50 | Affinity notation |
| pEC50 | Negative log10 of EC50 | Affinity notation |
| INN | International Nonproprietary Name | Drug naming |
| GPCR | G Protein-Coupled Receptor | Target class |
| NHR | Nuclear Hormone Receptor | Target class |
| SLC | Solute Carrier | Transporter family |
| ABC | ATP-Binding Cassette | Transporter family |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UniProt | Accession | Target proteins |
| PubChem | CID/SID | Chemical data |
| ChEMBL | ChEMBL ID | Bioactivity |
| HGNC | HGNC ID | Gene nomenclature |
| PDB | PDB ID | 3D structures |
| Ensembl | Gene ID | Gene annotation |

---

## Data Quality Notes

1. **Expert Curation:** IUPHAR/BPS curated pharmacological data
2. **Target Coverage:** ~3,000 targets in 800+ families
3. **Ligand Coverage:** ~13,000 ligands including 2,200 approved drugs
4. **Interaction Data:** 170,000+ target-ligand interactions
5. **Affinity Data:** Quantitative pX values with ranges
6. **Nomenclature:** IUPHAR-standard target naming
7. **Quarterly Updates:** Regular data releases
8. **API Access:** REST API with JSON/XML responses

