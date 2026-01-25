# DrugBank - Data Dictionary

## Overview

This data dictionary documents the schema for DrugBank comprehensive drug database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | drugbank |
| **Name** | DrugBank |
| **Parent** | 2.2.pharmaceuticals |
| **Total Fields** | 50+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Drug Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| drugbank_id | string | 1:1 | Yes | DrugBank identifier | DB00945 |
| name | string | 1:1 | Yes | Drug name | Aspirin |
| type | string | 1:1 | Yes | Drug type | small molecule, biotech |
| description | text | 1:1 | No | Drug description | Full text description |
| cas_number | string | 1:1 | No | CAS registry number | 50-78-2 |
| unii | string | 1:1 | No | FDA UNII identifier | R16CO5Y76E |
| state | string | 1:1 | No | Physical state | solid, liquid, gas |
| groups | array | 1:N | Yes | Approval status | approved, investigational |

### Pharmacology

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| indication | text | 1:1 | No | Therapeutic indications | Pain, fever, inflammation |
| pharmacodynamics | text | 1:1 | No | PD description | Effects on body |
| mechanism_of_action | text | 1:1 | No | MOA description | COX inhibition |
| toxicity | text | 1:1 | No | Toxicity information | Side effects |
| metabolism | text | 1:1 | No | Metabolic pathways | Hepatic metabolism |
| absorption | text | 1:1 | No | Absorption profile | Oral bioavailability |
| half_life | string | 1:1 | No | Elimination half-life | 15-20 minutes |
| protein_binding | string | 1:1 | No | Protein binding | 99% |
| volume_of_distribution | string | 1:1 | No | Vd value | 0.15 L/kg |
| clearance | string | 1:1 | No | Clearance value | 650 mL/min |

### Chemical Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| molecular_formula | string | 1:1 | No | Chemical formula | C9H8O4 |
| molecular_weight | decimal | 1:1 | No | MW in Daltons | 180.16 |
| monoisotopic_mass | decimal | 1:1 | No | Exact mass | 180.042 |
| smiles | string | 1:1 | No | SMILES string | CC(=O)Oc1ccccc1C(=O)O |
| inchi | string | 1:1 | No | InChI identifier | InChI=1S/C9H8O4/... |
| inchikey | string | 1:1 | No | InChI Key | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| logp | decimal | 1:1 | No | Partition coefficient | 1.2 |
| pka | string | 1:1 | No | Ionization constant | 3.5 |

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| target_id | string | 1:1 | Yes | Target identifier | BE0000017 |
| target_name | string | 1:1 | Yes | Target name | Prostaglandin G/H synthase 1 |
| organism | string | 1:1 | Yes | Species | Humans |
| gene_name | string | 1:1 | No | Gene symbol | PTGS1 |
| uniprot_id | string | 1:1 | No | UniProt accession | P23219 |
| actions | array | 1:N | Yes | Action types | inhibitor, agonist |

### Drug Interactions

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| interacting_drug_id | string | 1:1 | Yes | Interacting drug ID | DB00001 |
| description | text | 1:1 | Yes | Interaction description | Increased bleeding risk |
| severity | string | 1:1 | No | Severity level | major, moderate, minor |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| DrugBank ID | DB + 5 digits | DB00945 | Primary identifier |
| CAS Number | varies | 50-78-2 | Chemical registry |
| UNII | 10 characters | R16CO5Y76E | FDA substance ID |
| InChI Key | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | Structure hash |
| UniProt | Accession | P23219 | Protein target |

---

## Enumerations

### Drug Groups

| Group | Description |
|-------|-------------|
| approved | FDA/EMA approved |
| investigational | Clinical trials |
| experimental | Preclinical |
| withdrawn | Removed from market |
| nutraceutical | Nutritional products |
| vet_approved | Veterinary use |
| illicit | Controlled substances |

### Drug Types

| Type | Description |
|------|-------------|
| small molecule | Chemical compound |
| biotech | Biologic product |

### Target Actions

| Action | Description |
|--------|-------------|
| inhibitor | Blocks target activity |
| agonist | Activates receptor |
| antagonist | Blocks receptor |
| substrate | Metabolized by enzyme |
| inducer | Increases expression |
| binder | Binds without action |
| modulator | Modifies activity |
| activator | Enhances activity |

### Interaction Severity

| Severity | Description |
|----------|-------------|
| major | Potentially life-threatening |
| moderate | May require intervention |
| minor | Minimal clinical significance |

---

## Entity Relationships

### Drug to Targets
- **Cardinality:** N:M
- **Description:** Drugs interact with multiple protein targets
- **Key Fields:** drugbank_id, target_id, actions

### Drug to Interactions
- **Cardinality:** N:M
- **Description:** Drugs may interact with multiple other drugs
- **Key Fields:** drugbank_id, interacting_drug_id

### Drug to Pathways
- **Cardinality:** N:M
- **Description:** Drugs participate in metabolic/signaling pathways
- **Key Fields:** drugbank_id, smpdb_id

### Drug to Chemical Properties
- **Cardinality:** 1:1
- **Description:** Each drug has one set of chemical properties
- **Key Fields:** drugbank_id

### Target to Sequences
- **Cardinality:** 1:1
- **Description:** Each target has one protein sequence
- **Key Fields:** target_id, uniprot_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DDI | Drug-Drug Interaction | Interaction data |
| MOA | Mechanism of Action | Drug-target interaction |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity | PK/PD properties |
| PD | Pharmacodynamics | Drug effects |
| PK | Pharmacokinetics | Drug metabolism |
| Vd | Volume of Distribution | PK parameter |
| CAS | Chemical Abstracts Service | Chemical registry |
| UNII | Unique Ingredient Identifier | FDA identifier |
| SMPDB | Small Molecule Pathway Database | Pathway reference |
| SMILES | Simplified Molecular Input Line Entry System | Structure |
| InChI | International Chemical Identifier | IUPAC standard |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UniProt | Accession | Protein targets |
| PDB | PDB ID | 3D structures |
| ChEMBL | ChEMBL ID | Bioactivity data |
| PubChem | CID | Chemical data |
| KEGG | Drug ID | Pathway information |
| SMPDB | Pathway ID | Metabolic pathways |
| PharmGKB | Accession | Pharmacogenomics |
| RxNorm | RXCUI | Drug terminology |
| DailyMed | Set ID | Drug labels |

---

## Data Quality Notes

1. **Comprehensive Coverage:** 16,000+ drug entries including approved, investigational, and experimental
2. **Target Annotation:** 5,200+ protein targets with UniProt cross-references
3. **DDI Database:** 380,000+ drug-drug interactions with severity ratings
4. **ADMET Data:** Detailed pharmacokinetic and pharmacodynamic information
5. **License:** CC BY-NC 4.0 (non-commercial) or commercial license
6. **Update Frequency:** Quarterly updates with new drug entries
7. **Data Format:** Primary XML with JSON API available
