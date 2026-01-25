# PubChem - Data Dictionary

## Overview

This data dictionary documents the schema for PubChem compound database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pubchem |
| **Name** | PubChem |
| **Parent** | 2.6.chemical.ontology.classification |
| **Total Fields** | 50+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Record

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| cid | integer | 1:1 | Yes | Compound Identifier | 2244 |
| iupac_name | string | 1:1 | No | IUPAC systematic name | 2-acetyloxybenzoic acid |
| molecular_formula | string | 1:1 | Yes | Chemical formula | C9H8O4 |
| molecular_weight | decimal | 1:1 | Yes | MW in Daltons | 180.16 |
| canonical_smiles | string | 1:1 | Yes | Canonical SMILES | CC(=O)OC1=CC=CC=C1C(=O)O |
| isomeric_smiles | string | 1:1 | No | Stereochemistry-aware SMILES | - |
| inchi | string | 1:1 | No | InChI identifier | InChI=1S/C9H8O4/... |
| inchi_key | string | 1:1 | Yes | InChI Key hash | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| synonyms | array | 1:N | No | Alternative names | [Aspirin, acetylsalicylic acid] |

### Physicochemical Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| XLogP | decimal | 1:1 | No | Partition coefficient | 1.31 |
| ExactMass | decimal | 1:1 | No | Exact mass | 180.042 |
| MonoisotopicMass | decimal | 1:1 | No | Monoisotopic mass | 180.042 |
| TPSA | decimal | 1:1 | No | Topological polar surface area | 63.6 |
| Complexity | decimal | 1:1 | No | Molecular complexity | 212 |
| Charge | integer | 1:1 | No | Formal charge | 0 |
| HBondDonorCount | integer | 1:1 | No | H-bond donors | 1 |
| HBondAcceptorCount | integer | 1:1 | No | H-bond acceptors | 4 |
| RotatableBondCount | integer | 1:1 | No | Rotatable bonds | 3 |
| HeavyAtomCount | integer | 1:1 | No | Non-hydrogen atoms | 13 |

### Stereochemistry Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| IsotopeAtomCount | integer | 1:1 | No | Isotope-labeled atoms | 0 |
| AtomStereoCount | integer | 1:1 | No | Stereocenters | 0 |
| DefinedAtomStereoCount | integer | 1:1 | No | Defined stereocenters | 0 |
| UndefinedAtomStereoCount | integer | 1:1 | No | Undefined stereocenters | 0 |
| BondStereoCount | integer | 1:1 | No | E/Z bonds | 0 |
| CovalentUnitCount | integer | 1:1 | No | Covalent units | 1 |

### 3D Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Volume3D | decimal | 1:1 | No | Molecular volume | - |
| FeatureCount3D | integer | 1:1 | No | 3D features | - |
| FeatureAcceptorCount3D | integer | 1:1 | No | 3D acceptors | - |
| FeatureDonorCount3D | integer | 1:1 | No | 3D donors | - |
| ConformerCount3D | integer | 1:1 | No | Conformer count | - |
| Fingerprint2D | string | 1:1 | No | Structural fingerprint | - |

### Bioactivity Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| aid | integer | 1:1 | Yes | Assay Identifier | 504327 |
| source_name | string | 1:1 | Yes | Data source | ChEMBL |
| assay_type | string | 1:1 | No | Assay type | Binding |
| target_name | string | 1:1 | No | Target protein name | COX-2 |
| target_gene_id | integer | 1:1 | No | Entrez Gene ID | 5743 |
| activity | string | 1:1 | No | Activity outcome | Active |
| value | decimal | 1:1 | No | Activity value | 0.1 |
| unit | string | 1:1 | No | Value unit | uM |

### Safety Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ghs_hazards | array | 1:N | No | GHS hazard codes | [H302, H315] |
| precautionary_statements | array | 1:N | No | Safety precautions | [P264, P280] |
| ld50_species | string | 1:1 | No | LD50 test species | Rat |
| ld50_route | string | 1:1 | No | Administration route | Oral |
| ld50_value | decimal | 1:1 | No | LD50 value | 200 |
| ld50_unit | string | 1:1 | No | LD50 unit | mg/kg |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| CID | Integer | 2244 | Compound Identifier |
| SID | Integer | 123456 | Substance Identifier |
| AID | Integer | 504327 | Assay Identifier |
| InChIKey | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | Structure hash |
| CAS Number | varies | 50-78-2 | Chemical registry |

---

## Enumerations

### Activity Outcomes

| Outcome | Description |
|---------|-------------|
| Active | Compound shows activity |
| Inactive | No significant activity |
| Inconclusive | Results unclear |
| Unspecified | Activity not specified |
| Probe | Active as chemical probe |

### Assay Types

| Type | Description |
|------|-------------|
| Binding | Receptor binding assay |
| Functional | Functional activity assay |
| ADMET | ADME/Toxicity assay |
| Cell-based | Cellular assay |
| Biochemical | Biochemical assay |

### GHS Hazard Classes

| Code | Description |
|------|-------------|
| H300-H319 | Acute toxicity |
| H320-H373 | Health hazards |
| H400-H413 | Environmental hazards |

### Query Parameters

| Parameter | Description |
|-----------|-------------|
| exact | Exact structure match |
| similarity | Similarity search (threshold) |
| substructure | Substructure search |
| superstructure | Superstructure search |

---

## Entity Relationships

### Compound to Properties
- **Cardinality:** 1:1
- **Description:** Each compound has one property set
- **Key Fields:** cid

### Compound to Synonyms
- **Cardinality:** 1:N
- **Description:** Each compound has multiple names
- **Key Fields:** cid, synonyms

### Compound to Bioassays
- **Cardinality:** N:M
- **Description:** Compounds tested in multiple assays
- **Key Fields:** cid, aid

### Compound to Literature
- **Cardinality:** 1:N
- **Description:** Linked to PubMed references
- **Key Fields:** cid, pmid

### Compound to Patents
- **Cardinality:** 1:N
- **Description:** Linked to patent documents
- **Key Fields:** cid, patent_id

### Compound to Conformers
- **Cardinality:** 1:N
- **Description:** 3D conformer structures
- **Key Fields:** cid, conformer_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CID | Compound Identifier | Unique compound ID |
| SID | Substance Identifier | Submitted substance ID |
| AID | Assay Identifier | Bioassay ID |
| SMILES | Simplified Molecular-Input Line-Entry System | Structure notation |
| InChI | International Chemical Identifier | IUPAC standard |
| TPSA | Topological Polar Surface Area | Absorption predictor |
| GHS | Globally Harmonized System | Hazard classification |
| LD50 | Lethal Dose 50% | Toxicity measure |
| PUG | Power User Gateway | PubChem API system |
| MW | Molecular Weight | Mass in Daltons |
| XLogP | Computed LogP | Lipophilicity estimate |
| HBA | Hydrogen Bond Acceptor | Lipinski parameter |
| HBD | Hydrogen Bond Donor | Lipinski parameter |
| FTP | File Transfer Protocol | Bulk download |
| NCBI | National Center for Biotechnology Information | Host institution |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| ChEMBL | ChEMBL ID | Bioactivity data |
| DrugBank | DrugBank ID | Drug information |
| UniProt | Accession | Target proteins |
| PubMed | PMID | Literature |
| ChEBI | ChEBI ID | Chemical ontology |
| KEGG | Compound ID | Pathway context |

---

## Data Quality Notes

1. **Comprehensive Coverage:** 115M+ compounds
2. **Bioassay Data:** 2M+ assays, 500M+ activity records
3. **Literature Links:** 100M+ publication references
4. **Patent Coverage:** 50M+ patent references
5. **3D Structures:** 50M+ computed conformers
6. **Public Domain:** CC0 license
7. **Daily Updates:** Continuous incremental updates
8. **API Access:** RESTful PUG API with rate limits

