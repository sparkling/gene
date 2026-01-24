---
id: schema-bindingdb
title: "BindingDB Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, binding-affinity, drug-target, ic50, ki, pharmacology]
---

# BindingDB - Binding Affinity Database Schema

**Document ID:** SCHEMA-BINDINGDB
**Version:** 2024
**Source Version:** Current (weekly updates)

---

## TL;DR

BindingDB provides quantitative binding affinity data (IC50, Ki, Kd, EC50) for drug-target interactions. The schema links compounds to protein targets with measured affinities, experimental conditions, and literature references, supporting drug discovery, SAR analysis, and computational model development.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Binding Measurements | 2,900,000+ | Activity data |
| Distinct Compounds | 1,300,000+ | Ligand collection |
| Protein Targets | ~9,400 | Target proteins |
| IC50 Measurements | 1,800,000+ | Inhibition data |
| Ki Measurements | 560,000+ | Binding constants |
| Publications | 40,000+ | Literature sources |
| Patents | 35,000+ | Patent sources |

---

## Entity Relationship Overview

```
Ligands (1) ←→ (many) Binding_Data (many) ←→ (1) Targets
   ↓                       ↓                      ↓
MonomerID              Ki/IC50/Kd            UniProt/PDB

Binding_Data (many) ←→ (1) Assays
                           ↓
                    Assay conditions

Binding_Data (many) ←→ (1) References
                           ↓
                    PubMed/Patent
```

---

## Core Tables/Entities

### ligands

**Description:** Small molecule compounds with structural data.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| monomer_id | integer | Yes | BindingDB ligand ID |
| name | string | No | Compound name |
| smiles | string | Yes | SMILES structure |
| inchi | string | No | InChI identifier |
| inchi_key | string | Yes | InChI Key |
| molecular_weight | decimal | No | MW in Daltons |
| molecular_formula | string | No | Chemical formula |
| alogp | decimal | No | Calculated LogP |
| psa | decimal | No | Polar surface area |
| hbd | integer | No | H-bond donors |
| hba | integer | No | H-bond acceptors |
| rotatable_bonds | integer | No | Rotatable bond count |

### targets

**Description:** Protein targets with identifiers.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| target_id | integer | Yes | BindingDB target ID |
| name | string | Yes | Target name |
| organism | string | Yes | Species |
| uniprot_id | string | No | UniProt accession |
| gene_symbol | string | No | Gene name |
| pdb_ids | array | No | PDB structure IDs |
| ec_number | string | No | EC classification |
| target_type | string | No | Protein, enzyme, GPCR |

### binding_data

**Description:** Measured binding affinities.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| bindingdb_id | integer | Yes | Primary identifier |
| monomer_id | integer | Yes | Foreign key to ligands |
| target_id | integer | Yes | Foreign key to targets |
| affinity_type | string | Yes | IC50, Ki, Kd, EC50 |
| affinity_value | decimal | Yes | Numeric value |
| affinity_unit | string | Yes | nM, uM, mM |
| affinity_relation | string | No | =, <, >, ~ |
| ph | decimal | No | Assay pH |
| temperature | decimal | No | Temperature (C) |
| assay_description | text | No | Assay method |
| reference_id | integer | No | Literature source |
| source_type | string | No | Publication, Patent |

### references

**Description:** Literature and patent sources.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| reference_id | integer | Yes | Primary identifier |
| pmid | integer | No | PubMed ID |
| doi | string | No | DOI |
| patent_number | string | No | Patent number |
| title | string | No | Article/patent title |
| authors | string | No | Author list |
| year | integer | No | Publication year |

---

## Affinity Types

| Type | Description | Typical Range | Use Case |
|------|-------------|---------------|----------|
| IC50 | Half-maximal inhibitory concentration | nM - uM | Enzyme inhibition |
| Ki | Inhibition constant | nM - uM | True binding affinity |
| Kd | Dissociation constant | nM - uM | Binding equilibrium |
| EC50 | Half-maximal effective concentration | nM - uM | Functional response |
| Kon | Association rate | M-1s-1 | Kinetics |
| Koff | Dissociation rate | s-1 | Kinetics |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /getLigandsByUniprots | GET | Ligands for UniProt targets |
| /getLigandsByPDBs | GET | Ligands for PDB structures |
| /getTargetByCompound | GET | Targets for compound |
| /getAffinitiesByLigand | GET | All affinities for ligand |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | TSV (tab-separated) |
| Alternative | SDF (structures) |
| API Response | XML, JSON |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "bindingdb_reactant_id": 50000001,
  "ligand": {
    "monomer_id": 12345,
    "name": "Example Inhibitor",
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "molecular_weight": 180.16
  },
  "target": {
    "name": "Cyclooxygenase-2",
    "uniprot_id": "P35354",
    "organism": "Homo sapiens",
    "pdb_ids": ["1CX2", "3LN1"]
  },
  "affinity": {
    "type": "IC50",
    "value": 0.6,
    "unit": "nM",
    "relation": "=",
    "ph": 7.4,
    "temperature": 25
  },
  "reference": {
    "pmid": 12345678,
    "year": 2020
  }
}
```

---

## TSV File Format

```
BindingDB_ID  Ligand_SMILES  Ligand_InChI_Key  Target_Name  UniProt_ID  Ki_nM  IC50_nM  Kd_nM  EC50_nM  pH  Temp  PubMed_ID
50000001      CC(=O)...      BSYNR...          COX-2        P35354      -      0.6      -      -        7.4 25    12345678
```

---

## Glossary

| Term | Definition |
|------|------------|
| MonomerID | BindingDB ligand identifier |
| IC50 | Concentration for 50% inhibition |
| Ki | Inhibitor binding constant |
| Kd | Equilibrium dissociation constant |
| pIC50 | -log10(IC50 in M) |
| SAR | Structure-Activity Relationship |

---

## References

1. BindingDB: https://www.bindingdb.org
2. Liu T, et al. (2025) Nucleic Acids Res. 53(D1):D1633-D1641
3. REST API: https://bindingdb.org/rest/
