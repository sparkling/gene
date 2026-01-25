---
id: schema-gtopdb
title: "GtoPdb Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, pharmacology, iuphar, receptors, ion-channels, curated]
---

# GtoPdb - Guide to Pharmacology Schema

**Document ID:** SCHEMA-GTOPDB
**Version:** 2024.1
**Source Version:** 2024 (quarterly updates)

---

## TL;DR

GtoPdb provides expert-curated pharmacological data on drug targets (GPCRs, ion channels, enzymes) and their ligands with quantitative affinities in pX notation. The schema organizes targets by family with selective tool compounds, approved drugs, and IUPHAR-standard nomenclature.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Targets | ~3,000 | Target database |
| Target Families | 800+ | Classification |
| Ligands | ~13,000 | Compound collection |
| Approved Drugs | ~2,200 | Drug subset |
| Interactions | 170,000+ | Target-ligand pairs |
| References | 85,000+ | Literature citations |

---

## Entity Relationship Overview

```
Targets (1) ←→ (many) Interactions (many) ←→ (1) Ligands
   ↓                       ↓                      ↓
HGNC symbol           pKi/pIC50/pEC50      Ligand ID

Targets (1) ←→ (1) Target_Families
                      ↓
              GPCRs, Ion Channels, etc.

Ligands (1) ←→ (1) Ligand_Types
                    ↓
           Synthetic, Endogenous, Approved
```

---

## Core Tables/Entities

### targets

**Description:** Pharmacological targets with IUPHAR nomenclature.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| target_id | integer | Yes | GtoPdb target ID |
| name | string | Yes | Target name |
| abbreviation | string | No | Short name |
| systematic_name | string | No | IUPHAR systematic name |
| family_id | integer | Yes | Target family |
| type | string | Yes | receptor, ion_channel, etc. |
| hgnc_symbol | string | No | HGNC gene symbol |
| hgnc_id | integer | No | HGNC ID |
| uniprot_id | string | No | UniProt accession |
| species | string | Yes | Organism |

### ligands

**Description:** Compounds acting on targets.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ligand_id | integer | Yes | GtoPdb ligand ID |
| name | string | Yes | Ligand name |
| type | string | Yes | synthetic, metabolite, endogenous |
| approved | boolean | No | Approved drug status |
| radioactive | boolean | No | Radioligand |
| pubchem_sid | integer | No | PubChem SID |
| pubchem_cid | integer | No | PubChem CID |
| chembl_id | string | No | ChEMBL ID |
| inn | string | No | International nonproprietary name |
| smiles | string | No | Chemical structure |
| inchi | string | No | InChI |
| inchi_key | string | No | InChI Key |

### interactions

**Description:** Quantitative target-ligand interactions.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| interaction_id | integer | Yes | Primary identifier |
| target_id | integer | Yes | Foreign key to targets |
| ligand_id | integer | Yes | Foreign key to ligands |
| type | string | Yes | agonist, antagonist, etc. |
| action | string | No | Activation, Inhibition |
| selectivity | string | No | Selective, Non-selective |
| endogenous | boolean | No | Endogenous ligand flag |
| affinity_type | string | No | pKi, pIC50, pEC50, pKd |
| affinity_median | decimal | No | Median pX value |
| affinity_low | decimal | No | Low pX value |
| affinity_high | decimal | No | High pX value |
| species | string | No | Assay species |
| reference_ids | array | No | Literature references |

### target_families

**Description:** Hierarchical target classification.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| family_id | integer | Yes | Family identifier |
| name | string | Yes | Family name |
| parent_id | integer | No | Parent family |
| type | string | Yes | GPCR, ion_channel, etc. |

---

## Target Family Types

| Family | Description | Examples |
|--------|-------------|----------|
| GPCRs | G protein-coupled receptors | Adrenoceptors, opioid receptors |
| Ion channels | Voltage/ligand-gated channels | Nav, Cav, nAChR |
| NHRs | Nuclear hormone receptors | Estrogen receptor, PPAR |
| Kinases | Protein kinases | EGFR, BRAF |
| Enzymes | Other enzymes | COX, PDE |
| Transporters | SLC, ABC transporters | SERT, P-gp |
| Other proteins | Miscellaneous targets | Various |

---

## Affinity Notation (pX)

| Type | Description | Conversion |
|------|-------------|------------|
| pKi | -log10(Ki in M) | Ki (nM) = 10^(9-pKi) |
| pIC50 | -log10(IC50 in M) | IC50 (nM) = 10^(9-pIC50) |
| pEC50 | -log10(EC50 in M) | EC50 (nM) = 10^(9-pEC50) |
| pKd | -log10(Kd in M) | Kd (nM) = 10^(9-pKd) |

Example: pKi = 8 means Ki = 10 nM

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /targets | GET | List targets |
| /targets/{id} | GET | Get target details |
| /ligands | GET | List ligands |
| /ligands/{id} | GET | Get ligand details |
| /interactions | GET | Get interactions |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | PostgreSQL database dump |
| Alternative | CSV, TSV |
| API Response | JSON, XML |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "target": {
    "target_id": 290,
    "name": "Serine/threonine-protein kinase B-raf",
    "abbreviation": "BRAF",
    "type": "kinase",
    "hgnc_symbol": "BRAF",
    "uniprot_id": "P15056",
    "family": "Kinases"
  },
  "ligand": {
    "ligand_id": 5085,
    "name": "vemurafenib",
    "type": "synthetic",
    "approved": true,
    "chembl_id": "CHEMBL1667",
    "smiles": "CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2cc3ccc(Cl)cc3[nH]2)c1"
  },
  "interaction": {
    "type": "inhibitor",
    "action": "Inhibition",
    "selectivity": "Selective",
    "affinity_type": "pIC50",
    "affinity_median": 8.3,
    "species": "Human"
  }
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| IUPHAR | International Union of Basic and Clinical Pharmacology |
| BPS | British Pharmacological Society |
| pKi | Negative log10 of Ki in molar |
| Tool compound | Selective probe for target study |
| INN | International Nonproprietary Name |
| GPCR | G Protein-Coupled Receptor |

---

## References

1. GtoPdb: https://www.guidetopharmacology.org
2. Harding SD, et al. (2024) Nucleic Acids Res. 52(D1):D1438-D1449
3. API Documentation: https://www.guidetopharmacology.org/webServices.jsp
