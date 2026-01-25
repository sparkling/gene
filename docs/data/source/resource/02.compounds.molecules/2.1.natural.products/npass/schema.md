---
id: schema-npass
title: "NPASS Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, natural-products, bioactivity, species-source]
---

# NPASS - Natural Product Activity and Species Source Database Schema

**Document ID:** SCHEMA-NPASS
**Version:** 2.0
**Source Version:** 2023

---

## TL;DR

NPASS integrates natural product structures with quantitative bioactivity measurements and species source data. The schema links compounds to their biological targets (with IC50/Ki values) and to the organisms that produce them, enabling structure-activity-organism relationship studies.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Natural Products | 35,000+ | Compound collection |
| Quantitative Activities | 460,000+ | Bioactivity records |
| Source Species | 25,000+ | Organism database |
| Protein Targets | 5,000+ | Target annotations |
| Activity Types | IC50, EC50, Ki, etc. | Standardized measurements |

---

## Entity Relationship Overview

```
Compounds (1) ←→ (many) Activities (many) ←→ (1) Targets
     ↓                        ↓
  SMILES/InChI          IC50/Ki values

Compounds (1) ←→ (many) Compound_Species (many) ←→ (1) Species
                              ↓
                       NCBI Taxonomy links
```

---

## Core Tables/Entities

### compounds

**Description:** Natural product chemical structures with identifiers.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| npc_id | string | Yes | NPASS compound ID (NPC+digits) |
| name | string | No | Compound name |
| canonical_smiles | string | Yes | Canonical SMILES structure |
| inchi | string | No | InChI identifier |
| inchi_key | string | Yes | InChI Key for lookups |
| molecular_formula | string | No | Molecular formula |
| molecular_weight | decimal | No | Molecular weight (Da) |
| compound_class | string | No | Natural product class |

### activities

**Description:** Quantitative bioactivity measurements for compounds.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| activity_id | integer | Yes | Primary identifier |
| npc_id | string | Yes | Foreign key to compounds |
| target_id | integer | Yes | Foreign key to targets |
| activity_type | string | Yes | IC50, EC50, Ki, Kd, etc. |
| activity_value | decimal | Yes | Numeric value |
| activity_unit | string | Yes | nM, uM, mg/ml, etc. |
| activity_relation | string | No | =, <, >, ~, etc. |
| assay_type | string | No | Binding, functional, etc. |
| reference | string | No | PubMed ID or DOI |

### targets

**Description:** Biological targets of natural products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| target_id | integer | Yes | Primary identifier |
| target_name | string | Yes | Target name |
| uniprot_id | string | No | UniProt accession |
| gene_symbol | string | No | Gene symbol |
| organism | string | No | Target organism |
| target_type | string | No | Protein, enzyme, receptor |

### species

**Description:** Source organisms producing natural products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| species_id | integer | Yes | Primary identifier |
| scientific_name | string | Yes | Binomial name |
| ncbi_taxon_id | integer | No | NCBI Taxonomy ID |
| kingdom | string | No | Plant, fungi, bacteria, etc. |
| family | string | No | Taxonomic family |

### compound_species

**Description:** Links compounds to source species.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| npc_id | string | Yes | Foreign key to compounds |
| species_id | integer | Yes | Foreign key to species |
| part | string | No | Plant part if applicable |
| reference | string | No | Literature source |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /compound/{npc_id} | GET | Get compound by ID |
| /search/structure | POST | Substructure/similarity search |
| /search/activity | GET | Search by activity criteria |
| /species/{taxon_id}/compounds | GET | Compounds from species |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | TSV/CSV downloads |
| Alternative | SDF (structures with metadata) |
| Encoding | UTF-8 |
| Structure | SMILES, InChI |

---

## Sample Record

```json
{
  "npc_id": "NPC12345",
  "name": "Artemisinin",
  "smiles": "CC1CCC2C(C)C(=O)OC3OC4(C)CCC1C23OO4",
  "inchi_key": "BLUAFEHZUWYNDE-NNWCWBAJSA-N",
  "molecular_weight": 282.33,
  "activities": [
    {
      "target": "Plasmodium falciparum",
      "activity_type": "IC50",
      "value": 4.8,
      "unit": "nM",
      "assay_type": "Antimalarial"
    }
  ],
  "species": [
    {
      "name": "Artemisia annua",
      "ncbi_taxon_id": 35608,
      "kingdom": "Plantae"
    }
  ]
}
```

---

## Activity Value Standards

| Activity Type | Description | Typical Unit |
|---------------|-------------|--------------|
| IC50 | Half-maximal inhibitory concentration | nM, uM |
| EC50 | Half-maximal effective concentration | nM, uM |
| Ki | Inhibition constant | nM, uM |
| Kd | Dissociation constant | nM, uM |
| MIC | Minimum inhibitory concentration | ug/ml |

---

## Glossary

| Term | Definition |
|------|------------|
| NPC ID | NPASS compound identifier (NPC + digits) |
| IC50 | Concentration for 50% inhibition |
| Ki | Equilibrium inhibitor binding constant |
| Chemotaxonomy | Chemical patterns across taxonomic groups |

---

## References

1. Zeng X, et al. (2018) Nucleic Acids Res. 46(D1):D1217-D1222
2. NPASS Website: http://bidd.group/NPASS/
