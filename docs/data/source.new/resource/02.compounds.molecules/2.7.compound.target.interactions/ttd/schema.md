---
id: schema-ttd
title: "TTD Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, therapeutic-targets, drug-targets, druggability, clinical-trials]
---

# TTD - Therapeutic Target Database Schema

**Document ID:** SCHEMA-TTD
**Version:** 2024
**Source Version:** 2024

---

## TL;DR

TTD provides therapeutic target information organized by development status (successful, clinical trial, preclinical). The schema links targets to drugs, diseases (ICD-11), and pathways with druggability assessments across molecular, system, and cellular perspectives.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Total Targets | 3,131 | Target database |
| Successful Targets | 426 | FDA-approved drug |
| Clinical Trial Targets | 1,014 | In development |
| Total Drugs | 39,862 | Drug collection |
| Approved Drugs | 2,895 | Marketed drugs |
| Clinical Trial Drugs | 11,796 | In clinical development |

---

## Entity Relationship Overview

```
Targets (1) ←→ (many) Target_Drugs (many) ←→ (1) Drugs
    ↓                       ↓                    ↓
TTD Target ID         Development status    TTD Drug ID

Targets (1) ←→ (many) Target_Diseases (many) ←→ (1) Diseases
                            ↓
                      ICD-11 codes

Targets (1) ←→ (many) Target_Pathways
                           ↓
                    KEGG, Reactome
```

---

## Core Tables/Entities

### targets

**Description:** Therapeutic protein and nucleic acid targets.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ttd_target_id | string | Yes | TTD identifier (TTDT+5 digits) |
| name | string | Yes | Target name |
| uniprot_id | string | No | UniProt accession |
| gene_name | string | No | Gene symbol |
| target_type | string | Yes | Protein, DNA, RNA |
| biochemical_class | string | No | Enzyme, receptor, etc. |
| ecl_classification | string | No | EC number if enzyme |
| target_validation | string | No | Successful, clinical, literature |
| organism | string | No | Species |

### drugs

**Description:** Drugs and compounds targeting entries.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ttd_drug_id | string | Yes | TTD drug ID (6 chars) |
| name | string | Yes | Drug name |
| cas_number | string | No | CAS registry |
| pubchem_cid | integer | No | PubChem compound |
| chembl_id | string | No | ChEMBL ID |
| drugbank_id | string | No | DrugBank ID |
| drug_class | string | No | Small molecule, antibody, etc. |
| highest_status | string | No | Approved, clinical, discontinued |
| company | string | No | Developer company |

### target_drugs

**Description:** Target-drug relationships.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ttd_target_id | string | Yes | Foreign key to targets |
| ttd_drug_id | string | Yes | Foreign key to drugs |
| activity_type | string | No | Inhibitor, agonist, etc. |
| mechanism | string | No | Mechanism of action |
| development_status | string | Yes | Approved, Phase I-III, etc. |
| indication | string | No | Therapeutic indication |

### diseases

**Description:** Disease associations with ICD-11 coding.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| disease_id | string | Yes | TTD disease ID |
| name | string | Yes | Disease name |
| icd_11_code | string | No | ICD-11 classification |
| therapeutic_area | string | No | Oncology, CNS, etc. |

### pathways

**Description:** Biological pathway associations.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| pathway_id | string | Yes | Pathway identifier |
| pathway_name | string | Yes | Pathway name |
| source | string | Yes | KEGG, Reactome, WikiPathways |
| external_id | string | No | Source pathway ID |

---

## Target Validation Status

| Status | Count | Description |
|--------|-------|-------------|
| Successful | 426 | Approved drug exists |
| Clinical Trial | 1,014 | Phase I-III drugs |
| Preclinical/Patented | 212 | Early development |
| Literature-reported | 1,479 | Research targets |

---

## Drug Development Status

| Status | Description |
|--------|-------------|
| Approved | FDA/EMA approved |
| Phase III | Late-stage trials |
| Phase II | Efficacy trials |
| Phase I | Safety trials |
| Preclinical | Pre-IND |
| Discontinued | Development stopped |
| Withdrawn | Removed from market |

---

## Drug Types

| Type | Description | Examples |
|------|-------------|----------|
| Small Molecule | Chemical compounds | Imatinib |
| Antibody | Monoclonal antibodies | Trastuzumab |
| ADC | Antibody-drug conjugate | T-DM1 |
| Peptide | Therapeutic peptides | Exenatide |
| ASO | Antisense oligonucleotide | Nusinersen |
| siRNA | Small interfering RNA | Patisiran |
| mRNA | Messenger RNA | COVID vaccines |
| Cell Therapy | CAR-T, stem cells | Axicabtagene |
| Gene Therapy | Gene delivery | Zolgensma |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | TXT (tab-separated) |
| Alternative | Web search |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "target": {
    "ttd_target_id": "TTDT00001",
    "name": "Epidermal growth factor receptor",
    "uniprot_id": "P00533",
    "gene_name": "EGFR",
    "target_type": "Successful target",
    "biochemical_class": "Kinase"
  },
  "drugs": [
    {
      "ttd_drug_id": "D0A9YA",
      "name": "Erlotinib",
      "highest_status": "Approved",
      "drug_class": "Small Molecule",
      "mechanism": "Inhibitor"
    },
    {
      "ttd_drug_id": "D01B2X",
      "name": "Cetuximab",
      "highest_status": "Approved",
      "drug_class": "Antibody",
      "mechanism": "Antagonist"
    }
  ],
  "diseases": [
    {
      "name": "Non-small cell lung cancer",
      "icd_11_code": "2C25"
    }
  ]
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| TTD Target ID | TTDT + 5 digits identifier |
| TTD Drug ID | 6-character drug identifier |
| Successful target | Target with approved drug |
| Clinical target | Target in clinical development |
| ICD-11 | International Classification of Diseases 11 |
| Druggability | Likelihood of target being druggable |

---

## References

1. TTD: https://idrblab.net/ttd/
2. Zhou Y, et al. (2024) Nucleic Acids Res. 52(D1):D1465-D1477
