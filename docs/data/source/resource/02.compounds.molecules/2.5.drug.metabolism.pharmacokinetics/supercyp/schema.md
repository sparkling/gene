---
id: schema-supercyp
title: "SuperCYP Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, cytochrome-p450, drug-metabolism, cyp-enzymes, drug-interactions]
---

# SuperCYP - Cytochrome P450 Database Schema

**Document ID:** SCHEMA-SUPERCYP
**Version:** 2.0
**Source Version:** 2023

---

## TL;DR

SuperCYP provides comprehensive CYP450 enzyme-drug interaction data covering substrates, inhibitors, and inducers. The schema links compounds to specific CYP isoforms with interaction type and literature evidence, enabling drug-drug interaction prediction and metabolic liability assessment.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Compounds | 6,000+ | Drug collection |
| CYP-Compound Relations | 50,000+ | Interaction records |
| CYP Isoforms | 17 | Major human CYPs |
| Substrates | 3,500+ | Metabolized drugs |
| Inhibitors | 4,500+ | CYP inhibitors |
| Inducers | 1,200+ | CYP inducers |
| Literature References | 2,500+ | PubMed citations |

---

## Entity Relationship Overview

```
Compounds (1) ←→ (many) CYP_Interactions (many) ←→ (1) CYP_Enzymes
     ↓                        ↓                         ↓
 Drug/CAS ID           Interaction type            CYP1A2, CYP3A4

CYP_Interactions (many) ←→ (many) References
                              ↓
                         PubMed IDs

CYP_Enzymes (1) ←→ (1) Gene_Info
                         ↓
                    Gene symbol, UniProt
```

---

## Core Tables/Entities

### compounds

**Description:** Drug compounds with identifiers.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| compound_id | integer | Yes | Primary identifier |
| name | string | Yes | Drug name |
| cas_number | string | No | CAS registry number |
| drugbank_id | string | No | DrugBank identifier |
| pubchem_cid | integer | No | PubChem compound ID |
| smiles | string | No | Chemical structure |
| molecular_weight | decimal | No | MW in Daltons |
| drug_class | string | No | Therapeutic class |

### cyp_enzymes

**Description:** Cytochrome P450 enzyme isoforms.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| cyp_id | string | Yes | CYP identifier (e.g., CYP3A4) |
| gene_symbol | string | Yes | Gene symbol |
| uniprot_id | string | No | UniProt accession |
| chromosome | string | No | Chromosomal location |
| drug_metabolism_role | string | No | Percentage of drug metabolism |
| substrate_count | integer | No | Number of substrates |
| inhibitor_count | integer | No | Number of inhibitors |

### cyp_interactions

**Description:** Drug-CYP enzyme interactions.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| interaction_id | integer | Yes | Primary identifier |
| compound_id | integer | Yes | Foreign key to compounds |
| cyp_id | string | Yes | Foreign key to cyp_enzymes |
| interaction_type | string | Yes | substrate, inhibitor, inducer |
| strength | string | No | strong, moderate, weak |
| ki_value | decimal | No | Inhibition constant (uM) |
| ic50_value | decimal | No | IC50 (uM) |
| induction_fold | decimal | No | Fold induction |
| evidence_level | string | No | in_vitro, in_vivo, clinical |
| reference_ids | array | No | PubMed IDs |

### metabolic_pathways

**Description:** Metabolic transformations by CYP enzymes.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| pathway_id | integer | Yes | Primary identifier |
| compound_id | integer | Yes | Parent compound |
| cyp_id | string | Yes | Metabolizing enzyme |
| metabolite_name | string | No | Product name |
| reaction_type | string | No | Oxidation, hydroxylation |
| site_of_metabolism | string | No | Position on molecule |

---

## CYP Isoforms Covered

| CYP | Drug Metabolism % | Key Substrates |
|-----|-------------------|----------------|
| CYP3A4 | ~50% | Statins, CCBs, macrolides |
| CYP3A5 | ~50% | Similar to CYP3A4 |
| CYP2D6 | ~25% | Codeine, SSRIs, beta-blockers |
| CYP2C9 | ~15% | Warfarin, NSAIDs |
| CYP2C19 | ~10% | PPIs, clopidogrel |
| CYP1A2 | ~5% | Caffeine, theophylline |
| CYP2B6 | ~3% | Efavirenz, ketamine |
| CYP2E1 | ~2% | Ethanol, acetaminophen |

---

## Interaction Types

| Type | Description | Clinical Relevance |
|------|-------------|-------------------|
| Substrate | Metabolized by CYP | Affected by inhibitors |
| Inhibitor | Reduces CYP activity | May increase substrate levels |
| Inducer | Increases CYP expression | May decrease substrate levels |
| Activator | Enhances CYP activity | Rare, research interest |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | Web interface |
| Alternative | Download available |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "compound": {
    "name": "Ketoconazole",
    "cas_number": "65277-42-1",
    "drug_class": "Antifungal"
  },
  "cyp_interactions": [
    {
      "cyp_id": "CYP3A4",
      "interaction_type": "inhibitor",
      "strength": "strong",
      "ki_value": 0.015,
      "evidence_level": "clinical",
      "references": ["PMID:12345678"]
    },
    {
      "cyp_id": "CYP2C9",
      "interaction_type": "inhibitor",
      "strength": "moderate",
      "ic50_value": 2.5,
      "evidence_level": "in_vitro"
    }
  ]
}
```

---

## Inhibitor Strength Classification

| Strength | AUC Fold Change | Ki Range |
|----------|-----------------|----------|
| Strong | >= 5-fold | < 1 uM |
| Moderate | 2-5 fold | 1-10 uM |
| Weak | 1.25-2 fold | > 10 uM |

---

## Glossary

| Term | Definition |
|------|------------|
| CYP450 | Cytochrome P450 enzyme family |
| Ki | Inhibition constant |
| IC50 | Half-maximal inhibitory concentration |
| Substrate | Compound metabolized by enzyme |
| Inhibitor | Compound blocking enzyme activity |
| Inducer | Compound increasing enzyme expression |

---

## References

1. SuperCYP: http://bioinformatics.charite.de/supercyp
2. Preissner S, et al. (2010) Nucleic Acids Res. 38(Database issue):D237-D243
