# HIT 2.0 - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | hit.2.0 |
| **Name** | HIT 2.0 - Herbal Ingredients' Targets |
| **Total Fields** | 28 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| ingredient_id | String | Yes | HIT ingredient identifier | HIT001234 |
| ingredient_name | String | No | Chemical compound name | Curcumin |
| target_id | String | Yes | UniProt protein accession | P35354 |
| gene_symbol | String | No | HGNC gene symbol | PTGS2 |
| target_name | String | No | Full protein name | Prostaglandin G/H synthase 2 |
| interaction_id | String | No | HIT interaction identifier | HIT_INT_001234 |

---

## Evidence Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| evidence.evidence_type | Enum | Type of experimental validation |
| evidence.assay_type | String | Specific assay used |
| evidence.pubmed_id | Integer | PubMed reference |
| evidence.citation | String | Full literature citation |

---

## Binding Affinity Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| binding_affinity.affinity_type | Enum | Measurement type | IC50, Ki, Kd |
| binding_affinity.affinity_value | Number | Quantitative value | 2.5 |
| binding_affinity.affinity_unit | Enum | Unit | nM, uM, mM |
| binding_affinity.activity_type | Enum | Activity classification | inhibitor |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Ingredient ID | HIT###### | HIT compound ID | HIT001234 |
| Interaction ID | HIT_INT_###### | Interaction ID | HIT_INT_001234 |
| UniProt | [A-Z][0-9]{5} | Protein accession | P35354 |
| PubChem CID | Numeric | Compound ID | 969516 |
| PubMed ID | Numeric | Literature reference | 12345678 |

---

## Enumerations

### Evidence Types

| Value | Description |
|-------|-------------|
| binding_assay | Direct binding measurement |
| enzyme_inhibition | Enzyme activity assay |
| cell_based | Cellular functional assay |
| in_vivo | Animal model study |
| clinical | Human clinical study |

### Affinity Types

| Value | Description |
|-------|-------------|
| IC50 | Half-maximal inhibitory concentration |
| Ki | Inhibition constant |
| Kd | Dissociation constant |
| EC50 | Half-maximal effective concentration |

### Activity Types

| Value | Description |
|-------|-------------|
| inhibitor | Inhibits target activity |
| agonist | Activates target |
| antagonist | Blocks target activation |
| modulator | Modulates target activity |
| binder | Binds without functional effect |

### Traditional Systems

| Value | Description |
|-------|-------------|
| TCM | Traditional Chinese Medicine |
| Ayurveda | Indian traditional medicine |
| Kampo | Japanese traditional medicine |
| Western | Western herbal medicine |
| African | African traditional medicine |

---

## Entity Relationships

### Ingredient to Source Herbs
- **Cardinality:** N:M
- **Description:** Compounds found in multiple herbs
- **Key Fields:** ingredient_id, source_herbs

### Ingredient to Target
- **Cardinality:** N:M
- **Description:** Experimentally validated interactions
- **Key Fields:** ingredient_id, target_id

### Interaction to Evidence
- **Cardinality:** 1:N
- **Description:** Interactions supported by evidence
- **Key Fields:** interaction_id, evidence

---

## Cross-References

| Database | Purpose |
|----------|---------|
| BATMAN-TCM | TCM validation |
| KampoDB | Kampo validation |
| IMPPAT | Ayurveda validation |
| SymMap | Symptom mapping |
| ChEMBL | Activity data |
| DrugBank | Drug information |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HIT | Herbal Ingredients' Targets | Database name |
| TTI | Target-TCM Interaction | Interaction type |
| IC50 | Inhibitory Concentration 50% | Affinity metric |
| Ki | Inhibition constant | Affinity metric |
| Kd | Dissociation constant | Affinity metric |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
