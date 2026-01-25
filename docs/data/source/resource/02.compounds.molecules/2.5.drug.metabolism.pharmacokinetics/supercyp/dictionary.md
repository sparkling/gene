# SuperCYP - Data Dictionary

## Overview

This data dictionary documents the schema for SuperCYP cytochrome P450 enzyme-drug interaction database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | supercyp |
| **Name** | SuperCYP |
| **Parent** | 2.5.drug.metabolism.pharmacokinetics |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| compound_id | integer | 1:1 | Yes | Primary identifier | 1234 |
| name | string | 1:1 | Yes | Drug name | Ketoconazole |
| cas_number | string | 1:1 | No | CAS registry number | 65277-42-1 |
| drugbank_id | string | 1:1 | No | DrugBank identifier | DB01026 |
| pubchem_cid | integer | 1:1 | No | PubChem compound ID | 456201 |
| smiles | string | 1:1 | No | Chemical structure | CC(=O)Nc1ccc... |
| molecular_weight | decimal | 1:1 | No | MW in Daltons | 531.43 |
| drug_class | string | 1:1 | No | Therapeutic class | Antifungal |

### CYP Enzyme Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| cyp_id | string | 1:1 | Yes | CYP identifier | CYP3A4 |
| gene_symbol | string | 1:1 | Yes | Gene symbol | CYP3A4 |
| uniprot_id | string | 1:1 | No | UniProt accession | P08684 |
| chromosome | string | 1:1 | No | Chromosomal location | 7q22.1 |
| drug_metabolism_role | string | 1:1 | No | Percentage of drug metabolism | ~50% |
| substrate_count | integer | 1:1 | No | Number of substrates | 450 |
| inhibitor_count | integer | 1:1 | No | Number of inhibitors | 320 |

### CYP Interaction Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| interaction_id | integer | 1:1 | Yes | Primary identifier | 5678 |
| compound_id | integer | 1:1 | Yes | Foreign key to compounds | 1234 |
| cyp_id | string | 1:1 | Yes | Foreign key to cyp_enzymes | CYP3A4 |
| interaction_type | string | 1:1 | Yes | Type of interaction | substrate, inhibitor, inducer |
| strength | string | 1:1 | No | Interaction strength | strong, moderate, weak |
| ki_value | decimal | 1:1 | No | Inhibition constant (uM) | 0.015 |
| ic50_value | decimal | 1:1 | No | IC50 (uM) | 2.5 |
| induction_fold | decimal | 1:1 | No | Fold induction | 3.5 |
| evidence_level | string | 1:1 | No | Evidence type | in_vitro, in_vivo, clinical |
| reference_ids | array | 1:N | No | PubMed IDs | [12345678] |

### Metabolic Pathway Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| pathway_id | integer | 1:1 | Yes | Primary identifier | 101 |
| compound_id | integer | 1:1 | Yes | Parent compound | 1234 |
| cyp_id | string | 1:1 | Yes | Metabolizing enzyme | CYP3A4 |
| metabolite_name | string | 1:1 | No | Product name | N-dealkyl ketoconazole |
| reaction_type | string | 1:1 | No | Transformation type | N-dealkylation |
| site_of_metabolism | string | 1:1 | No | Position on molecule | N-1 imidazole |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Compound ID | Integer | 1234 | Internal identifier |
| CYP ID | CYP + family + subfamily + gene | CYP3A4 | Enzyme identifier |
| CAS Number | varies | 65277-42-1 | Chemical registry |
| DrugBank ID | DB + 5 digits | DB01026 | Drug database |
| PubChem CID | Integer | 456201 | Chemical database |
| UniProt ID | Accession | P08684 | Protein identifier |
| PubMed ID | PMID + digits | PMID:12345678 | Literature reference |

---

## Enumerations

### Interaction Types

| Type | Description | Clinical Relevance |
|------|-------------|-------------------|
| substrate | Metabolized by CYP | Affected by inhibitors/inducers |
| inhibitor | Reduces CYP activity | May increase substrate levels |
| inducer | Increases CYP expression | May decrease substrate levels |
| activator | Enhances CYP activity | Rare, research interest |

### Interaction Strength

| Strength | AUC Fold Change | Ki Range | Description |
|----------|-----------------|----------|-------------|
| strong | >= 5-fold | < 1 uM | Clinically significant |
| moderate | 2-5 fold | 1-10 uM | May require dose adjustment |
| weak | 1.25-2 fold | > 10 uM | Minimal clinical significance |

### Evidence Levels

| Level | Description |
|-------|-------------|
| in_vitro | Laboratory cell/microsome studies |
| in_vivo | Animal studies |
| clinical | Human clinical trials |

### CYP Isoforms

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

### Reaction Types

| Type | Description |
|------|-------------|
| Hydroxylation | Addition of hydroxyl group |
| N-dealkylation | Removal of alkyl from nitrogen |
| O-dealkylation | Removal of alkyl from oxygen |
| S-oxidation | Oxidation at sulfur |
| Epoxidation | Formation of epoxide |
| Reduction | Electron gain |

---

## Entity Relationships

### Compound to CYP Interactions
- **Cardinality:** 1:N
- **Description:** Each compound can have interactions with multiple CYP enzymes
- **Key Fields:** compound_id, cyp_id

### CYP Enzyme to Interactions
- **Cardinality:** 1:N
- **Description:** Each CYP enzyme has interactions with multiple compounds
- **Key Fields:** cyp_id, compound_id

### Interaction to References
- **Cardinality:** N:M
- **Description:** Interactions linked to multiple literature references
- **Key Fields:** interaction_id, reference_ids

### Compound to Metabolic Pathways
- **Cardinality:** 1:N
- **Description:** One compound may have multiple metabolic pathways
- **Key Fields:** compound_id, pathway_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CYP | Cytochrome P450 | Enzyme superfamily |
| CYP450 | Cytochrome P450 | Full name |
| Ki | Inhibition constant | Binding affinity |
| IC50 | Half-maximal inhibitory concentration | Potency measure |
| AUC | Area Under the Curve | Pharmacokinetic parameter |
| DDI | Drug-Drug Interaction | Clinical outcome |
| CCB | Calcium Channel Blocker | Drug class |
| SSRI | Selective Serotonin Reuptake Inhibitor | Drug class |
| NSAID | Non-Steroidal Anti-Inflammatory Drug | Drug class |
| PPI | Proton Pump Inhibitor | Drug class |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| DrugBank | DrugBank ID | Drug information |
| PubChem | CID | Chemical data |
| UniProt | Accession | Enzyme sequences |
| PubMed | PMID | Literature references |
| PharmGKB | Accession | Pharmacogenomics |
| KEGG | Drug ID | Pathway context |

---

## Data Quality Notes

1. **CYP Coverage:** 17 major human CYP isoforms
2. **Interaction Data:** 50,000+ CYP-compound relations
3. **Evidence Levels:** Data includes in vitro, in vivo, and clinical evidence
4. **Quantitative Data:** Ki and IC50 values when available
5. **Literature Support:** 2,500+ PubMed citations
6. **Clinical Focus:** Emphasis on drug metabolism prediction
7. **Update Frequency:** Periodic updates from literature

