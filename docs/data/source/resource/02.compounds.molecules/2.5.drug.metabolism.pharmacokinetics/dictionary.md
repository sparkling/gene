# 2.5 Drug Metabolism and Pharmacokinetics - Data Dictionary

## Overview

This data dictionary documents the unified schema for drug metabolism and pharmacokinetics (DMPK) data from two major databases: SuperCYP and SwissADME.

**Subcategory ID:** 2.5
**Subcategory Name:** Drug Metabolism and Pharmacokinetics
**Data Sources:** SuperCYP, SwissADME

---

## Unified Fields

### Core Identification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| compound_id | string/integer | 1:1, Optional | Primary identifier for compound | `1234`, `DB00945` | SuperCYP: drug_id, SwissADME: molecule_id |
| name | string | 1:1, Required | Drug/compound name | `Aspirin`, `Ibuprofen` | SuperCYP: drug_name, SwissADME: molecule_name |

### Structure Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| smiles | string | 1:1, Optional | Canonical SMILES structure | `CC(=O)Oc1ccccc1C(=O)O` | SuperCYP, SwissADME |
| molecular_weight | decimal | 1:1, Optional | Molecular weight in Daltons | `180.16` | SuperCYP, SwissADME |

### CYP Interaction Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| cyp_interactions | object[] | 1:N, Optional | CYP450 enzyme interactions | `[{"cyp_isoform": "CYP3A4", "interaction_type": "substrate", "strength": "strong"}]` | SuperCYP, SwissADME |

---

## Source-Specific Fields

### SuperCYP

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| cyp_id | string | 1:1, Optional | CYP isoform identifier | `CYP3A4`, `CYP2D6`, `CYP1A2`, `CYP2C9`, `CYP2C19` |
| interaction_type | string | 1:1, Optional | Type of CYP interaction | `substrate`, `inhibitor`, `inducer` |
| interaction_strength | string | 1:1, Optional | Interaction strength | `strong`, `moderate`, `weak` |
| ki_value | decimal | 1:1, Optional | Inhibition constant (uM) | `0.5`, `2.3` |
| ic50_value | decimal | 1:1, Optional | IC50 value (uM) | `1.2`, `15.0` |
| induction_fold | decimal | 1:1, Optional | Fold induction for inducers | `2.5`, `10.0` |
| evidence_level | string | 1:1, Optional | Level of evidence | `in_vitro`, `in_vivo`, `clinical` |
| drugbank_id | string | 1:1, Optional | DrugBank identifier for cross-reference | `DB00945` |

### SwissADME

#### Physicochemical Properties

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| inchi_key | string | 1:1, Optional | InChI Key identifier | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` |
| tpsa | decimal | 1:1, Optional | Topological Polar Surface Area (A^2) | `63.6` |
| hbd | integer | 1:1, Optional | Hydrogen bond donors | `1`, `2`, `3` |
| hba | integer | 1:1, Optional | Hydrogen bond acceptors | `3`, `4`, `5` |
| rotatable_bonds | integer | 1:1, Optional | Rotatable bond count | `2`, `5`, `8` |
| consensus_logp | decimal | 1:1, Optional | Average of 5 LogP methods (iLOGP, XLOGP3, WLOGP, MLOGP, SILICOS-IT) | `1.19` |

#### Drug-likeness Rules

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| lipinski_violations | integer | 1:1, Optional | Rule of 5 violations (0-4) | `0`, `1`, `2` |
| bioavailability_score | decimal | 1:1, Optional | Abbott bioavailability score (0-1) | `0.55`, `0.85` |

#### Pharmacokinetics Predictions

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| gi_absorption | string | 1:1, Optional | GI absorption prediction | `High`, `Low` |
| bbb_permeant | string | 1:1, Optional | Blood-brain barrier permeation | `Yes`, `No` |
| pgp_substrate | string | 1:1, Optional | P-glycoprotein substrate prediction | `Yes`, `No` |

#### CYP Inhibition Predictions

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| cyp1a2_inhibitor | string | 1:1, Optional | CYP1A2 inhibitor prediction | `Yes`, `No` |
| cyp2c19_inhibitor | string | 1:1, Optional | CYP2C19 inhibitor prediction | `Yes`, `No` |
| cyp2c9_inhibitor | string | 1:1, Optional | CYP2C9 inhibitor prediction | `Yes`, `No` |
| cyp2d6_inhibitor | string | 1:1, Optional | CYP2D6 inhibitor prediction | `Yes`, `No` |
| cyp3a4_inhibitor | string | 1:1, Optional | CYP3A4 inhibitor prediction | `Yes`, `No` |

#### Medicinal Chemistry

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| pains_alerts | integer | 1:1, Optional | PAINS (Pan-Assay Interference) pattern matches | `0`, `1`, `2` |
| synthetic_accessibility | decimal | 1:1, Optional | Synthetic accessibility score (1-10, lower=easier) | `2.5`, `4.8` |

---

## Metadata Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| _source.database | string | 1:1, Required | Source database name | `SuperCYP`, `SwissADME` |
| _source.version | string | 1:1, Optional | Database version | `2.1`, `2023` |
| _source.access_date | date | 1:1, Optional | Date data was accessed | `2026-01-24` |
| _source.original_id | string | 1:1, Optional | Original identifier from source | - |

---

## Field Mappings by Source

### SuperCYP to Unified Schema

| SuperCYP Field | Unified Field |
|----------------|---------------|
| drug_id | compound_id |
| drug_name | name |
| smiles | smiles |
| molecular_weight | molecular_weight |
| cyp_id | cyp_id |
| interaction_type | interaction_type |
| strength | interaction_strength |
| ki_value | ki_value |
| ic50_value | ic50_value |
| induction_fold | induction_fold |
| evidence_level | evidence_level |
| drugbank_id | drugbank_id |

### SwissADME to Unified Schema

| SwissADME Field | Unified Field |
|-----------------|---------------|
| molecule_name | name |
| smiles | smiles |
| inchikey | inchi_key |
| mw | molecular_weight |
| tpsa | tpsa |
| hbd | hbd |
| hba | hba |
| rotatable_bonds | rotatable_bonds |
| consensus_logp | consensus_logp |
| lipinski_violations | lipinski_violations |
| bioavailability_score | bioavailability_score |
| gi_absorption | gi_absorption |
| bbb_permeant | bbb_permeant |
| pgp_substrate | pgp_substrate |
| cyp1a2_inhibitor | cyp1a2_inhibitor |
| cyp2c19_inhibitor | cyp2c19_inhibitor |
| cyp2c9_inhibitor | cyp2c9_inhibitor |
| cyp2d6_inhibitor | cyp2d6_inhibitor |
| cyp3a4_inhibitor | cyp3a4_inhibitor |
| pains_alerts | pains_alerts |
| synthetic_accessibility | synthetic_accessibility |

---

## Data Quality Notes

- **Required Field:** Only `name` is strictly required across all sources
- **CYP Data:** SuperCYP provides experimental CYP interaction data; SwissADME provides computational predictions
- **Evidence Levels:** SuperCYP distinguishes between in vitro, in vivo, and clinical evidence
- **ADME Predictions:** SwissADME provides comprehensive ADME property predictions based on structure

### CYP450 Isoforms

| Isoform | Drug Metabolism Role | Common Substrates |
|---------|---------------------|-------------------|
| CYP3A4 | Most abundant; metabolizes ~50% of drugs | Midazolam, Nifedipine, Erythromycin |
| CYP2D6 | Metabolizes ~25% of drugs; highly polymorphic | Codeine, Dextromethorphan, Tamoxifen |
| CYP2C9 | Metabolizes ~15% of drugs | Warfarin, Phenytoin, NSAIDs |
| CYP2C19 | Metabolizes PPIs, clopidogrel | Omeprazole, Clopidogrel |
| CYP1A2 | Metabolizes caffeine, theophylline | Caffeine, Theophylline, Clozapine |
| CYP2B6 | Metabolizes bupropion, efavirenz | Bupropion, Efavirenz, Cyclophosphamide |
| CYP2E1 | Metabolizes ethanol, acetaminophen | Ethanol, Acetaminophen |

### Interaction Types (SuperCYP)

| Type | Description |
|------|-------------|
| substrate | Compound is metabolized by the CYP enzyme |
| inhibitor | Compound inhibits the CYP enzyme activity |
| inducer | Compound increases CYP enzyme expression |

### Interaction Strength (SuperCYP)

| Strength | Description |
|----------|-------------|
| strong | Significant clinical interaction expected |
| moderate | Moderate clinical interaction possible |
| weak | Minor clinical significance |

### Lipinski Rule of 5 (SwissADME)

| Rule | Threshold | Description |
|------|-----------|-------------|
| MW | <= 500 Da | Molecular weight |
| LogP | <= 5 | Lipophilicity |
| HBD | <= 5 | Hydrogen bond donors |
| HBA | <= 10 | Hydrogen bond acceptors |

Compounds with 0-1 violations are considered drug-like; >1 violation indicates potential oral bioavailability issues.

### Bioavailability Score (SwissADME)

| Score | Probability of F >= 10% in rat |
|-------|--------------------------------|
| 0.17 | Very low bioavailability expected |
| 0.55 | Moderate bioavailability expected |
| 0.85 | Good bioavailability expected |

### PAINS Alerts

PAINS (Pan-Assay Interference Compounds) are structural patterns known to cause false positives in high-throughput screens. Compounds with PAINS alerts should be evaluated carefully as potential artifacts.
