# SwissADME - Data Dictionary

## Overview

This data dictionary documents the schema for SwissADME ADME property prediction web service.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | swissadme |
| **Name** | SwissADME |
| **Parent** | 2.5.drug.metabolism.pharmacokinetics |
| **Total Fields** | 45+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Molecular Identity

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Molecule | string | 1:1 | Yes | Molecule name or index | Aspirin |
| Canonical_SMILES | string | 1:1 | Yes | Standardized SMILES | CC(=O)OC1=CC=CC=C1C(=O)O |
| InChIKey | string | 1:1 | Yes | InChI key identifier | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| Formula | string | 1:1 | Yes | Molecular formula | C9H8O4 |

### Physicochemical Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| MW | decimal | 1:1 | Yes | Molecular weight (Da) | 180.16 |
| Heavy_Atoms | integer | 1:1 | Yes | Non-hydrogen atom count | 13 |
| Aromatic_Heavy_Atoms | integer | 1:1 | Yes | Aromatic heavy atoms | 6 |
| Fraction_Csp3 | decimal | 1:1 | Yes | Fraction sp3 carbons (0-1) | 0.11 |
| Rotatable_Bonds | integer | 1:1 | Yes | Rotatable bond count | 3 |
| HBD | integer | 1:1 | Yes | H-bond donors | 1 |
| HBA | integer | 1:1 | Yes | H-bond acceptors | 4 |
| MR | decimal | 1:1 | Yes | Molar refractivity | 44.46 |
| TPSA | decimal | 1:1 | Yes | Topological polar surface area (A^2) | 63.60 |

### Lipophilicity Predictions

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| iLOGP | decimal | 1:1 | Yes | In-house physics-based LogP | 1.31 |
| XLOGP3 | decimal | 1:1 | Yes | Atomistic method LogP | 1.24 |
| WLOGP | decimal | 1:1 | Yes | Wildman-Crippen method | 1.02 |
| MLOGP | decimal | 1:1 | Yes | Moriguchi method | 0.88 |
| SILICOS_IT_LogP | decimal | 1:1 | Yes | SILICOS-IT method | 1.51 |
| Consensus_LogP | decimal | 1:1 | Yes | Average of 5 methods | 1.19 |

### Water Solubility Predictions

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ESOL_LogS | decimal | 1:1 | Yes | ESOL method log S | -2.17 |
| ESOL_Solubility | decimal | 1:1 | Yes | ESOL solubility (mg/mL) | 1.22 |
| ESOL_Class | string | 1:1 | Yes | Solubility class | Soluble |
| Ali_LogS | decimal | 1:1 | Yes | Ali method log S | -2.36 |
| Ali_Solubility | decimal | 1:1 | Yes | Ali solubility (mg/mL) | 0.79 |
| Ali_Class | string | 1:1 | Yes | Solubility class | Soluble |
| SILICOS_IT_LogS | decimal | 1:1 | Yes | SILICOS-IT log S | -1.89 |
| SILICOS_IT_Solubility | decimal | 1:1 | Yes | SILICOS-IT solubility | 2.31 |
| SILICOS_IT_Class | string | 1:1 | Yes | Solubility class | Soluble |

### Drug-likeness Assessment

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Lipinski_Violations | integer | 1:1 | Yes | Lipinski Rule of 5 violations | 0 |
| Ghose_Violations | integer | 1:1 | Yes | Ghose filter violations | 0 |
| Veber_Violations | integer | 1:1 | Yes | Veber filter violations | 0 |
| Egan_Violations | integer | 1:1 | Yes | Egan filter violations | 0 |
| Muegge_Violations | integer | 1:1 | Yes | Muegge filter violations | 0 |
| Bioavailability_Score | decimal | 1:1 | Yes | Abbott bioavailability score | 0.85 |
| Leadlikeness_Violations | integer | 1:1 | Yes | Lead-likeness violations | 0 |

### Pharmacokinetics Predictions

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| GI_Absorption | string | 1:1 | Yes | GI absorption prediction | High, Low |
| BBB_Permeant | string | 1:1 | Yes | BBB permeation prediction | Yes, No |
| Pgp_Substrate | string | 1:1 | Yes | P-glycoprotein substrate | Yes, No |
| CYP1A2_Inhibitor | string | 1:1 | Yes | CYP1A2 inhibition prediction | Yes, No |
| CYP2C19_Inhibitor | string | 1:1 | Yes | CYP2C19 inhibition prediction | Yes, No |
| CYP2C9_Inhibitor | string | 1:1 | Yes | CYP2C9 inhibition prediction | Yes, No |
| CYP2D6_Inhibitor | string | 1:1 | Yes | CYP2D6 inhibition prediction | Yes, No |
| CYP3A4_Inhibitor | string | 1:1 | Yes | CYP3A4 inhibition prediction | Yes, No |
| Log_Kp | decimal | 1:1 | Yes | Skin permeation (cm/s) | -6.26 |

### Medicinal Chemistry Assessment

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| PAINS_Alerts | integer | 1:1 | Yes | PAINS pattern matches | 0 |
| Brenk_Alerts | integer | 1:1 | Yes | Brenk structural alerts | 0 |
| Synthetic_Accessibility | decimal | 1:1 | Yes | SA score (1-10, lower=easier) | 1.34 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| InChIKey | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | Structure hash |
| SMILES | varies | CC(=O)OC1=CC=CC=C1C(=O)O | Line notation |
| Molecular Formula | varies | C9H8O4 | Atom composition |

---

## Enumerations

### Lipinski Rule of Five Criteria

| Criterion | Threshold | Violation if |
|-----------|-----------|--------------|
| Molecular Weight | <= 500 Da | MW > 500 |
| LogP | <= 5 | LogP > 5 |
| H-bond Donors | <= 5 | HBD > 5 |
| H-bond Acceptors | <= 10 | HBA > 10 |

### Veber Rules

| Criterion | Threshold | Violation if |
|-----------|-----------|--------------|
| TPSA | <= 140 A^2 | TPSA > 140 |
| Rotatable Bonds | <= 10 | Rotatable > 10 |

### Egan Rules

| Criterion | Threshold | Violation if |
|-----------|-----------|--------------|
| WLOGP | <= 5.88 | WLOGP > 5.88 |
| TPSA | <= 131.6 A^2 | TPSA > 131.6 |

### Ghose Rules

| Criterion | Range | Violation if |
|-----------|-------|--------------|
| MW | 160-480 Da | Outside range |
| WLOGP | -0.4 to 5.6 | Outside range |
| MR | 40-130 | Outside range |
| Heavy Atoms | 20-70 | Outside range |

### Muegge Rules

| Criterion | Range | Violation if |
|-----------|-------|--------------|
| MW | 200-600 Da | Outside range |
| XLOGP3 | -2 to 5 | Outside range |
| TPSA | <= 150 A^2 | > 150 |
| Rings | <= 7 | > 7 |
| HBD | <= 5 | > 5 |
| HBA | <= 10 | > 10 |
| Rotatable Bonds | <= 15 | > 15 |

### Solubility Classes

| Class | Log S Range | Interpretation |
|-------|-------------|----------------|
| Insoluble | < -10 | Very poor solubility |
| Poorly soluble | -10 to -6 | Poor solubility |
| Moderately soluble | -6 to -4 | Moderate solubility |
| Soluble | -4 to -2 | Good solubility |
| Very soluble | -2 to 0 | Very good solubility |
| Highly soluble | > 0 | Excellent solubility |

### BOILED-Egg Model Regions

| Region | TPSA Range | WLOGP Range | Interpretation |
|--------|------------|-------------|----------------|
| White (Egg White) | <= 142 | -2.3 to 6.8 | High GI absorption |
| Yellow (Yolk) | <= 79 | 0.4 to 6.0 | High BBB permeation |
| Gray (Outside) | >142 | Out of range | Poor absorption |

### GI/BBB Predictions

| Value | Description |
|-------|-------------|
| High | High absorption/permeation predicted |
| Low | Low absorption/permeation predicted |
| Yes | Substrate/inhibitor predicted |
| No | Not substrate/inhibitor predicted |

---

## Entity Relationships

### Input to Properties
- **Cardinality:** 1:1
- **Description:** Each input molecule produces one set of predictions
- **Key Fields:** Canonical_SMILES, InChIKey

### Properties to Drug-likeness
- **Cardinality:** 1:1
- **Description:** Physicochemical properties determine rule violations
- **Key Fields:** MW, HBD, HBA, LogP, TPSA

### Properties to Pharmacokinetics
- **Cardinality:** 1:1
- **Description:** Properties determine PK predictions
- **Key Fields:** TPSA, WLOGP, structural features

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ADME | Absorption, Distribution, Metabolism, Excretion | Core prediction areas |
| TPSA | Topological Polar Surface Area | Sum of polar atom surfaces |
| HBD | Hydrogen Bond Donor | Atoms donating H-bonds |
| HBA | Hydrogen Bond Acceptor | Atoms accepting H-bonds |
| LogP | Partition Coefficient | Octanol/water ratio |
| LogS | Solubility | Log of mol/L |
| BBB | Blood-Brain Barrier | CNS barrier |
| GI | Gastrointestinal | Absorption site |
| Pgp | P-glycoprotein | Efflux transporter |
| CYP | Cytochrome P450 | Metabolizing enzymes |
| PAINS | Pan-Assay Interference Compounds | False positive patterns |
| SA | Synthetic Accessibility | Synthesis difficulty |
| MR | Molar Refractivity | Volume/polarizability |
| MW | Molecular Weight | Mass in Daltons |
| BOILED-Egg | Brain Or Intestinal EstimateD permeation | Visualization model |
| Csp3 | sp3-hybridized Carbon | Saturation indicator |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubChem | CID (via InChIKey) | Structure lookup |
| ChEMBL | ChEMBL ID (via InChIKey) | Bioactivity data |
| DrugBank | DrugBank ID | Drug properties |
| SuperCYP | Internal | CYP interaction data |

---

## Data Quality Notes

1. **On-Demand Prediction:** All properties computed on-the-fly
2. **Multiple Methods:** 5 LogP methods, 3 solubility methods
3. **Drug-likeness Rules:** 5 established rule sets
4. **CYP Coverage:** Predictions for 5 major isoforms
5. **Batch Processing:** Up to 100 compounds per query
6. **PAINS Patterns:** 480+ structural alert patterns
7. **Open Access:** Free web tool for academic use
8. **Visual Output:** BOILED-Egg plot for absorption/permeation

