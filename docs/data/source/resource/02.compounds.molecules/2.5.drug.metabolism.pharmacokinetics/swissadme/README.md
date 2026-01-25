---
id: swissadme
title: "SwissADME - ADME Property Prediction"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: drug.metabolism.pharmacokinetics
tags:
  - adme
  - pharmacokinetics
  - drug-likeness
  - physicochemical
  - prediction
---

# SwissADME - ADME Property Prediction

## Overview

SwissADME is a free web tool for computing physicochemical descriptors and predicting ADME (Absorption, Distribution, Metabolism, Excretion) parameters, pharmacokinetics, drug-likeness, and medicinal chemistry friendliness of small molecules. Developed by the Swiss Institute of Bioinformatics, it provides rapid assessment of compound properties essential for drug discovery.

The tool offers multiple drug-likeness rules (Lipinski, Veber, Egan, Ghose, Muegge), bioavailability scoring, P-glycoprotein substrate prediction, CYP inhibition prediction, and synthetic accessibility estimation. The intuitive BOILED-Egg model provides visual representation of gastrointestinal absorption and brain penetration potential.

SwissADME enables medicinal chemists and drug discovery scientists to quickly assess the drug development potential of compounds and prioritize candidates for synthesis and testing.

## Key Statistics

| Metric | Value |
|--------|-------|
| Physicochemical Descriptors | 20+ |
| Drug-likeness Filters | 5 rules |
| CYP Inhibition Predictions | 5 isoforms |
| Compounds per Query | Up to 100 |
| Response Time | Seconds |

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| SMILES | String | CC(=O)OC1=CC=CC=C1C(=O)O |
| InChI | InChI=... | InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12) |
| InChIKey | 27-char | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| Canonical SMILES | String | Standardized SMILES notation |

## Primary Use Cases

1. **Lead Prioritization** - Rank compounds by drug-likeness
2. **ADME Assessment** - Predict absorption and metabolism properties
3. **CNS Penetration** - Assess blood-brain barrier permeability
4. **Virtual Screening Filter** - Eliminate poor drug candidates early
5. **Medicinal Chemistry** - Guide structural modifications

## Predicted Properties

| Category | Properties |
|----------|------------|
| Physicochemical | MW, LogP, TPSA, HBD, HBA, rotatable bonds |
| Drug-likeness | Lipinski, Veber, Egan, Ghose, Muegge |
| Pharmacokinetics | GI absorption, BBB permeation, Pgp substrate |
| Metabolism | CYP1A2, 2C19, 2C9, 2D6, 3A4 inhibition |
| Medicinal Chemistry | PAINS, Brenk alerts, synthetic accessibility |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://www.swissadme.ch | Interactive tool |
| Batch Input | SMILES list | Up to 100 compounds |
| Results Export | CSV | Download predictions |

## Drug-Likeness Rules

| Rule | Criteria |
|------|----------|
| Lipinski | MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10 |
| Veber | TPSA <= 140, Rotatable bonds <= 10 |
| Egan | LogP <= 5.88, TPSA <= 131.6 |
| Ghose | 160 <= MW <= 480, -0.4 <= LogP <= 5.6 |
| Muegge | 200 <= MW <= 600, -2 <= LogP <= 5 |

## Limitations

- Predictions are computational estimates, not experimental values
- Batch processing limited to 100 compounds
- No bulk download of training data or models
- Predictions may be less accurate for chemotypes outside training set

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [SuperCYP](../supercyp/README.md) - CYP interactions
- [ChEMBL](../../2.2.pharmaceuticals/chembl/README.md) - Bioactivity data
- [DrugBank](../../2.2.pharmaceuticals/drugbank/README.md) - Drug properties

## References

1. Daina A, Michielin O, Zoete V. (2017) "SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules." Sci Rep. 7:42717. DOI: 10.1038/srep42717
