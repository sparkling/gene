---
id: supercyp
title: "SuperCYP - Cytochrome P450 Database"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: drug.metabolism.pharmacokinetics
tags:
  - cytochrome-p450
  - drug-metabolism
  - cyp-enzymes
  - drug-interactions
  - pharmacokinetics
---

# SuperCYP - Cytochrome P450 Database

## Overview

SuperCYP is a comprehensive database dedicated to cytochrome P450 (CYP) enzyme interactions with drugs and other xenobiotics. CYP enzymes are responsible for metabolizing approximately 75% of all drugs, making this database essential for understanding drug metabolism, predicting drug-drug interactions, and designing safer drugs.

The database provides detailed information about which CYP isoforms metabolize specific drugs, which compounds inhibit or induce CYP enzymes, and the metabolic pathways involved. This information is crucial for identifying potential drug-drug interactions that can lead to adverse events or therapeutic failure.

SuperCYP supports drug development by providing early assessment of metabolic liabilities and interaction potential, enabling medicinal chemists to design compounds with improved metabolic profiles.

## Key Statistics

| Metric | Value |
|--------|-------|
| Compounds | 6,000+ |
| CYP-Compound Relations | 50,000+ |
| CYP Isoforms | 17 |
| Substrates | 3,500+ |
| Inhibitors | 4,500+ |
| Inducers | 1,200+ |
| Literature References | 2,500+ |

## Primary Use Cases

1. **Drug Interaction Prediction** - Identify potential CYP-mediated interactions
2. **Metabolic Profiling** - Determine which CYPs metabolize compounds
3. **Lead Optimization** - Design compounds with favorable CYP profiles
4. **Clinical Risk Assessment** - Evaluate interaction potential of drug combinations
5. **Regulatory Submission** - Support CYP interaction data for DDI guidance

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Drug Name | Text | Ketoconazole |
| CAS Number | Variable | 65277-42-1 |
| CYP Isoform | CYP + code | CYP3A4 |
| DrugBank ID | DB + 5 digits | Cross-referenced |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://bioinformatics.charite.de/supercyp | Search and browse |
| Download | Available | Contact for access |

## Major CYP Isoforms

| CYP | Drug Metabolism Role |
|-----|---------------------|
| CYP3A4/5 | ~50% of drugs |
| CYP2D6 | ~25% of drugs |
| CYP2C9 | ~15% of drugs |
| CYP2C19 | ~10% of drugs |
| CYP1A2 | Various drugs |
| CYP2B6 | Various drugs |
| CYP2E1 | Alcohol, solvents |

## Interaction Types

| Type | Description |
|------|-------------|
| Substrate | Drug metabolized by CYP |
| Inhibitor | Drug reduces CYP activity |
| Inducer | Drug increases CYP expression |

## Limitations

- Commercial use requires contacting database maintainers
- In vitro CYP data may not predict in vivo interactions
- Some older literature data uses superseded nomenclature
- Genetic polymorphism effects not comprehensively covered

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [SwissADME](../swissadme/README.md) - ADMET prediction
- [DrugBank](../../2.2.pharmaceuticals/drugbank/README.md) - Drug information
- [ChEMBL](../../2.2.pharmaceuticals/chembl/README.md) - Bioactivity

## References

1. Preissner S, et al. (2010) "SuperCYP: a comprehensive database on Cytochrome P450 enzymes including a tool for analysis of CYP-drug interactions." Nucleic Acids Res. 38(Database issue):D237-D243.
