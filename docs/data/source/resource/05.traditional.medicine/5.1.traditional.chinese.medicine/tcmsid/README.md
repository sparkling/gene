---
id: tcmsid
title: "TCMSID - TCM Systematic Pharmacology Database"
type: source
parent: ../README.md
tier: 2
status: active
category: traditional.medicine
subcategory: traditional.chinese.medicine
tags:
  - tcm
  - traditional-chinese-medicine
  - adme
  - pharmacokinetics
  - systems-pharmacology
---

# TCMSID - TCM Systematic Pharmacology Database

**Category:** [Traditional Medicine](../../_index.md) > [Traditional Chinese Medicine](../_index.md)

## Overview

TCMSID (Traditional Chinese Medicine Systematic pharmacology database and analysis platform) provides comprehensive ADME/T (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties for TCM compounds. It focuses on the pharmacokinetic and pharmacodynamic properties that determine whether natural compounds can become effective drugs.

The database emphasizes systems pharmacology approaches, integrating compound properties with target networks and pathway analysis. TCMSID uses validated ADME prediction models to assess oral bioavailability, drug-likeness, and Caco-2 permeability for TCM ingredients.

A key feature is the drug-likeness scoring that helps prioritize compounds for experimental validation and drug development efforts.

## Key Statistics

| Metric | Value |
|--------|-------|
| TCM Herbs | 499 |
| Compounds | 29,384 |
| Target Proteins | 3,311 |
| Diseases | 837 |
| Compound-Target Pairs | 98,215 |
| Oral Bioavailability Data | All compounds |
| Drug-likeness Scores | All compounds |

## Primary Use Cases

1. ADME property prediction for TCM compounds
2. Filtering compounds by oral bioavailability
3. Drug-likeness assessment of natural products
4. Systems pharmacology network analysis
5. Prioritizing compounds for experimental validation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Compound ID | MOL prefix | MOL000001 |
| Herb ID | TCMSID internal | TCMSID001 |
| Target ID | UniProt / Gene Symbol | P12345 / TP53 |
| PubChem CID | Integer | 5280343 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://lsp.nwu.edu.cn/tcmsid.php | Interactive search |
| Download | Bulk data available | CSV/Excel files |
| Analysis Tools | Integrated platform | Network visualization |

## Data Model

```
Herbs (499)
    |
    v
Compounds (29,384) ---> ADME Properties
    |                   - Oral Bioavailability (OB)
    |                   - Drug-Likeness (DL)
    |                   - Caco-2 Permeability
    |                   - BBB Permeability
    v
Targets (3,311) -----> Diseases (837)
    |
    v
Pathway Networks
```

## ADME Properties

| Property | Description | Threshold |
|----------|-------------|-----------|
| OB (Oral Bioavailability) | Fraction reaching systemic circulation | >= 30% |
| DL (Drug-Likeness) | Similarity to known drugs | >= 0.18 |
| Caco-2 | Intestinal permeability | >= -0.4 |
| BBB | Blood-brain barrier penetration | >= -0.3 |
| HL (Half-Life) | Pharmacokinetic half-life | Long/Short |
| Lipinski | Rule of five compliance | 0-4 violations |

## Filtering Workflow

| Step | Filter | Purpose |
|------|--------|---------|
| 1 | OB >= 30% | Adequate absorption |
| 2 | DL >= 0.18 | Drug-like properties |
| 3 | Caco-2 >= -0.4 | Intestinal permeability |
| 4 | Target validation | Confirmed interactions |

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

## Limitations

- ADME predictions are computational models
- Not all herbs have comprehensive compound data
- OB thresholds may exclude valid candidates

## See Also

- [BATMAN-TCM](../batman.tcm/_index.md)
- [TCMBank](../tcmbank/_index.md)
- [HERB Database](../herb/_index.md)
