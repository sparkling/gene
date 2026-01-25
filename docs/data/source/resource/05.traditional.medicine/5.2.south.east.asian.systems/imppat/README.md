---
id: imppat
title: "IMPPAT 2.0 - Indian Medicinal Plants, Phytochemistry And Therapeutics"
type: source
parent: ../README.md
tier: 1
status: active
category: traditional.medicine
subcategory: south.east.asian.systems
tags:
  - ayurveda
  - indian-medicine
  - phytochemistry
  - drug-likeness
  - admet
  - natural-products
---

# IMPPAT 2.0 - Indian Medicinal Plants, Phytochemistry And Therapeutics

## Overview

IMPPAT 2.0 (Indian Medicinal Plants, Phytochemistry And Therapeutics) is the largest digital database on phytochemicals of Indian medicinal plants. It integrates data from Ayurveda, Siddha, Unani, and Homeopathy traditions, providing comprehensive information on plants, their chemical constituents, and therapeutic uses.

A distinguishing feature is the granular plant-part-level associations, specifying which compounds are found in roots, leaves, bark, seeds, or other tissues. IMPPAT 2.0 provides extensive ADMET predictions from SwissADME, 1,875 molecular descriptors per compound, and drug-likeness assessments across six standard filters.

The database includes predicted human target proteins from STITCH with high-confidence interactions (score >= 700), enabling network pharmacology analysis of Indian traditional medicine.

## Key Statistics

| Metric | Value |
|--------|-------|
| Indian Medicinal Plants | 4,010 |
| Phytochemicals | 17,967 |
| Therapeutic Uses | 1,095 |
| Plant-Phytochemical Associations | 189,386 |
| Plant-Therapeutic Associations | 89,733 |
| Predicted Human Targets | 5,042 |
| Compound-Target Interactions | 27,365 |
| Molecular Descriptors/Compound | 1,875 |
| Drug-like Compounds | 1,335 |

## Primary Use Cases

1. Discovering active compounds from Ayurvedic plants
2. Drug-likeness screening of natural products
3. ADMET property assessment for herbal compounds
4. Target prediction for Indian medicinal plant constituents
5. Ethnobotanical data mining by therapeutic use

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Plant ID | IMPPAT_PLANT_XXXXX | IMPPAT_PLANT_03985 |
| Compound ID | IMPPAT_CHEM_XXXXX | IMPPAT_CHEM_15234 |
| InChIKey | 27-character string | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| HGNC ID | HGNC:XXXXX | HGNC:11998 |
| PubChem CID | Via UniChem | 5280343 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://cb.imsc.res.in/imppat/ | Interactive browsing |
| TSV Export | Via search results | Bulk associations |
| Structure Download | SDF, MOL, MOL2, PDB, PDBQT | Per compound |
| GitHub | https://github.com/asamallab/IMPPAT2 | Analysis scripts |

## Data Model

```
Plants (4,010) --- Traditional Systems (Ayurveda, Siddha, Unani, Homeopathy)
    |
    v
Plant Parts (root, leaf, bark, seed, etc.)
    |
    v
Phytochemicals (17,967) ---> Drug-Likeness (6 filters)
    |                        ADMET Properties
    |                        1,875 Descriptors
    v
Target Proteins (5,042) ---> STITCH score >= 700
```

## Compound Data Tabs (6 per compound)

| Tab | Content |
|-----|---------|
| 1. Summary | IDs, structure, SMILES, InChI, source plants |
| 2. Physicochemical | MW, logP, TPSA, HBD, HBA, rotatable bonds |
| 3. Drug-Likeness | Lipinski, Ghose, Veber, Egan, Pfizer, GSK, QED |
| 4. ADMET | GI absorption, BBB, CYP inhibition, toxicity |
| 5. Descriptors | 1,875 2D/3D chemical descriptors |
| 6. Targets | Predicted proteins (STITCH >= 700) |

## Drug-Likeness Filters

| Filter | Criteria | Drug-like Pass |
|--------|----------|----------------|
| Lipinski RO5 | MW<=500, logP<=5, HBD<=5, HBA<=10 | 1,335 compounds |
| Ghose | MW:160-480, logP:-0.4-5.6 | - |
| Veber | RotBonds<=10, TPSA<=140 | - |
| QED Score | 0-1 scale | >0.5 favorable |

## Limitations

- No public REST API (manual export required)
- 493 compounds lack ADMET due to SMILES length
- Target predictions only (no experimental validation)
- Confidence threshold excludes lower-confidence interactions

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
