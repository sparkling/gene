---
id: download-tcmbank
title: "TCMBank Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# TCMBank Download Instructions

## Quick Start

```bash
# Visit the web interface for structure downloads
# URL: http://tcmbank.cn/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ with RDKit for structure processing (optional)
- Sufficient storage for SDF files (~1GB)

## Download Methods

### Primary: Web Interface

Access TCMBank through the main portal:
- **URL**: http://tcmbank.cn/
- **Data Access**: Search and download by compound

### Structure Downloads

TCMBank provides molecular structure files:

| Format | Extension | Description |
|--------|-----------|-------------|
| SDF | .sdf | Standard structure-data format |
| MOL | .mol | Single molecule format |
| SMILES | Text | Linear notation |

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| TCM Herbs | 9,000+ | Herb metadata |
| Compounds | 60,000+ | Chemical structures with properties |
| Prescriptions/Formulas | 75,000+ | Traditional formulations |
| Target Proteins | 15,000+ | Predicted targets |
| Diseases | 8,000+ | Disease associations |
| Compound-Target Pairs | 250,000+ | Interaction predictions |

## Data Structure

### Compound Data Fields
- Compound ID (TCMBank internal)
- Compound name
- SMILES notation
- InChI / InChIKey
- Molecular formula
- Molecular weight (MW)
- LogP (lipophilicity)
- TPSA (polar surface area)
- Hydrogen bond donors/acceptors
- Rotatable bonds
- Drug-likeness scores (Lipinski, Veber, PAINS)
- ADMET predictions
- Source herbs
- Predicted targets

### Herb Data Fields
- Herb ID
- Chinese name
- Pinyin name
- English name
- Latin botanical name
- Compound count

## Download Steps

1. Navigate to http://tcmbank.cn/
2. Use search function (by herb, compound, target, or structure)
3. Browse compound pages
4. Download SDF files for structures
5. Export search results for metadata

## Bulk Download Strategy

For comprehensive compound library:

```python
# Recommended approach for structure collection
# 1. Export compound list via search
# 2. Download SDF files in batches
# 3. Combine into master SDF

from rdkit import Chem

# Load and verify structures
suppl = Chem.SDMolSupplier('tcmbank_compounds.sdf')
valid_mols = [m for m in suppl if m is not None]
print(f"Valid structures: {len(valid_mols)}")
```

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Herbs | 9,000+ |
| Compounds | 60,000+ |
| Formulas | 75,000+ |
| Targets | 15,000+ |

## Drug-Likeness Filters Available

| Filter | Criteria |
|--------|----------|
| Lipinski RO5 | MW <= 500, logP <= 5, HBD <= 5, HBA <= 10 |
| Veber | RotBonds <= 10, TPSA <= 140 |
| PAINS | No pan-assay interference compounds |
| Lead-like | MW 250-350, logP <= 3.5 |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Periodic |
| Notification | Check website |

---

## Dataset Versions

### Current Release: TCMBank

| Property | Value |
|----------|-------|
| Version | 1.0 |
| Release Date | 2022-01-01 |
| Total Size | ~1 GB |
| Focus | Molecular Structures |

### Version Contents

| Component | Records | Description |
|-----------|---------|-------------|
| TCM Herbs | 9,000+ | Herb metadata |
| Compounds | 60,000+ | Chemical structures with properties |
| Prescriptions/Formulas | 75,000+ | Traditional formulations |
| Target Proteins | 15,000+ | Predicted targets |
| Diseases | 8,000+ | Disease associations |
| Compound-Target Pairs | 250,000+ | Interaction predictions |

### Structure Data

| Format | Description |
|--------|-------------|
| SDF | Standard structure-data format |
| MOL | Single molecule format |
| SMILES | Linear notation |

---

## Notes

- Largest TCM compound structure collection
- ADMET predictions are computational
- Structure files ideal for virtual screening
- Chinese literature focus may limit some coverage
- Contact maintainers for commercial use
