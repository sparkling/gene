---
id: download-imppat
title: "IMPPAT 2.0 Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# IMPPAT 2.0 Download Instructions

## Quick Start

```bash
# Visit the web interface - no REST API available
# URL: https://cb.imsc.res.in/imppat/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ with RDKit for structure processing (optional)
- Sufficient storage for structure files (~200MB)

## Download Methods

### Primary: Web Interface

Access IMPPAT through the main portal:
- **URL**: https://cb.imsc.res.in/imppat/
- **Data Export**: TSV via search results
- **Structure Download**: SDF, MOL, MOL2, PDB, PDBQT per compound

### GitHub Repository

Analysis scripts (not raw data):
- **URL**: https://github.com/asamallab/IMPPAT2

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| Indian Medicinal Plants | 4,010 | Plant metadata |
| Phytochemicals | 17,967 | Compounds with 1,875 descriptors |
| Therapeutic Uses | 1,095 | Traditional indications |
| Plant-Phytochemical Links | 189,386 | Part-specific associations |
| Plant-Therapeutic Links | 89,733 | Usage associations |
| Predicted Targets | 5,042 | Human proteins (STITCH) |
| Compound-Target Links | 27,365 | High-confidence interactions |

## Data Structure

### Compound Data (6 Tabs per compound)

| Tab | Content |
|-----|---------|
| 1. Summary | IDs, structure, SMILES, InChI, source plants |
| 2. Physicochemical | MW, logP, TPSA, HBD, HBA, rotatable bonds |
| 3. Drug-Likeness | Lipinski, Ghose, Veber, Egan, Pfizer, GSK, QED |
| 4. ADMET | GI absorption, BBB, CYP inhibition, toxicity |
| 5. Descriptors | 1,875 2D/3D chemical descriptors |
| 6. Targets | Predicted proteins (STITCH >= 700) |

### Structure File Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| SDF | .sdf | Standard 2D/3D |
| MOL | .mol | Single molecule |
| MOL2 | .mol2 | With atom types |
| PDB | .pdb | Docking input |
| PDBQT | .pdbqt | AutoDock ready |

## Download Steps

1. Navigate to https://cb.imsc.res.in/imppat/
2. Use search (plant name, compound, therapeutic use)
3. Export search results as TSV
4. Download structure files for compounds of interest
5. Use GitHub scripts for analysis

## Bulk Download Strategy

```python
# Recommended pipeline for comprehensive data collection

# 1. Search by plant family or therapeutic use
# 2. Export association tables (TSV)
# 3. Download SDF files for compound library
# 4. Map to PubChem via InChIKey
# 5. Link targets via HGNC to UniProt

import pandas as pd
from rdkit import Chem

# Load exported associations
associations = pd.read_csv('imppat_plant_phytochem.tsv', sep='\t')
print(f"Plant-compound associations: {len(associations)}")

# Load structures
suppl = Chem.SDMolSupplier('imppat_compounds.sdf')
valid_mols = [m for m in suppl if m is not None]
print(f"Valid structures: {len(valid_mols)}")
```

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Plants | 4,010 |
| Phytochemicals | 17,967 |
| Therapeutic Uses | 1,095 |
| Plant-Compound Links | 189,386 |
| Targets | 5,042 |

## Drug-Likeness Subset

1,335 compounds pass all 6 drug-likeness filters:
- Lipinski Rule of 5
- Ghose filter
- Veber filter
- Egan filter
- Pfizer 3/75 rule
- GSK 4/400 rule

## Update Schedule

| Aspect | Value |
|--------|-------|
| Current Version | IMPPAT 2.0 (June 2022) |
| Update Frequency | Periodic major updates |
| License | CC BY-NC 4.0 |

---

## Dataset Versions

### Current Release: IMPPAT 2.0

| Property | Value |
|----------|-------|
| Version | 2.0 |
| Release Date | 2022-06-01 |
| Total Size | ~200 MB |
| Descriptors | 1,875 per compound |

### Version Contents

| Component | Records | Description |
|-----------|---------|-------------|
| Indian Medicinal Plants | 4,010 | Plant metadata |
| Phytochemicals | 17,967 | Compounds with 1,875 descriptors |
| Therapeutic Uses | 1,095 | Traditional indications |
| Plant-Phytochemical Links | 189,386 | Part-specific associations |
| Plant-Therapeutic Links | 89,733 | Usage associations |
| Predicted Targets | 5,042 | Human proteins (STITCH) |

### Previous Versions

| Version | Release | Status |
|---------|---------|--------|
| IMPPAT 1.0 | 2018-01-01 | Archived |

### Drug-Likeness Subset

| Criteria | Compounds Passing |
|----------|-------------------|
| All 6 filters (Lipinski, Ghose, Veber, Egan, Pfizer, GSK) | 1,335 |

---

## Notes

- No public REST API (manual export required)
- 493 compounds lack ADMET due to SMILES length
- Target predictions only (STITCH score >= 700)
- Plant-part-level associations unique to IMPPAT
- Comprehensive coverage of Indian traditional medicine
