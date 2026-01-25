# PDB - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pdb |
| **Name** | Protein Data Bank |
| **Total Fields** | 200+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| entry_id | string | 1:1 | Yes | PDB ID (4 chars) | `1TUP` |
| title | string | 1:1 | Yes | Structure title | `Crystal structure of...` |
| deposition_date | date | 1:1 | No | Deposition date | `2024-01-15` |
| release_date | date | 1:1 | No | Release date | `2024-03-01` |

### Experimental Details

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| experimental_method | string | 1:1 | Yes | Method | `X-RAY DIFFRACTION` |
| resolution | number | 1:1 | No | Resolution (Angstroms) | `1.7` |
| r_work | number | 1:1 | No | R-work value | `0.194` |
| r_free | number | 1:1 | No | R-free value | `0.256` |

### Coordinate Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| atom_id | integer | 1:N | Yes | Atom serial number | `1` |
| type_symbol | string | 1:1 | Yes | Element symbol | `C`, `N`, `O` |
| Cartn_x | float | 1:1 | Yes | X coordinate | `10.123` |
| B_iso | float | 1:1 | Yes | B-factor | `15.23` |
| occupancy | float | 1:1 | Yes | Occupancy (0-1) | `1.00` |

---

## Enumerations

### Experimental Methods

| Value | Description |
|-------|-------------|
| X-RAY DIFFRACTION | X-ray crystallography |
| ELECTRON MICROSCOPY | Cryo-EM |
| SOLUTION NMR | Solution NMR |
| SOLID-STATE NMR | Solid-state NMR |
| NEUTRON DIFFRACTION | Neutron diffraction |

### Entity Types

| Value | Description |
|-------|-------------|
| polymer | Proteins, nucleic acids |
| non-polymer | Ligands, ions |
| water | Water molecules |

### Quality Metrics

| Metric | Good | Acceptable |
|--------|------|------------|
| Resolution | < 2.0 A | < 3.0 A |
| R-free | < 0.25 | < 0.30 |
| Clashscore | < 5 | < 20 |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PDB | Protein Data Bank | Archive |
| wwPDB | Worldwide PDB | Consortium |
| RCSB | Research Collaboratory for Structural Bioinformatics | US partner |
| mmCIF | macromolecular CIF | Primary format |
| EM | Electron Microscopy | Method |
| NMR | Nuclear Magnetic Resonance | Method |
| FSC | Fourier Shell Correlation | Resolution method |

---

## Data Quality Notes

1. **mmCIF primary format**: PDB format limited, prefer mmCIF
2. **Validation reports**: Available for all entries
3. **Resolution varies**: X-ray typically < 3A, EM improving
4. **Weekly updates**: New structures released Wednesday

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
