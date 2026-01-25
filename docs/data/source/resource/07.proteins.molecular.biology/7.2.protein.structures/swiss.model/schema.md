---
id: schema-swiss-model
title: "SWISS-MODEL Repository Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: final
tags: [schema, database, structures, homology, modeling, templates]
---

# SWISS-MODEL Repository Schema Documentation

**Document ID:** SCHEMA-SWISS-MODEL
**Version:** 2026.01
**Source Version:** SWISS-MODEL Repository 2024

---

## TL;DR

SWISS-MODEL Repository provides pre-computed homology models for UniProt sequences built using template-based modeling. Each model includes quality estimates (QMEAN, QMEANDisCo), template information, and coverage statistics. The modeling server allows custom model building for sequences without pre-computed structures.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Pre-computed models | 5,000,000+ | SWISS-MODEL Stats |
| UniProt coverage | Extensive | SWISS-MODEL Stats |
| Template structures | All PDB | SWISS-MODEL Stats |
| Average sequence identity | >30% | SWISS-MODEL Stats |
| Model formats | PDB, mmCIF | SWISS-MODEL Docs |

---

## Entity Relationship Overview

```
SWISS-MODEL Entry
  ├── Target (UniProt sequence)
  │     ├── Accession
  │     ├── Sequence
  │     └── Coverage range
  ├── Model
  │     ├── Coordinates
  │     ├── QMEAN scores
  │     └── Per-residue quality
  ├── Template
  │     ├── PDB ID + chain
  │     ├── Sequence identity
  │     └── Coverage
  └── Alignment
        ├── Target sequence
        ├── Template sequence
        └── Secondary structure
```

---

## Core Tables/Entities

### Model Entry

**Description:** Pre-computed homology model for a UniProt sequence.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| uniprot_ac | string | Yes | UniProt accession |
| sequence_id | string | Yes | Internal sequence ID |
| model_id | string | Yes | Unique model identifier |
| from | integer | Yes | Start residue in UniProt |
| to | integer | Yes | End residue in UniProt |
| template | string | Yes | Template PDB ID + chain |
| identity | float | Yes | Sequence identity (0-1) |
| similarity | float | No | Sequence similarity (0-1) |
| coverage | float | Yes | Sequence coverage (0-1) |
| qmean | float | Yes | QMEAN Z-score |
| qmean_disco | float | No | QMEANDisCo score (0-1) |
| created_date | date | Yes | Model generation date |
| method | string | Yes | Modeling method |

### Quality Metrics

**Description:** Model quality assessment scores.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| qmean | float | Yes | Global QMEAN Z-score (-4 to 0) |
| qmean_disco | float | No | Distance-constraint QMEANDisCo (0-1) |
| qmean_torsion | float | No | Torsion angle score |
| qmean_solvation | float | No | Solvation score |
| qmean_interaction | float | No | Atom interaction score |
| qmean_cbeta | float | No | Cbeta interaction score |
| qmean_all_atom | float | No | All-atom packing score |

### Template Information

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| template_id | string | Yes | PDB ID (e.g., 1tup.1.A) |
| pdb_id | string | Yes | 4-character PDB ID |
| chain | string | Yes | Chain identifier |
| resolution | float | No | Template resolution (A) |
| method | string | No | Experimental method |
| oligo_state | string | No | Oligomeric state |
| ligands | array | No | Bound ligands |

### Per-Residue Quality

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| residue_number | integer | Yes | Residue position |
| residue_name | string | Yes | 3-letter amino acid |
| qmean_local | float | Yes | Local QMEAN (0-1) |
| secondary_structure | string | No | H, E, C |
| template_residue | string | No | Aligned template residue |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/repository/uniprot/{accession}` | GET | Models for UniProt entry |
| `/repository/uniprot/{accession}.pdb` | GET | Download model (PDB) |
| `/automodel` | POST | Submit modeling job |
| `/project/{id}` | GET | Get project status |
| `/project/{id}/models` | GET | Get project models |

### Repository API Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| provider | Model source | swissmodel |
| template | Filter by template | 1tup |
| sort | Sort order | qmean |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | PDB format |
| Alternative | mmCIF, SMTL (project files) |
| Quality | JSON (API) |
| Encoding | UTF-8 |

---

## Sample Record

### API Response (JSON)

```json
{
  "result": {
    "uniprot_ac": "P04637",
    "sequence_length": 393,
    "models": [
      {
        "model_id": "5a1b2c3d4e",
        "from": 94,
        "to": 289,
        "template": "1tup.1.A",
        "identity": 0.95,
        "coverage": 0.50,
        "qmean": -0.85,
        "qmean_disco": 0.78,
        "oligo_state": "monomer",
        "created_date": "2024-06-15",
        "coordinate_url": "https://swissmodel.expasy.org/repository/uniprot/P04637/1.pdb"
      },
      {
        "model_id": "6b2c3d4e5f",
        "from": 319,
        "to": 360,
        "template": "3sak.1.A",
        "identity": 1.00,
        "coverage": 0.10,
        "qmean": -1.23,
        "qmean_disco": 0.65,
        "oligo_state": "tetramer",
        "created_date": "2024-06-15"
      }
    ],
    "coverage_image": "https://swissmodel.expasy.org/repository/uniprot/P04637.coverage.png"
  }
}
```

### PDB Format Header (Model File)

```
HEADER    SWISS-MODEL REPOSITORY
TITLE     Model of CELLULAR TUMOR ANTIGEN P53
EXPDTA    THEORETICAL MODEL (SWISS-MODEL SERVER)
AUTHOR    SWISS-MODEL SERVER
REMARK   1
REMARK   1 MODEL INFORMATION
REMARK   1  UNIPROT     P04637
REMARK   1  RANGE       94-289
REMARK   1  TEMPLATE    1tup.1.A
REMARK   1  IDENTITY    0.95
REMARK   1  COVERAGE    0.50
REMARK   1  QMEAN       -0.85
REMARK   1
REMARK   3  GLOBAL QMEAN SCORE
REMARK   3   QMEAN                 -0.85
REMARK   3   QMEAN_DISCO            0.78
REMARK   3
REMARK   3  PER-RESIDUE QMEAN SCORES IN B-FACTOR COLUMN
ATOM      1  N   MET A  94      10.123  20.456  30.789  1.00  0.72           N
ATOM      2  CA  MET A  94      11.234  21.567  31.890  1.00  0.72           C
...
```

---

## Quality Score Interpretation

### QMEAN Z-Score

| Score | Quality | Interpretation |
|-------|---------|----------------|
| > -1 | Excellent | High-confidence model |
| -1 to -2 | Good | Reliable for most uses |
| -2 to -3 | Moderate | Use with caution |
| < -3 | Poor | Low reliability |

### QMEANDisCo Score (0-1)

| Score | Quality | Interpretation |
|-------|---------|----------------|
| > 0.7 | High | Very reliable |
| 0.5-0.7 | Medium | Generally acceptable |
| < 0.5 | Low | Use with caution |

### Local Quality (Per-residue)

| Score | Quality | Description |
|-------|---------|-------------|
| > 0.7 | High | Reliable local structure |
| 0.4-0.7 | Medium | Moderate confidence |
| < 0.4 | Low | Potentially unreliable |

---

## Modeling Server Project

### Project Submission

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| target_sequence | string | Yes | FASTA sequence |
| project_name | string | No | User-defined name |
| email | string | Optional | Notification email |

### Project Status

| Status | Description |
|--------|-------------|
| PENDING | In queue |
| RUNNING | Processing |
| COMPLETED | Models available |
| FAILED | Error occurred |

### Project Results

| Field | Type | Description |
|-------|------|-------------|
| models | array | Generated models |
| templates | array | Templates identified |
| alignments | array | Sequence alignments |
| errors | array | Processing errors |

---

## Template Selection Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Sequence identity | High | >30% preferred |
| Coverage | High | Maximize coverage |
| Resolution | Medium | Lower is better |
| QMEAN | Medium | Template quality |
| Oligomeric state | Low | Match expected state |

---

## Glossary

| Term | Definition |
|------|------------|
| Homology model | 3D structure based on related template |
| Template | Experimental structure used for modeling |
| QMEAN | Qualitative Model Energy ANalysis score |
| QMEANDisCo | QMEAN with distance constraints |
| Sequence identity | Percentage identical residues |
| Coverage | Fraction of target modeled |
| Oligomeric state | Biological assembly form |

---

## References

1. Waterhouse A, et al. (2018). SWISS-MODEL: homology modelling of protein structures and complexes. Nucleic Acids Research. https://doi.org/10.1093/nar/gky427
2. Studer G, et al. (2020). QMEANDisCo - distance constraints applied on model quality estimation. Bioinformatics. https://doi.org/10.1093/bioinformatics/btz828
3. SWISS-MODEL Documentation: https://swissmodel.expasy.org/docs/
