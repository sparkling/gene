# SWISS-MODEL - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | swiss-model |
| **Name** | SWISS-MODEL Repository |
| **Total Fields** | 50+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| uniprot_ac | string | 1:1 | Yes | UniProt accession | `P04637` |
| sequence_id | string | 1:1 | No | Internal sequence ID | `SM00001` |
| model_id | string | 1:N | Yes | Unique model identifier | `5a1b2c3d4e` |
| sequence_length | integer | 1:1 | No | Full protein length | `393` |

### Model Coverage

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| from | integer | 1:1 | Yes | Start residue | `94` |
| to | integer | 1:1 | Yes | End residue | `289` |
| coverage | float | 1:1 | Yes | Sequence coverage (0-1) | `0.50` |

### Template Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| template | string | 1:1 | Yes | Template PDB ID + chain | `1tup.1.A` |
| identity | float | 1:1 | Yes | Sequence identity (0-1) | `0.95` |
| similarity | float | 1:1 | No | Sequence similarity (0-1) | `0.98` |
| oligo_state | string | 1:1 | No | Oligomeric state | `tetramer` |

### Quality Metrics

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| qmean | float | 1:1 | Yes | QMEAN Z-score | `-0.85` |
| qmean_disco | float | 1:1 | No | QMEANDisCo score (0-1) | `0.78` |
| qmean_torsion | float | 1:1 | No | Torsion angle score | `-0.5` |
| qmean_solvation | float | 1:1 | No | Solvation score | `-0.3` |
| qmean_interaction | float | 1:1 | No | Atom interaction score | `-0.4` |
| qmean_local | float | 1:N | No | Per-residue QMEAN (0-1) | `0.72` |

---

## Enumerations

### Oligomeric States

| Value | Description |
|-------|-------------|
| monomer | Single chain |
| homodimer | Two identical chains |
| homotrimer | Three identical chains |
| homotetramer | Four identical chains |
| heterodimer | Two different chains |
| heteromultimer | Multiple different chains |

### Project Status

| Value | Description |
|-------|-------------|
| PENDING | In queue |
| RUNNING | Processing |
| COMPLETED | Models available |
| FAILED | Error occurred |

### Template Methods

| Value | Description |
|-------|-------------|
| X-ray | X-ray crystallography |
| EM | Electron microscopy |
| NMR | Nuclear magnetic resonance |

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

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| QMEAN | Qualitative Model Energy ANalysis | Quality score |
| QMEANDisCo | QMEAN with Distance Constraints | Enhanced metric |
| SMTL | SWISS-MODEL Template Library | Template database |
| PDB | Protein Data Bank | Template source |
| GMQE | Global Model Quality Estimate | Alternative metric |

---

## Data Quality Notes

1. **Template-based**: Models depend on template quality and identity
2. **Minimum identity**: 30% sequence identity recommended
3. **Coverage varies**: Some proteins have multiple partial models
4. **Weekly updates**: New templates incorporated from PDB

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
