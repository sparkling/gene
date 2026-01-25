---
id: swiss.model
title: "SWISS-MODEL Repository"
type: data-source
category: proteins
subcategory: protein.structures
parent: ../README.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [modeling, homology, templates, 3d, prediction]
---

# SWISS-MODEL Repository

**Category:** [Proteins](../../../README.md) > [Protein Structures](../README.md)

## Overview

SWISS-MODEL is a fully automated protein structure homology-modeling server and repository. It provides 3D protein models built using comparative modeling techniques, where structures are predicted based on evolutionary-related proteins with known experimental structures.

The SWISS-MODEL Repository contains pre-computed models for UniProt sequences, annotated with quality estimates (QMEAN scores). The associated modeling server allows users to build custom models for sequences of interest.

SWISS-MODEL complements experimental structures (PDB) and AI predictions (AlphaFold) by providing template-based models with explicit evolutionary context.

## Key Statistics

| Metric | Value |
|--------|-------|
| Pre-computed Models | 5M+ |
| UniProt Coverage | Extensive |
| Quality Scores | QMEAN, QMEANDisCo |
| Template Database | PDB-derived |
| Last Update | Continuous |

## Primary Use Cases

1. Homology model building
2. Template-based structure prediction
3. Model quality assessment
4. Comparative structural analysis
5. Structure-function relationship studies

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Model ID | Internal | 5a1a8b2c3d |
| UniProt Accession | Standard | P04637 |
| Template PDB | `[0-9][A-Z0-9]{3}` | 1TUP |
| QMEAN Score | -4 to 0 (0 = best) | -1.23 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://swissmodel.expasy.org | Interactive modeling |
| Repository | https://swissmodel.expasy.org/repository | Pre-computed models |
| REST API | https://swissmodel.expasy.org/docs/api | Programmatic access |
| Workspace | Project-based | User projects |

## License

| Aspect | Value |
|--------|-------|
| License | Free for academic use |
| Commercial Use | Requires license |
| Attribution | Required |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [PDB](../pdb/README.md) - Experimental structures (templates)
- [AlphaFold DB](../alphafold.db/README.md) - AI predictions
- [UniProt](../../7.1.protein.sequences.annotations/uniprot/README.md) - Sequences
