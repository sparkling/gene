---
id: pdb
title: "PDB - Protein Data Bank"
type: data-source
category: proteins
subcategory: protein.structures
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [structures, crystallography, cryo-em, nmr, experimental, 3d]
---

# PDB - Protein Data Bank

**Category:** [Proteins](../../../_index.md) > [Protein Structures](../_index.md)

## Overview

The Protein Data Bank (PDB) is the single worldwide archive of experimentally determined three-dimensional structures of proteins, nucleic acids, and complex assemblies. It is managed by the Worldwide Protein Data Bank (wwPDB) partnership of RCSB PDB (US), PDBe (Europe), PDBj (Japan), and BMRB.

PDB contains structures determined by X-ray crystallography, cryo-electron microscopy (cryo-EM), and NMR spectroscopy. Each entry includes atomic coordinates, experimental details, and extensive validation information.

PDB is the authoritative source for experimental protein structures and is fundamental to structural biology, drug discovery, and molecular biophysics.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Structures | 220,000+ |
| Protein Structures | 200,000+ |
| Nucleic Acid Structures | 15,000+ |
| Cryo-EM Structures | 25,000+ |
| Last Update | Weekly |

## Primary Use Cases

1. Structure-based drug design
2. Molecular dynamics simulations
3. Protein engineering
4. Structure validation and quality assessment
5. Educational and visualization resources

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| PDB ID | `[0-9][A-Z0-9]{3}` | 1TUP |
| Extended PDB ID | `pdb_[0-9]{8}` | pdb_00001tup |
| Chain ID | `[A-Za-z0-9]` | A |
| Entity ID | Numeric | 1 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| RCSB PDB | https://www.rcsb.org | US portal |
| PDBe | https://www.ebi.ac.uk/pdbe | European portal |
| PDBj | https://pdbj.org | Japan portal |
| REST API | https://data.rcsb.org | Programmatic access |
| FTP | https://ftp.wwpdb.org | Bulk downloads |

## License

| Aspect | Value |
|--------|-------|
| License | CC0 (Public Domain) |
| Commercial Use | Yes |
| Attribution | Recommended (cite original publication) |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [AlphaFold DB](../alphafold.db/README.md) - Predicted structures
- [SWISS-MODEL](../swiss.model/README.md) - Homology models
