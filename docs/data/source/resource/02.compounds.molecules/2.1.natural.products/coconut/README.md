---
id: coconut
title: "COCONUT - COlleCtion of Open Natural prodUcTs"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: natural.products
tags:
  - natural-products
  - compounds
  - chemistry
  - drug-discovery
  - phytochemicals
---

# COCONUT - COlleCtion of Open Natural prodUcTs

## Overview

COCONUT (COlleCtion of Open Natural prodUcTs) is one of the largest open databases of natural products, aggregating compound data from over 50 source databases. It provides comprehensive structural, physicochemical, and biological information about naturally-occurring compounds from plants, microorganisms, marine organisms, and other biological sources.

The database serves as a key resource for natural product drug discovery, enabling researchers to identify bioactive scaffolds, explore chemical space diversity, and connect compounds to their biological source organisms. COCONUT is particularly valuable for identifying drug-like natural products through integrated Lipinski Rule of 5 analysis and QED drug-likeness scoring.

COCONUT integrates with major chemistry databases including PubChem, ChEMBL, ZINC, and ChEBI, providing extensive cross-references that support multi-database research workflows.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Natural Products | ~450,000 |
| Unique Organisms | ~35,000 |
| Data Sources Aggregated | 52 |
| Average Molecular Weight | 312 Da |
| Drug-like Compounds (Lipinski 0) | ~38% |
| Last Update | 2025 |

## Primary Use Cases

1. **Natural Product Drug Discovery** - Identify bioactive scaffolds and lead compounds from natural sources
2. **Structure-Activity Relationship Studies** - Explore chemical diversity and scaffold analysis
3. **Compound Dereplication** - Check novelty of isolated natural products against known compounds
4. **Organism-Compound Mapping** - Link natural products to their biological source species
5. **Drug-Likeness Screening** - Filter compounds by physicochemical properties

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| COCONUT ID | CNP + 7 digits | CNP0123456 |
| InChI Key | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| Canonical SMILES | Variable | CC(C)CCCC(C)C |
| CAS Number | Variable | 50-99-7 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://coconut.naturalproducts.net | Interactive search and browsing |
| REST API | https://coconut.naturalproducts.net/api | Programmatic access |
| Bulk Download | https://coconut.naturalproducts.net/download | SDF, CSV, JSON formats |
| GitHub | https://github.com/Steinbeck-Lab/coconut | Source code and releases |

## Data Formats

| Format | Size (compressed) | Description |
|--------|-------------------|-------------|
| SDF | ~500 MB | Complete structures with properties |
| CSV | ~300 MB | Tabular format |
| JSON | ~600 MB | Full metadata |
| SMILES | ~200 MB | Structure strings only |

## Limitations

- Aggregated data quality varies across source databases
- Some entries may lack complete structural validation
- Organism-compound links not always experimentally verified
- Updates may lag behind primary literature

## Related Resources

- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Bioactivity data
- [PubChem](../../2.6.chemical.ontology.classification/pubchem/_index.md) - Chemical information
- [LOTUS](../lotus/_index.md) - Structure-organism pairs
- [ChEBI](../../2.6.chemical.ontology.classification/chebi/_index.md) - Chemical ontology

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## References

1. Sorokina M, Merseburger P, Rajan K, Yirik MA, Steinbeck C. (2021) "COCONUT online: Collection of Open Natural Products database." Journal of Cheminformatics 13, 2. DOI: 10.1186/s13321-020-00478-9
