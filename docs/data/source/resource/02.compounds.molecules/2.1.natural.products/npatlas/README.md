---
id: npatlas
title: "NPAtlas - Natural Products Atlas"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: natural.products
tags:
  - natural-products
  - microbial
  - marine
  - chemical-structures
  - drug-discovery
---

# NPAtlas - Natural Products Atlas

## Overview

The Natural Products Atlas (NPAtlas) is a curated database of microbial natural products, with a particular emphasis on compounds from bacteria, fungi, and marine microorganisms. The database provides high-quality structural data with comprehensive literature references, organism sources, and compound classification.

NPAtlas distinguishes itself through rigorous manual curation and a focus on microbial secondary metabolites, which represent a major source of bioactive compounds for drug discovery. Each entry includes validated chemical structures, original isolation references, and taxonomic information about the producing organism.

The database is particularly valuable for natural product dereplication, biosynthetic gene cluster analysis, and exploring the chemical diversity of microbial metabolomes.

## Key Statistics

| Metric | Value |
|--------|-------|
| Natural Products | 33,000+ |
| Bacterial Compounds | 15,000+ |
| Fungal Compounds | 18,000+ |
| Marine-Derived | 5,000+ |
| Literature References | 25,000+ |
| Last Update | 2024 |

## Primary Use Cases

1. **Microbial Natural Product Research** - Explore bacterial and fungal metabolites
2. **Dereplication** - Identify known compounds from mass spectrometry data
3. **Biosynthetic Gene Cluster Analysis** - Connect structures to BGC predictions
4. **Marine Natural Products** - Access marine-derived compound data
5. **Taxonomic Profiling** - Study chemical production across microbial taxa

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| NPAtlas ID | NPA + digits | NPA012345 |
| InChI Key | 27 characters | Standard format |
| Canonical SMILES | Variable | Structure string |
| DOI | 10.xxxx/... | Literature reference |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.npatlas.org | Search and browse |
| API | https://www.npatlas.org/api/v1 | REST API |
| Bulk Download | https://www.npatlas.org/download | JSON, SDF formats |

## Organism Categories

| Category | Description |
|----------|-------------|
| Bacteria | Actinomycetes, cyanobacteria, others |
| Fungi | Ascomycetes, basidiomycetes |
| Marine | Marine-derived microorganisms |

## Limitations

- Focus on microbial sources; plant compounds less represented
- Bioactivity data not included for all compounds
- Some historical compounds may have revised structures
- Manual curation means update lag for newly discovered compounds

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [COCONUT](../coconut/_index.md) - All natural products
- [LOTUS](../lotus/_index.md) - Structure-organism pairs
- [MIBiG](../../_index.md) - Biosynthetic gene clusters

## References

1. van Santen JA, et al. (2022) "The Natural Products Atlas 2.0: a database of microbially-derived natural products." Nucleic Acids Res. 50(D1):D1317-D1323.
