---
id: npass
title: "NPASS - Natural Product Activity and Species Source Database"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: natural.products
tags:
  - natural-products
  - bioactivity
  - species-source
  - drug-discovery
---

# NPASS - Natural Product Activity and Species Source Database

## Overview

NPASS (Natural Product Activity and Species Source Database) is a comprehensive resource that integrates natural product structural data with quantitative bioactivity information and species source data. Unlike databases that focus solely on structure or occurrence, NPASS emphasizes the connection between natural products, their biological activities, and the organisms that produce them.

The database provides curated bioactivity data including IC50, EC50, and Ki values for natural products tested against various biological targets. This quantitative activity information is particularly valuable for structure-activity relationship studies and prioritizing compounds for drug development.

NPASS serves as a bridge between natural product chemistry and pharmacology, enabling researchers to identify bioactive scaffolds from specific taxonomic groups and explore the relationship between chemical structure and biological activity.

## Key Statistics

| Metric | Value |
|--------|-------|
| Natural Products | 35,000+ |
| Quantitative Activities | 460,000+ |
| Source Species | 25,000+ |
| Protein Targets | 5,000+ |
| Last Update | 2023 |

## Primary Use Cases

1. **Bioactivity Mining** - Find natural products with specific activity profiles
2. **Species-Activity Correlation** - Identify organisms producing bioactive compounds
3. **Target Fishing** - Discover potential targets for natural product hits
4. **Lead Optimization** - Explore SAR using quantitative activity data
5. **Chemotaxonomy** - Study chemical patterns across taxonomic groups

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| NPASS ID | NPC + digits | NPC12345 |
| InChI Key | 27 characters | Standard format |
| NCBI Taxon ID | Integer | 3702 |
| Target Name | Text | Cyclooxygenase-2 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://bidd.group/NPASS/ | Search and browse |
| Bulk Download | http://bidd.group/NPASS/download.html | SDF, TSV formats |
| API | Limited | Contact for access |

## Data Categories

| Category | Description |
|----------|-------------|
| Compound Structure | SMILES, InChI, molecular properties |
| Bioactivity | IC50, EC50, Ki values with units |
| Species Source | Organism taxonomy and occurrence |
| Target Information | Protein targets with UniProt links |

## Limitations

- Commercial use requires contact with database maintainers
- Limited API access compared to other databases
- Last update was in 2023; may not include recent compounds
- Bioactivity data quality varies by source

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [COCONUT](../coconut/_index.md) - Comprehensive natural products
- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Bioactivity data
- [BindingDB](../../2.7.compound.target.interactions/bindingdb/_index.md) - Binding affinities

## References

1. Zeng X, et al. (2018) "NPASS: natural product activity and species source database for natural product research, discovery and tool development." Nucleic Acids Res. 46(D1):D1217-D1222.
