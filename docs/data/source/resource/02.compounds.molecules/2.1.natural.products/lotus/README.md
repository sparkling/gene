---
id: lotus
title: "LOTUS - Linked Open Total Unified Structure-organism"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: natural.products
tags:
  - natural-products
  - wikidata
  - sparql
  - structure-organism
  - linked-data
---

# LOTUS - Linked Open Total Unified Structure-organism

## Overview

LOTUS (Linked Open Total Unified Structure-organism) is the largest publicly available database linking natural product chemical structures to their biological source organisms with literature evidence. The database is fully integrated into Wikidata, making it accessible through SPARQL queries and enabling seamless integration with the broader Wikidata knowledge graph.

Each entry in LOTUS represents a documented structure-organism-reference triplet: a chemical structure (identified by InChIKey), found in a specific organism (linked to NCBI Taxonomy), with supporting literature citation. This rigorous provenance tracking enables researchers to trace the origin of natural product occurrence data back to primary literature.

LOTUS follows a single-source-of-truth model with PostgreSQL as the canonical database, synchronized to both Wikidata and the LNPN (LOTUS Natural Products Network) web interface. The CC0 license makes it freely usable for any purpose without restrictions.

## Key Statistics

| Metric | Value |
|--------|-------|
| Structure-Organism Pairs | 750,000+ |
| Unique Chemical Structures | 290,000+ |
| Distinct Organisms | 40,000+ |
| Literature References | 75,000+ |
| Validation Accuracy | 97% true positives |
| Last Update | Continuous (Wikidata) |

## Primary Use Cases

1. **Natural Product Origin Verification** - Validate which organisms produce specific compounds
2. **Chemotaxonomy Research** - Study chemical patterns across taxonomic groups
3. **Biosynthetic Pathway Discovery** - Identify organisms with related compound production
4. **Literature Mining** - Access curated references for structure-organism associations
5. **Wikidata Integration** - Federated queries combining LOTUS with other knowledge bases

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Wikidata QID | Q + digits | Q312879 |
| InChI Key | 27 characters | KZJWDPNRJALLNS-VJSFXXLFSA-N |
| NCBI Taxon ID | Integer | 3702 |
| Open Tree of Life ID | Integer | 1054861 |
| DOI | 10.xxxx/... | 10.1016/j.phytochem.2018.01.001 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| SPARQL Query | https://query.wikidata.org | Primary access via Wikidata |
| Web Interface | https://lotus.naturalproducts.net | LNPN interactive search |
| Zenodo Archive | https://doi.org/10.5281/zenodo.5794106 | Bulk TSV/JSON exports |
| Substructure Search | https://idsm.elixir-czech.cz/sparql/endpoint/wikidata | IDSM/Sachem endpoint |

## Core Wikidata Properties

| Property | ID | Description |
|----------|-----|-------------|
| InChIKey | P235 | Chemical structure identifier |
| found in taxon | P703 | Links compound to organism |
| stated in | P248 | Literature reference |
| NCBI taxon ID | P685 | Organism taxonomy |

## Limitations

- Data inherited from source databases may contain errors
- Some structure-organism pairs may be computationally inferred
- Wikidata query interface has rate limits for large queries
- Not all entries have validated taxonomic identifiers

## Related Resources

- [COCONUT](../coconut/README.md) - Comprehensive natural products database
- [ChEBI](../../2.6.chemical.ontology.classification/chebi/README.md) - Chemical ontology
- [PubChem](../../2.6.chemical.ontology.classification/pubchem/README.md) - Chemical information

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## References

1. Rutz A, et al. (2022) "The LOTUS initiative for open knowledge management in natural products research." eLife 11:e70780. DOI: 10.7554/eLife.70780
