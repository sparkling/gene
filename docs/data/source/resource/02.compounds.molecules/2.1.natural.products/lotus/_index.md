---
id: lotus
title: "LOTUS - Linked Open Total Unified Structure-organism"
type: data-source
category: compounds
subcategory: natural-products
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [natural-products, wikidata, sparql, structure-organism, linked-data]
---

# LOTUS - Linked Open Total Unified Structure-organism

**Category:** [Compounds & Molecules](../../_index.md) > [Natural Products](../_index.md)

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

## License

| Aspect | Value |
|--------|-------|
| License | CC0 (Public Domain) |
| Commercial Use | Yes (unrestricted) |
| Attribution Required | No (but appreciated) |

## Related Resources

- [COCONUT](../coconut/_index.md) - Comprehensive natural products database
- [ChEBI](../../2.6.chemical.ontology.classification/chebi/_index.md) - Chemical ontology
- [PubChem](../../2.6.chemical.ontology.classification/pubchem/_index.md) - Chemical information

## See Also

- [Schema Documentation](./schema.md)

## References

1. Rutz A, et al. (2022) "The LOTUS initiative for open knowledge management in natural products research." eLife 11:e70780. DOI: 10.7554/eLife.70780
