---
id: uniprot.id.mapping
title: "UniProt ID Mapping Service"
type: data-source
category: literature
subcategory: identifier.mapping
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [identifiers, proteins, genes, mapping, cross-reference]
---

# UniProt ID Mapping Service

**Category:** [Literature](../../../README.md) > [Identifier Mapping](../README.md)

## Overview

The UniProt ID Mapping service converts identifiers between UniProt accessions and over 100 external database identifiers. It supports mapping to/from gene databases (NCBI Gene, Ensembl), protein databases (RefSeq, PDB), pathway databases (KEGG, Reactome), and many others.

The service handles both individual queries and batch conversions of thousands of identifiers. It is implemented as both an interactive web tool and a programmatic REST API.

UniProt ID Mapping is the primary tool for converting between protein and gene identifiers across major biological databases.

## Key Statistics

| Metric | Value |
|--------|-------|
| External Databases | 100+ |
| Batch Size | Up to 100,000 IDs |
| Response Formats | TSV, FASTA, JSON |
| Mapping Direction | Bidirectional |
| Last Update | Monthly |

## Primary Use Cases

1. Converting gene IDs to protein IDs
2. Mapping across database systems
3. Batch identifier standardization
4. Cross-database data integration
5. Annotation pipeline support

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| UniProt AC | Standard | P04637 |
| Gene ID | Numeric | 7157 |
| Ensembl Gene | `ENSG[0-9]+` | ENSG00000141510 |
| RefSeq Protein | `[NXY]P_[0-9]+` | NP_000537.3 |
| PDB ID | `[0-9][A-Z0-9]{3}` | 1TUP |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.uniprot.org/id-mapping | Interactive |
| REST API | https://rest.uniprot.org/idmapping | Programmatic |
| Documentation | https://www.uniprot.org/help/id_mapping | Full guide |

## License

| Aspect | Value |
|--------|-------|
| License | Creative Commons Attribution 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [UniProt](../../../07.proteins.molecular.biology/7.1.protein.sequences.annotations/uniprot/README.md) - Full database
- [NCBI E-Link](../ncbi.elink/README.md) - NCBI linking
