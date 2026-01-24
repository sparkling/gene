---
id: ncbi.elink
title: "NCBI E-Link"
type: data-source
category: literature
subcategory: identifier.mapping
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [identifiers, mapping, ncbi, entrez, cross-reference]
---

# NCBI E-Link

**Category:** [Literature](../../../_index.md) > [Identifier Mapping](../_index.md)

## Overview

NCBI E-Link (ELink) is part of the Entrez E-utilities suite that discovers links between records in different NCBI databases. It can find connections between genes, proteins, literature, structures, and other biological data stored across NCBI resources.

ELink supports several link types including neighbor links (similar records), computational links (sequence similarity), and curated links (manually annotated relationships). It enables powerful cross-database queries within the NCBI ecosystem.

ELink is essential for programmatic integration of NCBI databases and discovering relationships between different types of biological data.

## Key Statistics

| Metric | Value |
|--------|-------|
| Linked Databases | 40+ NCBI databases |
| Link Types | Neighbor, Related, Computed |
| Query Rate | 3 requests/second |
| API Key Rate | 10 requests/second |
| Last Update | Real-time |

## Primary Use Cases

1. Cross-database record linking
2. Finding related records
3. Literature-gene associations
4. Protein-structure connections
5. Batch ID conversion within NCBI

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Database Name | NCBI db codes | pubmed, gene, protein |
| UID | Database-specific | 12345678 (PMID) |
| Link Name | `{db}_{db}` | pubmed_gene |
| API Key | String | Optional for higher rate |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| E-Link | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi | Core API |
| EDirect | Command-line | elink command |
| Documentation | https://www.ncbi.nlm.nih.gov/books/NBK25499 | Full guide |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [PubMed](../../8.1.scientific.literature/pubmed/_index.md) - Literature database
- [PMC ID Converter](../pmc.id.converter/_index.md) - Publication IDs
