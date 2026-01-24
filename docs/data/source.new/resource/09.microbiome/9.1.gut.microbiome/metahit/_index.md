---
id: metahit
title: "MetaHIT - Metagenomics of the Human Intestinal Tract"
type: data-source
category: microbiome
subcategory: gut.microbiome
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [metagenomics, gut, reference, gene-catalog, european]
---

# MetaHIT - Metagenomics of the Human Intestinal Tract

**Category:** [Microbiome](../../../_index.md) > [Gut Microbiome](../_index.md)

## Overview

MetaHIT (Metagenomics of the Human Intestinal Tract) was a major European project that established foundational resources for gut microbiome research. It generated reference gene catalogs, developed analysis methods, and studied associations between gut microbiota and disease.

The project produced landmark publications on gut microbial gene richness, enterotypes, and associations with obesity, IBD, and metabolic conditions. MetaHIT data continues to serve as a reference for gut metagenome studies.

MetaHIT resources are essential for comparative metagenomics, understanding gut microbiome diversity, and as reference data for new studies.

## Key Statistics

| Metric | Value |
|--------|-------|
| Reference Genes | 10M+ (IGC catalog) |
| Samples | 1,267 (original) |
| Countries | 8 European |
| Publications | 60+ |
| Last Update | Project completed (legacy) |

## Primary Use Cases

1. Reference gene catalog queries
2. Metagenomic read mapping
3. Comparative gut microbiome analysis
4. Enterotype classification
5. Method benchmarking

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| IGC Gene ID | Internal | MH0001_GL0000001 |
| Sample ID | Project-specific | MH0001 |
| NCBI SRA | `[SED]RR[0-9]+` | ERR123456 |
| ENA Project | `PRJEB[0-9]+` | PRJEB1220 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Gene Catalog | http://meta.genomics.cn | IGC download |
| ENA | https://www.ebi.ac.uk/ena | Raw sequences |
| Publications | PubMed | Searchable |

## License

| Aspect | Value |
|--------|-------|
| License | Public (Fort Lauderdale) |
| Commercial Use | Yes (with attribution) |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [HMP](../../9.2.body.site.microbiomes/hmp/_index.md) - US Human Microbiome Project
- [GMrepo](../gmrepo/_index.md) - Curated repository
