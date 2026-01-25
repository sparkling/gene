---
id: gmrepo
title: "GMrepo - Gut Microbiota Repository"
type: data-source
category: microbiome
subcategory: gut.microbiome
parent: ../README.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [gut, microbiome, metagenomics, curated, phenotypes]
---

# GMrepo - Gut Microbiota Repository

**Category:** [Microbiome](../../../README.md) > [Gut Microbiome](../README.md)

## Overview

GMrepo (Gut Microbiota Repository) is a manually curated database of human gut metagenomes. It integrates data from multiple sources and provides standardized annotations including phenotype metadata, taxonomic profiles, and quality metrics.

The database focuses on associating gut microbiome composition with health conditions and phenotypes. Each project and sample is annotated with detailed metadata enabling cross-study comparisons and meta-analyses.

GMrepo is valuable for gut microbiome research, biomarker discovery, and understanding the relationship between intestinal microbiota and human health.

## Key Statistics

| Metric | Value |
|--------|-------|
| Samples | 70,000+ |
| Projects | 400+ |
| Species | 7,000+ |
| Phenotypes | 100+ |
| Last Update | 2023 |

## Primary Use Cases

1. Gut microbiome composition analysis
2. Microbiome-phenotype associations
3. Cross-study meta-analysis
4. Reference taxonomic profiles
5. Biomarker discovery

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Run Accession | `[SED]RR[0-9]+` | SRR1234567 |
| Project Accession | `PRJNA[0-9]+` | PRJNA123456 |
| NCBI Taxon ID | Numeric | 562 (E. coli) |
| Phenotype | Controlled | Healthy, IBD, T2D |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://gmrepo.humangut.info | Search and browse |
| REST API | https://gmrepo.humangut.info/api | Programmatic access |
| Downloads | https://gmrepo.humangut.info/Downloads | Bulk data |

## License

| Aspect | Value |
|--------|-------|
| License | Free for academic use |
| Commercial Use | Contact maintainers |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [MetaHIT](../metahit/README.md) - Reference metagenomes
- [gutMGene](../gutmgene/README.md) - Microbiota-gene links
