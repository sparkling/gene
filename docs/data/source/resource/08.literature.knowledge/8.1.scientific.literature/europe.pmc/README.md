---
id: europe.pmc
title: "Europe PMC"
type: data-source
category: literature
subcategory: scientific.literature
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [literature, publications, open-access, europe, full-text]
---

# Europe PMC

**Category:** [Literature](../../../README.md) > [Scientific Literature](../README.md)

## Overview

Europe PMC is a free, full-text literature database extending beyond PubMed to include patents, clinical guidelines, and research grant abstracts from major European funders. It provides advanced text mining, citation networks, and semantic annotations.

The database integrates content from PubMed, PubMed Central, and additional sources including preprints, patents, and agricultural research. Europe PMC offers unique features like SciLite (text-mined annotations), citation networks, and funder manuscript compliance tracking.

Europe PMC is particularly valuable for European researchers, open access advocacy, and text mining applications.

## Key Statistics

| Metric | Value |
|--------|-------|
| Articles | 44M+ |
| Full-Text Articles | 9M+ |
| Preprints | 750,000+ |
| Patents | 4M+ |
| Last Update | Daily |

## Primary Use Cases

1. Literature search and discovery
2. Text mining and NLP research
3. Citation network analysis
4. Open access compliance tracking
5. Grant-publication linking

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| PMID | Numeric | 12345678 |
| PMCID | `PMC[0-9]+` | PMC1234567 |
| DOI | Standard | 10.1000/xyz123 |
| EuropePMC ID | Internal | PPR123456 (preprints) |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://europepmc.org | Search and browse |
| REST API | https://www.ebi.ac.uk/europepmc/webservices/rest | Programmatic access |
| Annotations API | https://www.ebi.ac.uk/europepmc/annotations_api | Text-mined data |
| OAI-PMH | https://www.ebi.ac.uk/europepmc/oai | Metadata harvesting |
| FTP | https://europepmc.org/ftp | Bulk downloads |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access (varies by article) |
| Commercial Use | Depends on article license |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [PubMed](../pubmed/README.md) - NCBI literature database
- [PubMed Central](../pubmed.central/README.md) - Full-text archive
