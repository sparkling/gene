---
id: openalex
title: "OpenAlex"
type: data-source
category: literature
subcategory: scientific.literature
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [literature, open-data, knowledge-graph, citations, metadata]
---

# OpenAlex

**Category:** [Literature](../../../_index.md) > [Scientific Literature](../_index.md)

## Overview

OpenAlex is a fully open catalog of the global research system, created by OurResearch as an open alternative to proprietary scholarly databases. It indexes scholarly works, authors, institutions, venues, concepts, and their relationships.

Built as the successor to Microsoft Academic Graph, OpenAlex provides a comprehensive knowledge graph of scholarly literature with rich metadata, citation links, and concept tagging. All data is freely available via API and bulk download.

OpenAlex is particularly valuable for bibliometric analysis, research assessment, and building open scholarly infrastructure.

## Key Statistics

| Metric | Value |
|--------|-------|
| Works | 250M+ |
| Authors | 90M+ |
| Institutions | 110,000+ |
| Concepts | 65,000+ |
| Last Update | Daily |

## Primary Use Cases

1. Bibliometric analysis
2. Research impact assessment
3. Author disambiguation
4. Institution analytics
5. Open scholarly infrastructure development

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| OpenAlex ID | `W[0-9]+` (works) | W2741809807 |
| Author ID | `A[0-9]+` | A5023888391 |
| Institution ID | `I[0-9]+` | I136199984 |
| DOI | Standard | 10.1000/xyz123 |
| ORCID | Standard | 0000-0002-1825-0097 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://openalex.org | Search and explore |
| REST API | https://api.openalex.org | Free, rate-limited |
| Snapshot | https://docs.openalex.org/download-snapshot | Bulk download |
| Data Dump | S3/AWS | Full database |

## License

| Aspect | Value |
|--------|-------|
| License | CC0 (Public Domain) |
| Commercial Use | Yes |
| Attribution | Appreciated but not required |

## See Also

- [Schema Documentation](./schema.md)
- [Semantic Scholar](../semantic.scholar/_index.md) - AI-powered discovery
- [PubMed](../pubmed/_index.md) - Biomedical focus
