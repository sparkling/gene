# Semantic Scholar - Data Dictionary

## Overview

This data dictionary documents the schema for Semantic Scholar paper records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | semantic.scholar |
| **Name** | Semantic Scholar |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| paperId | string | Yes | S2 unique ID (40-char hex) | `649def34f8be52c8b66281af98ae884c09aef38b` |
| corpusId | integer | No | Corpus ID | `12345678` |
| externalIds.DOI | string | No | DOI | `10.1038/nature12345` |
| externalIds.PubMed | string | No | PMID | `12345678` |
| externalIds.ArXiv | string | No | arXiv ID | `2201.12345` |

### Bibliographic

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| title | string | Yes | Paper title |
| abstract | string | No | Abstract text |
| venue | string | No | Publication venue |
| year | integer | No | Publication year |

### Citation Metrics

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| citationCount | integer | Total citations |
| referenceCount | integer | References cited |
| influentialCitationCount | integer | Influential citations |

---

## AI-Generated Content

### TLDR (Summary)

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| tldr.model | string | Model used (e.g., "tldr@v2.0.0") |
| tldr.text | string | AI-generated summary |

### Embedding

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| embedding.model | string | Model name |
| embedding.vector | array | 768-dim vector |

---

## Fields of Study

| Category | Description |
|----------|-------------|
| Medicine | Medical research |
| Biology | Life sciences |
| Computer Science | CS research |
| Chemistry | Chemical sciences |
| Physics | Physical sciences |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| S2 | Semantic Scholar | AI2 literature database |
| TLDR | Too Long; Didn't Read | AI summary |
| MAG | Microsoft Academic Graph | Predecessor |
| AI2 | Allen Institute for AI | Organization |

---

## Cross-References

| Database | ID Type | Description |
|----------|---------|-------------|
| PubMed | PMID | via externalIds |
| DOI | DOI | via externalIds |
| arXiv | arXiv ID | via externalIds |
| DBLP | DBLP Key | CS papers |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
