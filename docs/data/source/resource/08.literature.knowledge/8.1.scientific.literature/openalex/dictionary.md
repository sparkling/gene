# OpenAlex - Data Dictionary

## Overview

This data dictionary documents the schema for OpenAlex scholarly works.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | openalex |
| **Name** | OpenAlex |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| id | string | Yes | OpenAlex ID | `https://openalex.org/W2741809807` |
| doi | string | No | DOI | `10.1038/nature12345` |

### Bibliographic

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| display_name | string | Yes | Work title |
| publication_year | integer | No | Year published |
| publication_date | string | No | Full date |
| type | string | No | Work type (article, book, etc.) |

### Open Access

| Field Name | Data Type | Description | Allowed Values |
|------------|-----------|-------------|----------------|
| open_access.is_oa | boolean | Open access status | true/false |
| open_access.oa_status | string | OA type | gold, green, hybrid, bronze, closed |
| open_access.oa_url | string | Free access URL | URL |

---

## Authorship Object

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| author_position | string | first, middle, last |
| author.id | string | OpenAlex author ID |
| author.display_name | string | Author name |
| author.orcid | string | ORCID identifier |
| institutions | array | Affiliated institutions |

---

## Concepts (Topics)

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| id | string | Concept ID |
| display_name | string | Concept name |
| level | integer | Hierarchy level (0-5) |
| score | number | Relevance score (0-1) |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| OA | Open Access | Free to read |
| MAG | Microsoft Academic Graph | Predecessor dataset |
| DOI | Digital Object Identifier | Persistent link |

---

## Cross-References

| Database | ID Type | Description |
|----------|---------|-------------|
| DOI | DOI | Publisher link |
| PubMed | PMID | Via external IDs |
| ORCID | ORCID ID | Author identifier |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
