# PubMed Central - Data Dictionary

## Overview

This data dictionary documents the schema for PubMed Central full-text articles.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pubmed.central |
| **Name** | PubMed Central (PMC) |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| pmc_id | string | Yes | PMC unique identifier | `PMC7654321` |
| pmid | integer | No | PubMed identifier | `12345678` |
| doi | string | No | Digital Object Identifier | `10.1038/nature12345` |

### Content

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| title | string | Yes | Article title |
| abstract | string | No | Abstract text |
| body.sections | array | No | Full-text sections |

### Bibliographic

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| authors | array | Author list with affiliations |
| journal.title | string | Journal name |
| journal.publisher | string | Publisher name |
| pub_date | object | Publication date |

---

## License Object

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| type | string | License type (CC-BY, CC-BY-NC, etc.) |
| text | string | License statement |
| url | string | License URL |

---

## Supplementary Materials

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| id | string | Supplementary item identifier |
| label | string | Display label |
| media_type | string | MIME type |
| url | string | Download URL |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PMC | PubMed Central | Full-text archive |
| OA | Open Access | Freely available articles |
| JATS | Journal Article Tag Suite | XML format |
| NXML | NLM XML | Article format |

---

## Cross-References

| Database | ID Type | Description |
|----------|---------|-------------|
| PubMed | PMID | Citation metadata |
| DOI | DOI | Publisher link |
| ORCID | ORCID ID | Author identifier |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
