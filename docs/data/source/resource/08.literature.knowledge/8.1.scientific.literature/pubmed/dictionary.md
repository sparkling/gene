# PubMed - Data Dictionary

## Overview

This data dictionary documents the schema for PubMed biomedical literature citations.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pubmed |
| **Name** | PubMed |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| pmid | integer | Yes | PubMed unique identifier | `12345678` |
| doi | string | No | Digital Object Identifier | `10.1038/nature12345` |
| pmc_id | string | No | PubMed Central ID | `PMC7654321` |

### Bibliographic

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| title | string | Yes | Article title | `Gene expression in cancer` |
| abstract | string | No | Abstract text | Full abstract content |
| authors | array | No | Author list | Array of author objects |
| journal.title | string | No | Journal name | `Nature` |
| journal.iso_abbreviation | string | No | ISO abbreviation | `Nat.` |
| pagination | string | No | Page numbers | `123-145` |

### Classification

| Field Name | Data Type | Required | Description | Allowed Values |
|------------|-----------|----------|-------------|----------------|
| publication_types | array | No | Article types | Journal Article, Review, Clinical Trial |
| mesh_terms | array | No | MeSH descriptors | Array with descriptor_ui (D-prefixed) |
| keywords | array | No | Author keywords | Free-text keywords |

---

## Author Object

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| last_name | string | Yes | Family name |
| fore_name | string | No | Given name |
| initials | string | No | Name initials |
| affiliation | string | No | Institution |
| orcid | string | No | ORCID identifier |

---

## MeSH Term Object

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| descriptor_name | string | MeSH term name |
| descriptor_ui | string | MeSH unique identifier (D-prefixed) |
| major_topic | boolean | True if major topic of article |
| qualifiers | array | Subheading qualifiers |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PMID | PubMed Identifier | Unique 8-digit integer |
| PMC | PubMed Central | Full-text archive |
| DOI | Digital Object Identifier | Persistent link |
| MeSH | Medical Subject Headings | NLM controlled vocabulary |
| ORCID | Open Researcher and Contributor ID | Author identifier |
| MEDLINE | Medical Literature Analysis and Retrieval System Online | NLM database |

---

## Cross-References

| Database | ID Type | Description |
|----------|---------|-------------|
| PMC | PMC ID | Full-text article |
| DOI | DOI | Publisher link |
| MeSH | Descriptor UI | Subject classification |
| ORCID | ORCID ID | Author identification |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
