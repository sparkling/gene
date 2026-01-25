# Wikipedia - Data Dictionary

## Overview

This data dictionary documents the schema for Wikipedia biomedical articles.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | wikipedia |
| **Name** | Wikipedia |
| **Total Fields** | 15+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| pageid | integer | Yes | Wikipedia page ID | `12345678` |
| title | string | Yes | Article title | `BRCA1` |
| wikidata_id | string | No | Wikidata Q-ID | `Q227339` |

### Content

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| extract | string | Plain text summary |
| content | string | Full wikitext/HTML |
| categories | array | Category membership |
| infobox | object | Structured infobox data |

### References

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| references | array | Citation list |
| references.pmid | string | PubMed ID in reference |
| references.doi | string | DOI in reference |
| external_links | array | External URLs |

---

## Infobox Fields (Gene Articles)

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| Name | string | Gene name |
| Symbol | string | Gene symbol |
| HGNC | string | HGNC ID |
| EntrezGene | string | NCBI Gene ID |
| RefSeq | string | RefSeq accession |
| UniProt | string | UniProt ID |
| Chromosome | string | Chromosomal location |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| WP | Wikipedia | Encyclopedia |
| WMF | Wikimedia Foundation | Organization |
| CC | Creative Commons | License family |

---

## Cross-References

| Database | ID Type | Description |
|----------|---------|-------------|
| Wikidata | Q-ID | Structured data link |
| PubMed | PMID | In references |
| NCBI Gene | Gene ID | In infoboxes |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
