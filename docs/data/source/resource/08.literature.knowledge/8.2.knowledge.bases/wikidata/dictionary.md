# Wikidata - Data Dictionary

## Overview

This data dictionary documents the schema for Wikidata entity records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | wikidata |
| **Name** | Wikidata |
| **Total Fields** | 10+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| id | string | Yes | Q-identifier | `Q227339` |
| type | string | Yes | Entity type | `item`, `property` |

### Labels and Descriptions

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| labels | object | Names by language |
| descriptions | object | Descriptions by language |
| aliases | object | Alternative names |

### Claims (Statements)

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| claims | object | Property statements |
| claims.P*.mainsnak | object | Main value |
| claims.P*.rank | string | preferred/normal/deprecated |
| claims.P*.references | array | Source references |
| claims.P*.qualifiers | object | Statement qualifiers |

---

## Common Properties (Biomedical)

| Property ID | Label | Description |
|-------------|-------|-------------|
| P31 | instance of | Type classification |
| P703 | found in taxon | Species association |
| P351 | Entrez Gene ID | NCBI Gene identifier |
| P352 | UniProt protein ID | Protein identifier |
| P353 | HGNC gene symbol | Gene symbol |
| P354 | HGNC ID | Gene identifier |
| P486 | MeSH descriptor ID | Medical subject |
| P492 | OMIM ID | Disease/gene |
| P493 | ICD-9-CM | Disease code |
| P494 | ICD-10 | Disease code |
| P2892 | UMLS CUI | Concept identifier |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| WD | Wikidata | Knowledge base |
| QID | Q-Identifier | Item ID (Q12345) |
| PID | Property ID | Property (P12345) |
| SPARQL | SPARQL Protocol and RDF Query Language | Query language |

---

## Cross-References

| Database | Property | Description |
|----------|----------|-------------|
| NCBI Gene | P351 | Gene ID |
| UniProt | P352 | Protein ID |
| OMIM | P492 | Disease/Gene |
| PubMed | P698 | Citations |
| DOI | P356 | Publications |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
