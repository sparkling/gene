# PMC ID Converter - Data Dictionary

## Overview

This data dictionary documents the PMC ID Converter mapping service.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pmc.id.converter |
| **Name** | PMC ID Converter |
| **Total Fields** | 7 |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| pmcid | string | No | PMC identifier | `PMC7654321` |
| pmid | string | No | PubMed identifier | `12345678` |
| doi | string | No | DOI | `10.1038/nature12345` |
| live | boolean | No | Currently in PMC | `true` |
| release_date | string | No | PMC release date | `2022-08-18` |
| errmsg | string | No | Error message | `PMID not found` |

---

## Conversion Directions

| From | To | Description |
|------|-----|-------------|
| PMID | PMCID, DOI | PubMed to PMC/DOI |
| PMCID | PMID, DOI | PMC to PubMed/DOI |
| DOI | PMID, PMCID | DOI to PubMed/PMC |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PMC | PubMed Central | Full-text archive |
| PMID | PubMed Identifier | Citation ID |
| DOI | Digital Object Identifier | Persistent link |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
