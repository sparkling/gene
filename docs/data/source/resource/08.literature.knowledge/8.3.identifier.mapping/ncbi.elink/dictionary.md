# NCBI ELink - Data Dictionary

## Overview

This data dictionary documents the NCBI ELink cross-database linking service.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | ncbi.elink |
| **Name** | NCBI ELink |
| **Total Fields** | 8 |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| dbfrom | string | Source database |
| id | string | Source record ID |
| linksets | array | Link result sets |
| linksetdbs | array | Scored link results |
| cmd | string | Command type |

---

## Supported Databases

| Database | Description |
|----------|-------------|
| pubmed | PubMed citations |
| gene | NCBI Gene |
| protein | NCBI Protein |
| nuccore | Nucleotide sequences |
| clinvar | Clinical variants |
| omim | OMIM |
| mesh | MeSH terms |
| biosample | BioSamples |
| sra | Sequence Read Archive |

---

## Link Commands

| Command | Description |
|---------|-------------|
| neighbor | Related records |
| neighbor_score | With relevance scores |
| acheck | All available links |
| llinks | LinkOut URLs |
| prlinks | Provider links |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
