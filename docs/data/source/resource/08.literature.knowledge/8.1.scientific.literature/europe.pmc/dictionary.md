# Europe PMC - Data Dictionary

## Overview

This data dictionary documents the schema for Europe PMC scientific literature.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | europe.pmc |
| **Name** | Europe PMC |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| id | string | Yes | Europe PMC identifier | `MED/12345678` |
| source | string | Yes | Source database | `MED`, `PMC`, `PPR` |
| pmid | string | No | PubMed ID | `12345678` |
| pmcid | string | No | PMC ID | `PMC7654321` |
| doi | string | No | DOI | `10.1038/nature12345` |

### Source Database Codes

| Code | Description |
|------|-------------|
| MED | MEDLINE/PubMed |
| PMC | PubMed Central |
| PAT | Patents |
| AGR | Agricola |
| CBA | Chinese Biological Abstracts |
| ETH | EThOS UK theses |
| HIR | NHS Health Research |
| CTX | ClinicalTrials.gov |
| PPR | Preprints |

### Bibliographic

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| title | string | Article title |
| authorString | string | Formatted author list |
| abstractText | string | Abstract text |
| journalInfo | object | Journal details |

---

## Annotations (Text Mining)

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| annotations.genes | array | Gene mentions |
| annotations.diseases | array | Disease mentions |
| annotations.chemicals | array | Chemical mentions |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| EPMC | Europe PMC | European literature database |
| EMBL-EBI | European Molecular Biology Laboratory - European Bioinformatics Institute | Host |
| OA | Open Access | Free full text |
| PPR | Preprints | Non-peer-reviewed |

---

## Cross-References

| Database | ID Type | Description |
|----------|---------|-------------|
| PubMed | PMID | Citation link |
| PMC | PMCID | Full text |
| UniProt | Accession | Protein annotations |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
