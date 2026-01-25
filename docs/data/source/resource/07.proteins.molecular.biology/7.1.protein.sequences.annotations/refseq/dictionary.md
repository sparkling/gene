# RefSeq - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | refseq |
| **Name** | NCBI Reference Sequence Database |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| accession | string | 1:1 | Yes | RefSeq ID with version | `NP_000537.3` |
| gi | integer | 1:1 | No | GenInfo ID (deprecated) | `4557757` |
| gene | string | 1:1 | No | Gene symbol | `TP53` |
| gene_id | integer | 1:1 | No | NCBI Gene ID | `7157` |
| definition | string | 1:1 | Yes | Protein description | `cellular tumor antigen p53` |

### Sequence

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| sequence | string | 1:1 | Yes | Amino acid sequence | `MEEPQSDPSV...` |
| length | integer | 1:1 | Yes | Sequence length | `393` |
| coded_by | string | 1:1 | No | Source transcript | `NM_000546.6:203..1384` |

### Organism

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| organism | string | 1:1 | Yes | Scientific name | `Homo sapiens` |
| tax_id | integer | 1:1 | No | NCBI Taxonomy ID | `9606` |
| taxonomy | array | 1:N | No | Taxonomic lineage | `Eukaryota, Metazoa...` |

---

## Enumerations

### Accession Prefixes

| Prefix | Description | Status |
|--------|-------------|--------|
| NP_ | Protein (from NM_) | Curated |
| XP_ | Predicted protein | Model |
| YP_ | Protein (bacterial/viral) | Curated |
| WP_ | Non-redundant protein | Prokaryotic |
| AP_ | Annotated on AC_ | Third-party |

### Annotation Status

| Value | Description |
|-------|-------------|
| VALIDATED | Full manual review |
| REVIEWED | Curator reviewed |
| PROVISIONAL | Automated + limited review |
| PREDICTED | Computational only |
| MODEL | Gene model based |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| RefSeq | Reference Sequence | Database name |
| NCBI | National Center for Biotechnology Information | Provider |
| GI | GenInfo Identifier | Deprecated |
| CDS | Coding Sequence | Gene region |
| MANE | Matched Annotation from NCBI and EBI | Standard |
| CCDS | Consensus CDS | Curated coding sequences |

---

## Data Quality Notes

1. **Curated vs predicted**: NP_ records are curated; XP_ are predicted
2. **Version tracking**: Version number increments with sequence changes
3. **E-utilities access**: Use Entrez API for programmatic access
4. **Public domain**: US government work, freely available

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
