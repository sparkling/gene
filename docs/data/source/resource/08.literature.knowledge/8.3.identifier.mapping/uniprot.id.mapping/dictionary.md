# UniProt ID Mapping - Data Dictionary

## Overview

This data dictionary documents UniProt ID Mapping service.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | uniprot.id.mapping |
| **Name** | UniProt ID Mapping |
| **Total Fields** | 10+ |
| **Last Updated** | 2026-01 |

---

## Supported ID Types

| Database | ID Example | Direction |
|----------|------------|-----------|
| UniProtKB AC | P38398 | From/To |
| UniProtKB ID | BRCA1_HUMAN | From/To |
| Gene Name | BRCA1 | From |
| NCBI Gene ID | 672 | From/To |
| RefSeq Protein | NP_009225 | From/To |
| Ensembl | ENSG00000012048 | From/To |
| PDB | 1JM7 | To |
| HGNC | HGNC:1100 | From/To |

---

## Response Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| from | string | Input identifier |
| to.primaryAccession | string | UniProt accession |
| to.uniProtkbId | string | Entry name |
| to.organism | object | Species info |
| to.proteinDescription | object | Protein name |
| to.genes | array | Gene names |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| AC | Accession |
| ID | Identifier/Entry Name |
| KB | KnowledgeBase |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
