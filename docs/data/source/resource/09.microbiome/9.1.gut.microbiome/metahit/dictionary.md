# MetaHIT - Data Dictionary

## Overview

This data dictionary documents MetaHIT Integrated Gene Catalog (IGC) entries.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | metahit |
| **Name** | MetaHIT |
| **Total Fields** | 15+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| gene_id | string | Yes | IGC gene identifier (MH####_GL#######) |
| sequence | string | No | Nucleotide sequence |
| length | integer | Yes | Gene length in base pairs |
| taxonomy | object | No | Taxonomic assignment |
| kegg_ko | string | No | KEGG Ortholog ID |
| eggnog | string | No | eggNOG annotation |
| cog | string | No | COG category |

---

## Gene ID Format

```
Format: MH{sample_number}_GL{gene_number}
Example: MH0001_GL0000001

Components:
- MH: MetaHIT prefix
- {sample_number}: 4-digit sample origin
- GL: Gene Locus
- {gene_number}: 7-digit gene index
```

---

## Enterotypes

| Enterotype | Dominant Genus | Diet Association |
|------------|---------------|------------------|
| Type 1 | Bacteroides | High protein/fat |
| Type 2 | Prevotella | High fiber/carb |
| Type 3 | Ruminococcus | Variable, Firmicutes |

---

## eggNOG Categories

| Category | Description |
|----------|-------------|
| J | Translation |
| K | Transcription |
| L | Replication |
| C | Energy metabolism |
| G | Carbohydrate metabolism |
| E | Amino acid metabolism |
| S | Function unknown |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| MetaHIT | Metagenomics of the Human Intestinal Tract |
| IGC | Integrated Gene Catalog |
| COG | Clusters of Orthologous Groups |
| KO | KEGG Ortholog |
| eggNOG | evolutionary genealogy of genes: Non-supervised Orthologous Groups |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
