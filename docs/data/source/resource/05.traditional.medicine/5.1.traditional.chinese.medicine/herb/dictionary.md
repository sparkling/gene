# HERB Database - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | herb |
| **Name** | HERB Database |
| **Total Fields** | 32 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| herb_id | String | Yes | HERB database identifier | HERB001234 |
| herb_name | String | Yes | Herb common name | Ginseng |
| chinese_name | String | No | Chinese name | 人参 |
| pinyin_name | String | No | Pinyin romanization | renshen |
| latin_name | String | No | Latin botanical name | Panax ginseng C.A.Mey. |
| ingredient_id | String | No | Ingredient identifier | HBIN001234 |
| ingredient_name | String | No | Chemical compound name | Ginsenoside Rg1 |
| target_id | String | No | Target identifier | HBTAR001234 |
| gene_name | String | No | Target gene symbol | TP53 |

---

## Integration Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| pubchem_cid | Integer | PubChem Compound ID | 441923 |
| drugbank_id | String | DrugBank ID | DB00000 |
| kegg_compound_id | String | KEGG Compound ID | C00001 |
| uniprot_id | String | UniProt accession | P04637 |
| omim_id | String | OMIM disease ID | 612555 |
| disease_name | String | Associated disease | Cancer |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Herb ID | HERB###### | Herb identifier | HERB001234 |
| Ingredient ID | HBIN###### | Compound identifier | HBIN001234 |
| Target ID | HBTAR###### | Target identifier | HBTAR001234 |
| Disease ID | HBDIS###### | Disease identifier | HBDIS001234 |

---

## Entity Relationships

### Herb to Ingredient
- **Cardinality:** 1:N
- **Description:** Herbs contain multiple ingredients
- **Key Fields:** herb_id, ingredient_id

### Ingredient to Target
- **Cardinality:** N:M
- **Description:** Ingredients affect multiple targets
- **Key Fields:** ingredient_id, target_id

### Target to Disease
- **Cardinality:** N:M
- **Description:** Targets associated with diseases
- **Key Fields:** target_id, disease_name

---

## Cross-Reference Databases

| Database | ID Type | Description |
|----------|---------|-------------|
| PubChem | CID | Compound structures |
| DrugBank | DB ID | Drug information |
| KEGG | C ID | Compound/pathway data |
| UniProt | Accession | Protein information |
| OMIM | MIM | Genetic disease data |
| TTD | TTDID | Therapeutic targets |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HERB | High-throughput Experiment- and Reference-guided Database | Full name |
| TCM | Traditional Chinese Medicine | Medical system |
| HBI | HERB Ingredient | Ingredient prefix |
| TTD | Therapeutic Target Database | External database |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
