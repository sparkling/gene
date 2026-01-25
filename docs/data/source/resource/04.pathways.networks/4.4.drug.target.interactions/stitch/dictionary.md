# STITCH - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | stitch |
| **Name** | STITCH Chemical-Protein Interactions |
| **Total Fields** | 22 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| chemical | String | Yes | STITCH chemical ID (CIDm/CIDs format) | CIDm00002244 |
| protein | String | Yes | STRING protein ID (taxid.ENSP format) | 9606.ENSP00000269305 |
| preferredName_chemical | String | No | Chemical common name | aspirin |
| preferredName_protein | String | No | Gene symbol | PTGS2 |
| ncbiTaxonId | Integer | No | NCBI taxonomy ID | 9606 |
| combined_score | Integer | Yes | Combined confidence score (0-1000) | 900 |

---

## Evidence Channel Scores

| Field Name | Data Type | Description | Range |
|------------|-----------|-------------|-------|
| experimental | Integer | Binding assays, crystal structures | 0-1000 |
| prediction | Integer | Binding site similarity, docking | 0-1000 |
| database | Integer | DrugBank, KEGG, pathway DBs | 0-1000 |
| textmining | Integer | Literature co-occurrence | 0-1000 |
| experimental_direct | Integer | Direct experimental evidence | 0-1000 |
| experimental_transferred | Integer | Evidence from homologs | 0-1000 |
| database_direct | Integer | Direct database annotation | 0-1000 |
| database_transferred | Integer | Transferred annotation | 0-1000 |
| textmining_direct | Integer | Direct text mining | 0-1000 |
| textmining_transferred | Integer | Transferred text mining | 0-1000 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| STITCH Chemical | CID[ms]######## | PubChem-based with stereo prefix | CIDm00002244 |
| STRING Protein | {taxid}.ENSP#### | STRING protein format | 9606.ENSP00000269305 |
| PubChem CID | Numeric | Compound ID | 2244 |
| Ensembl Protein | ENSP########### | Protein ID | ENSP00000269305 |
| NCBI Taxonomy | Numeric | Species ID | 9606 |

### STITCH Chemical ID Prefixes

| Prefix | Description |
|--------|-------------|
| CIDm | Stereo-specific compound (major form) |
| CIDs | Stereo-specific compound (specific stereoisomer) |
| CID0 | Flat/2D structure (no stereochemistry) |

---

## Enumerations

### Action Modes

| Value | Description |
|-------|-------------|
| activation | Chemical activates protein |
| inhibition | Chemical inhibits protein |
| binding | Chemical binds protein |
| catalysis | Enzymatic relationship |
| reaction | Biochemical reaction |
| expression | Affects expression |
| ptmod | Post-translational modification |

---

## Score Interpretation

| Score Range | Confidence Level | Description |
|-------------|-----------------|-------------|
| >= 900 | Highest | Very high confidence |
| 700 - 899 | High | High confidence |
| 400 - 699 | Medium | Medium confidence |
| 150 - 399 | Low | Low confidence |
| < 150 | Lowest | Very low confidence |

---

## Entity Relationships

### Chemical-Protein Interaction
- **Cardinality:** N:M
- **Description:** Chemical compounds interact with proteins
- **Key Fields:** chemical, protein, combined_score

### Chemical to Actions
- **Cardinality:** 1:N
- **Description:** Multiple action modes per interaction
- **Key Fields:** chemical, protein, actions

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| STITCH | Search Tool for Interactions of Chemicals | Database name |
| CID | Compound Identifier | PubChem format |
| ENSP | Ensembl Protein | Protein identifier prefix |
| CPI | Chemical-Protein Interaction | Interaction type |
| DTI | Drug-Target Interaction | Alternative term |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
