# STRING - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | string |
| **Name** | STRING Protein Interactions |
| **Total Fields** | 18 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| stringId_A | String | Yes | STRING ID for protein A (taxid.ENSP format) | 9606.ENSP00000269305 |
| stringId_B | String | Yes | STRING ID for protein B | 9606.ENSP00000417404 |
| preferredName_A | String | No | Gene symbol for protein A | TP53 |
| preferredName_B | String | No | Gene symbol for protein B | MDM2 |
| ncbiTaxonId | Integer | No | NCBI taxonomy identifier | 9606 |
| score | Number | Yes | Combined confidence score (0-1) | 0.999 |

---

## Evidence Channel Scores

| Field Name | Data Type | Description | Range |
|------------|-----------|-------------|-------|
| nscore | Number | Neighborhood score - gene proximity in genome | 0-1 |
| fscore | Number | Fusion score - gene fusion events | 0-1 |
| pscore | Number | Phylogenetic profile score - co-occurrence | 0-1 |
| ascore | Number | Co-expression score - expression correlation | 0-1 |
| escore | Number | Experimental score - laboratory evidence | 0-1 |
| dscore | Number | Database score - curated pathway sources | 0-1 |
| tscore | Number | Text-mining score - literature co-mentions | 0-1 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| STRING ID | {taxid}.ENSP{numbers} | Composite identifier | 9606.ENSP00000269305 |
| Ensembl Protein | ENSP########### | Ensembl protein ID | ENSP00000269305 |
| NCBI Taxonomy | Numeric | Species identifier | 9606 |
| Gene Symbol | Alphanumeric | HGNC gene symbol | TP53 |

---

## Enumerations

### Enrichment Categories

| Value | Description |
|-------|-------------|
| Process | GO Biological Process |
| Function | GO Molecular Function |
| Component | GO Cellular Component |
| KEGG | KEGG Pathway |
| Pfam | Protein family |
| InterPro | Protein domain |
| SMART | Domain annotation |
| Reactome | Reactome pathway |
| WikiPathways | WikiPathways pathway |
| HPO | Human Phenotype Ontology |
| DisGeNET | Disease-gene associations |

---

## Score Interpretation

| Score Range | Confidence Level | Description |
|-------------|-----------------|-------------|
| >= 0.9 | Highest | Very high confidence |
| 0.7 - 0.9 | High | High confidence |
| 0.4 - 0.7 | Medium | Medium confidence |
| 0.15 - 0.4 | Low | Low confidence |
| < 0.15 | Lowest | Very low confidence |

---

## Entity Relationships

### Protein-Protein Interaction
- **Cardinality:** N:M
- **Description:** Symmetric interaction between two proteins
- **Key Fields:** stringId_A, stringId_B, score

### Protein to Enrichment
- **Cardinality:** 1:N
- **Description:** Proteins associated with functional terms
- **Key Fields:** stringId, enrichment

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| STRING | Search Tool for Retrieval of Interacting Genes/Proteins | Database name |
| ENSP | Ensembl Protein | Protein identifier prefix |
| PPI | Protein-Protein Interaction | Interaction type |
| GO | Gene Ontology | Functional annotation |
| FDR | False Discovery Rate | Statistical correction |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
