# BioGRID - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | biogrid |
| **Name** | BioGRID Interactions |
| **Total Fields** | 28 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| INTERACTION_ID | Integer | Yes | Unique BioGRID interaction identifier | 103 |
| OFFICIAL_SYMBOL_A | String | Yes | HGNC symbol for interactor A | TP53 |
| OFFICIAL_SYMBOL_B | String | Yes | HGNC symbol for interactor B | MDM2 |
| EXPERIMENTAL_SYSTEM | String | No | Detection method used | Two-hybrid |
| EXPERIMENTAL_SYSTEM_TYPE | Enum | Yes | Interaction type category | physical, genetic |
| PUBMED_ID | Integer | No | Literature reference | 9054499 |
| ORGANISM_A | Integer | No | NCBI Taxonomy ID for A | 9606 |
| ORGANISM_B | Integer | No | NCBI Taxonomy ID for B | 9606 |
| THROUGHPUT | Enum | No | Experiment scale | High Throughput, Low Throughput |

---

## Identifier Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| ENTREZ_GENE_A | String | NCBI Gene ID for A | 7157 |
| ENTREZ_GENE_B | String | NCBI Gene ID for B | 4193 |
| BIOGRID_A | Integer | BioGRID internal ID for A | 106639 |
| BIOGRID_B | Integer | BioGRID internal ID for B | 108267 |
| SWISS_PROT_A | String | UniProt Swiss-Prot accession for A | P04637 |
| SWISS_PROT_B | String | UniProt Swiss-Prot accession for B | Q00987 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| BioGRID Interaction | Numeric | Unique interaction ID | 103 |
| BioGRID ID | Numeric | Internal protein ID | 106639 |
| Entrez Gene | Numeric | NCBI Gene ID | 7157 |
| NCBI Taxonomy | Numeric | Species ID | 9606 |
| PubMed | Numeric | Publication ID | 9054499 |

---

## Enumerations

### Experimental System Types

| Value | Description |
|-------|-------------|
| physical | Protein-protein physical interaction |
| genetic | Genetic interaction (epistasis) |

### Experimental Systems (Physical)

| Value | Description |
|-------|-------------|
| Two-hybrid | Yeast/mammalian two-hybrid |
| Affinity Capture-MS | Mass spectrometry after purification |
| Affinity Capture-Western | Western blot after purification |
| Co-crystal Structure | X-ray crystallography |
| Co-fractionation | Biochemical co-fractionation |
| Co-localization | Microscopy co-localization |
| Co-purification | Biochemical co-purification |
| FRET | Fluorescence resonance energy transfer |
| PCA | Protein-fragment complementation |
| Reconstituted Complex | In vitro reconstitution |

### Experimental Systems (Genetic)

| Value | Description |
|-------|-------------|
| Synthetic Lethality | Combined deletion lethal |
| Synthetic Growth Defect | Combined deletion impairs growth |
| Synthetic Rescue | Second mutation rescues first |
| Phenotypic Enhancement | Combined mutations enhance phenotype |
| Phenotypic Suppression | Second mutation suppresses first |
| Dosage Lethality | Overexpression lethal with mutation |
| Dosage Rescue | Overexpression rescues mutation |

### Throughput

| Value | Description |
|-------|-------------|
| High Throughput | Large-scale systematic screen |
| Low Throughput | Focused individual study |

---

## Entity Relationships

### Interaction to Proteins
- **Cardinality:** 1:2
- **Description:** Binary interaction between two proteins
- **Key Fields:** INTERACTION_ID, BIOGRID_A, BIOGRID_B

### Interaction to Publication
- **Cardinality:** N:1
- **Description:** Interactions supported by literature
- **Key Fields:** INTERACTION_ID, PUBMED_ID

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| BioGRID | Biological General Repository for Interaction Datasets | Database name |
| MS | Mass Spectrometry | Detection method |
| FRET | Fluorescence Resonance Energy Transfer | Detection method |
| PCA | Protein-fragment Complementation Assay | Detection method |
| Y2H | Yeast Two-Hybrid | Detection method |
| TAB | Tab-delimited format | File format |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
