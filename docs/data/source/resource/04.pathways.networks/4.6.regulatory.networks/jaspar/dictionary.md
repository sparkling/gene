# JASPAR - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | jaspar |
| **Name** | JASPAR Transcription Factor Binding Profiles |
| **Total Fields** | 24 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| matrix_id | String | Yes | JASPAR matrix identifier | MA0106.3 |
| name | String | Yes | Transcription factor name | TP53 |
| version | Integer | No | Matrix version number | 3 |
| collection | Enum | No | JASPAR collection | CORE |
| tax_group | String | No | Taxonomic group | vertebrates |
| class | String | No | TF structural class | Zinc-coordinating |
| family | String | No | TF family | C2H2 zinc finger factors |
| consensus | String | No | Consensus binding sequence | RRRCWWGYYY |

---

## Matrix Data Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| pfm.A | Array | Position Frequency Matrix - Adenine counts |
| pfm.C | Array | Position Frequency Matrix - Cytosine counts |
| pfm.G | Array | Position Frequency Matrix - Guanine counts |
| pfm.T | Array | Position Frequency Matrix - Thymine counts |
| pwm.A | Array | Position Weight Matrix - Adenine log-odds |
| pwm.C | Array | Position Weight Matrix - Cytosine log-odds |
| pwm.G | Array | Position Weight Matrix - Guanine log-odds |
| pwm.T | Array | Position Weight Matrix - Thymine log-odds |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Matrix ID | MA####.# | JASPAR matrix ID with version | MA0106.3 |
| UniProt | Alphanumeric | Protein accession | P04637 |
| NCBI Taxonomy | Numeric | Species ID | 9606 |
| PubMed | Numeric | Publication ID | 12345678 |

---

## Enumerations

### Collections

| Value | Description |
|-------|-------------|
| CORE | Curated high-quality profiles |
| CNE | Conserved non-coding elements |
| PHYLOFACTS | Phylogenetically derived |
| SPLICE | Splice site binding |
| POLII | RNA Polymerase II |
| FAM | Family profiles |

### Taxonomic Groups

| Value | Description |
|-------|-------------|
| vertebrates | Vertebrate species |
| plants | Plant species |
| insects | Insect species |
| fungi | Fungal species |
| nematodes | Nematode species |
| urochordates | Urochordate species |

### TF Structural Classes

| Value | Description |
|-------|-------------|
| Zinc-coordinating | Zinc finger proteins |
| Helix-turn-helix | HTH domain proteins |
| Basic domains | bHLH, bZIP proteins |
| Beta-Scaffold | Beta sheet structures |
| Other | Other structural classes |

### Data Types (Experimental Sources)

| Value | Description |
|-------|-------------|
| ChIP-seq | Chromatin immunoprecipitation sequencing |
| SELEX | Systematic evolution of ligands |
| HT-SELEX | High-throughput SELEX |
| PBM | Protein binding microarray |
| DAP-seq | DNA affinity purification sequencing |
| SMiLE-seq | Selective microfluidics-based SELEX |

---

## Information Content

| Field | Description | Range |
|-------|-------------|-------|
| ic_values | IC per position | 0-2 bits |
| total_ic | Sum of positional IC | 0-2*length |
| max_possible_ic | Maximum possible IC | 2*length |

---

## Entity Relationships

### Matrix to Species
- **Cardinality:** 1:N
- **Description:** One matrix can apply to multiple species
- **Key Fields:** matrix_id, species

### Matrix to Binding Sites
- **Cardinality:** 1:N
- **Description:** Known genomic binding sites
- **Key Fields:** matrix_id, binding_sites

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| JASPAR | Joint Analysis of Specificity of TFs | Database name |
| TF | Transcription Factor | Regulatory protein |
| PFM | Position Frequency Matrix | Raw count matrix |
| PWM | Position Weight Matrix | Log-odds matrix |
| PPM | Position Probability Matrix | Normalized frequency |
| IC | Information Content | Bits per position |
| ChIP-seq | Chromatin Immunoprecipitation Sequencing | Method |
| SELEX | Systematic Evolution of Ligands | Method |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
