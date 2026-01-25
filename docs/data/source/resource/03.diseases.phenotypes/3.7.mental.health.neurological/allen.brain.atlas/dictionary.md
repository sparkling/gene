# Allen Brain Atlas - Data Dictionary

## Overview

This data dictionary documents the schema for Allen Brain Atlas resources.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | allen.brain.atlas |
| **Name** | Allen Brain Atlas |
| **Parent** | 3.7.mental.health.neurological |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Structure Record

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| structure_id | integer | 1:1 | Yes | Structure identifier | 4020 |
| acronym | string | 1:1 | Yes | Structure abbreviation | HIP |
| name | string | 1:1 | Yes | Full structure name | Hippocampus |
| ontology_id | integer | 1:1 | Yes | Ontology reference | 1 |
| parent_id | integer | 1:1 | No | Parent structure | 4008 |
| color_hex_triplet | string | 1:1 | No | Display color | 7ED04B |
| graph_order | integer | 1:1 | No | Tree order | 384 |
| depth | integer | 1:1 | No | Hierarchy depth | 5 |

### Gene Expression

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_symbol | string | 1:1 | Yes | Gene symbol | BDNF |
| gene_id | integer | 1:1 | Yes | Entrez Gene ID | 627 |
| structure_id | integer | 1:1 | Yes | Structure ID | 4020 |
| expression_level | float | 1:1 | Yes | Expression value | 8.45 |
| z_score | float | 1:1 | No | Normalized score | 2.34 |
| probe_id | integer | 1:1 | No | Probe identifier | 1015382 |
| donor_id | string | 1:1 | Yes | Donor identifier | H0351.2001 |

### Donor Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| donor_id | string | 1:1 | Yes | Donor identifier | H0351.2001 |
| name | string | 1:1 | Yes | Donor name | Donor 1 |
| species | string | 1:1 | Yes | Species | Homo sapiens |
| age | integer | 1:1 | No | Age at death | 55 |
| sex | enum | 1:1 | No | Biological sex | Male |
| race | string | 1:1 | No | Race/ethnicity | White |
| hemisphere | string | 1:1 | No | Brain hemisphere | left |

### Single-Cell Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| cell_id | string | 1:1 | Yes | Cell identifier | AAGCGATCAGTC |
| cluster_id | integer | 1:1 | Yes | Cell type cluster | 45 |
| cell_type | string | 1:1 | Yes | Cell type name | L2/3 IT |
| class | string | 1:1 | Yes | Major class | Glutamatergic |
| subclass | string | 1:1 | No | Subclass | IT |
| region | string | 1:1 | Yes | Brain region | VISp |
| umi_counts | integer | 1:1 | No | UMI count | 15000 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Structure ID | Integer | 4020 | Allen structure ID |
| Structure Acronym | Text | HIP | Abbreviation |
| Donor ID | H[0-9]+\.[0-9]+ | H0351.2001 | Human donor |
| Gene Symbol | Text | BDNF | HGNC symbol |
| Entrez ID | Integer | 627 | NCBI Gene ID |
| Probe ID | Integer | 1015382 | Microarray probe |
| Cell ID | Barcode | AAGCGATCAGTC | Single-cell barcode |

---

## Enumerations

### Brain Atlases

| Atlas | Species | Data Types |
|-------|---------|------------|
| Human Brain Atlas | Human | Microarray, ISH |
| BrainSpan | Human | RNA-seq (development) |
| Mouse Brain Atlas | Mouse | ISH, connectivity |
| Cell Types Database | Human/Mouse | Single-cell RNA-seq |
| Allen Mouse Brain Connectivity | Mouse | Tract tracing |

### Expression Data Types

| Type | Description |
|------|-------------|
| Microarray | Gene expression (Agilent) |
| RNA-seq | Bulk RNA sequencing |
| scRNA-seq | Single-cell RNA-seq |
| ISH | In situ hybridization |
| FISH | Fluorescent ISH |

### Cell Type Classes

| Class | Description |
|-------|-------------|
| Glutamatergic | Excitatory neurons |
| GABAergic | Inhibitory neurons |
| Non-Neuronal | Non-neuronal cells |
| Astrocyte | Astrocytes |
| Oligodendrocyte | Oligodendrocytes |
| Microglia | Microglia |
| Endothelial | Endothelial cells |
| Pericyte | Pericytes |

### Brain Regions (Major)

| Acronym | Name | Description |
|---------|------|-------------|
| CTX | Cerebral cortex | Cortical regions |
| HIP | Hippocampus | Memory structure |
| AMY | Amygdala | Emotion processing |
| TH | Thalamus | Relay center |
| HY | Hypothalamus | Homeostasis |
| CB | Cerebellum | Motor coordination |
| BS | Brainstem | Basic functions |
| STR | Striatum | Basal ganglia |

### Donor Demographics

| Category | Values |
|----------|--------|
| Sex | Male, Female |
| Age Range | 18-68 years |
| Hemisphere | Left, Right, Both |
| Species | Homo sapiens, Mus musculus |

### Normalization Methods

| Method | Description |
|--------|-------------|
| Z-score | Structure-wise z-score |
| Log2 | Log2 transformation |
| TPM | Transcripts per million |
| CPM | Counts per million |
| FPKM | Fragments per kilobase million |

---

## Entity Relationships

### Structure Hierarchy
- **Cardinality:** N:1
- **Description:** Brain structure ontology tree
- **Key Fields:** structure_id, parent_id

### Gene to Expression
- **Cardinality:** 1:N
- **Description:** Gene expression across structures
- **Key Fields:** gene_symbol, structure_id

### Donor to Samples
- **Cardinality:** 1:N
- **Description:** Donor provides multiple samples
- **Key Fields:** donor_id, structure_id

### Cell to Cluster
- **Cardinality:** N:1
- **Description:** Cells assigned to cell types
- **Key Fields:** cell_id, cluster_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ISH | In Situ Hybridization | Expression method |
| scRNA-seq | Single-cell RNA Sequencing | Transcriptomics |
| UMI | Unique Molecular Identifier | Count method |
| TPM | Transcripts Per Million | Normalization |
| FPKM | Fragments Per Kilobase Million | Normalization |
| CTX | Cerebral Cortex | Brain region |
| HIP | Hippocampus | Brain region |
| VISp | Primary Visual Cortex | Brain region |
| IT | Intratelencephalic | Neuron type |
| PT | Pyramidal Tract | Neuron type |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| NCBI Gene | Entrez ID | Gene annotation |
| Ensembl | Gene ID | Gene annotation |
| HGNC | Gene symbol | Gene nomenclature |
| Uberon | Anatomy ID | Structure mapping |
| Cell Ontology | CL ID | Cell type mapping |
| PubMed | PMID | Literature |

---

## Data Quality Notes

1. **Human Donors:** 6 donors (2 hemispheres)
2. **Brain Structures:** ~900 anatomically-defined
3. **Genes Profiled:** 20,000+ genes
4. **Single-Cell Profiles:** 1M+ cells
5. **Cell Types:** 100+ transcriptomic types
6. **Mouse Connectivity:** 469 injection experiments
7. **REST API:** Full programmatic access
8. **CC BY:** Free for research and commercial use

