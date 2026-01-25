# MSigDB - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | msigdb |
| **Name** | Molecular Signatures Database |
| **Total Fields** | 20 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| geneSetName | String | Yes | Gene set identifier name | HALLMARK_APOPTOSIS |
| systematicName | String | No | MSigDB systematic ID | M5930 |
| geneSetCollection | Enum | Yes | Collection identifier | H, C2, C5 |
| geneSetSubcollection | String | No | Subcollection identifier | CP:KEGG, GO:BP |
| description | String | No | Brief description | Genes involved in apoptosis |
| geneSymbols | Array | No | Gene symbols in the set | ["CASP3", "BAX", "BCL2"] |
| numGenes | Integer | No | Number of genes | 161 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Gene Set Name | UPPERCASE_NAME | Descriptive identifier | HALLMARK_APOPTOSIS |
| Systematic Name | M##### | MSigDB internal ID | M5930 |
| Entrez ID | Numeric | NCBI Gene ID | 836 (CASP3) |
| GEO ID | GSE##### | GEO dataset ID | GSE1234 |
| PubMed | Numeric | Publication ID | 12345678 |

---

## Enumerations

### Gene Set Collections

| Collection | Name | Description |
|------------|------|-------------|
| H | Hallmark | 50 coherent biological states |
| C1 | Positional | Cytogenetic band position |
| C2 | Curated | Curated gene sets (pathways) |
| C3 | Regulatory | Motifs and binding sites |
| C4 | Computational | Computational signatures |
| C5 | Ontology | GO and HPO terms |
| C6 | Oncogenic | Cancer-related signatures |
| C7 | Immunologic | Immunological signatures |
| C8 | Cell Type | Cell type signatures |

### C2 Subcollections

| Subcollection | Description |
|---------------|-------------|
| CP:BIOCARTA | BioCarta pathways |
| CP:KEGG | KEGG pathways |
| CP:REACTOME | Reactome pathways |
| CP:WIKIPATHWAYS | WikiPathways |
| CGP | Chemical/genetic perturbations |

### C5 Subcollections

| Subcollection | Description |
|---------------|-------------|
| GO:BP | GO Biological Process |
| GO:CC | GO Cellular Component |
| GO:MF | GO Molecular Function |
| HPO | Human Phenotype Ontology |

### C3 Subcollections

| Subcollection | Description |
|---------------|-------------|
| TFT:GTRD | GTRD transcription factor targets |
| TFT:TFT_Legacy | Legacy TF targets |
| MIR:MIRDB | miRDB microRNA targets |
| MIR:MIR_Legacy | Legacy miRNA targets |

---

## GSEA Result Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| ES | Number | Enrichment Score |
| NES | Number | Normalized Enrichment Score |
| nominalPValue | Number | Nominal p-value |
| fdrQValue | Number | FDR q-value |
| fwerPValue | Number | Family-wise error rate |
| rankAtMax | Integer | Rank at maximum ES |
| leadingEdge.coreGenes | Array | Core enrichment genes |

---

## Entity Relationships

### Gene Set to Genes
- **Cardinality:** 1:N
- **Description:** Each gene set contains multiple genes
- **Key Fields:** geneSetName, geneSymbols

### Gene Set to Source
- **Cardinality:** N:1
- **Description:** Gene sets derived from publications or databases
- **Key Fields:** geneSetName, pmid, exactSource

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| MSigDB | Molecular Signatures Database | Database name |
| GSEA | Gene Set Enrichment Analysis | Analysis method |
| GMT | Gene Matrix Transposed | File format |
| GRP | Gene set RankedProduct | File format |
| ES | Enrichment Score | GSEA metric |
| NES | Normalized Enrichment Score | Normalized metric |
| FDR | False Discovery Rate | Statistical correction |
| FWER | Family-Wise Error Rate | Statistical correction |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
