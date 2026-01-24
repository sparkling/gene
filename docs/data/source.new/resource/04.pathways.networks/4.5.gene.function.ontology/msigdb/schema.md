---
id: schema-msigdb
title: "MSigDB Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: final
tags: [schema, database, gene-sets, enrichment]
---

# MSigDB - Schema Documentation

## TL;DR

MSigDB provides gene sets in multiple formats optimized for enrichment analysis. The primary format is GMT (Gene Matrix Transposed), with additional support for XML, JSON, and GRP formats.

## GMT Format (Primary)

### Structure

Tab-delimited file with variable columns:

```
GENE_SET_NAME<TAB>DESCRIPTION<TAB>GENE1<TAB>GENE2<TAB>GENE3<TAB>...
```

### Example

```
HALLMARK_APOPTOSIS	http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_APOPTOSIS	CASP3	CASP8	BAX	BCL2	FADD	BID	CYCS	APAF1
HALLMARK_HYPOXIA	http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_HYPOXIA	HIF1A	VEGFA	LDHA	PGK1	ENO1	PKM	SLC2A1
```

### Columns

| Position | Content |
|----------|---------|
| 1 | Gene set name |
| 2 | Description/URL |
| 3+ | Gene symbols (one per column) |

## Gene Set Object (JSON/API)

```json
{
  "geneSetName": "HALLMARK_APOPTOSIS",
  "systematicName": "M5930",
  "geneSetCollection": "H",
  "geneSetSubcollection": null,
  "description": "Genes mediating programmed cell death (apoptosis) by activation of caspases.",
  "contributor": "Arthur Liberzon",
  "exactSource": "http://www.broadinstitute.org/gsea/msigdb/",
  "externalDetailsUrl": "http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_APOPTOSIS",
  "pmid": null,
  "geoid": null,
  "geneSymbols": ["CASP3", "CASP8", "BAX", "BCL2", ...],
  "numGenes": 161,
  "organisms": ["Homo sapiens"],
  "chipPlatforms": null,
  "tags": ["hallmark", "apoptosis"]
}
```

## Gene Set Collections

### Collection Identifiers

| ID | Name | Description | Gene Sets |
|----|------|-------------|-----------|
| H | Hallmark | Well-defined biological states | 50 |
| C1 | Positional | Chromosome cytogenetic bands | 299 |
| C2 | Curated | Curated pathway databases | 6,300+ |
| C3 | Motif | TF and miRNA targets | 3,700+ |
| C4 | Computational | Cancer-related signatures | 850+ |
| C5 | GO | Gene Ontology terms | 15,000+ |
| C6 | Oncogenic | Oncogenic pathway signatures | 189 |
| C7 | Immunologic | Immunologic signatures | 5,200+ |
| C8 | Cell Type | Cell type markers | 700+ |

### Subcollections

```
C2:
  CP:BIOCARTA  - BioCarta pathways
  CP:KEGG      - KEGG pathways
  CP:PID       - NCI-PID pathways
  CP:REACTOME  - Reactome pathways
  CP:WIKIPATHWAYS - WikiPathways
  CGP          - Chemical/Genetic Perturbations

C3:
  TFT:GTRD     - GTRD TF targets
  TFT:TFT_Legacy - Legacy TF targets
  MIR:MIRDB    - miRNA targets (mirDB)
  MIR:MIR_Legacy - Legacy miRNA targets

C5:
  GO:BP        - Biological Process
  GO:MF        - Molecular Function
  GO:CC        - Cellular Component
  HPO          - Human Phenotype Ontology
```

## Gene Set Card (Web)

```json
{
  "card": {
    "geneSetName": "HALLMARK_APOPTOSIS",
    "systematicName": "M5930",
    "collection": {
      "name": "H: hallmark gene sets",
      "id": "H"
    },
    "description": "Genes mediating programmed cell death...",
    "fullDescription": "Detailed description from original source...",
    "sourcePublication": {
      "pmid": null,
      "title": null
    },
    "geneSetSize": 161,
    "msigdbUrl": "http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_APOPTOSIS",
    "genes": [
      {"symbol": "CASP3", "title": "caspase 3", "entrezId": 836},
      {"symbol": "CASP8", "title": "caspase 8", "entrezId": 841}
    ],
    "overlappingGeneSets": [...],
    "enrichmentAnalysis": {...}
  }
}
```

## GMX Format

Alternative column-oriented format:

```
GeneSet1	GeneSet2	GeneSet3
Description1	Description2	Description3
GENE1A	GENE2A	GENE3A
GENE1B	GENE2B	GENE3B
GENE1C	GENE2C
```

## GRP Format (Single Gene Set)

Simple list of genes:

```
# HALLMARK_APOPTOSIS
CASP3
CASP8
BAX
BCL2
FADD
BID
```

## XML Format

```xml
<?xml version="1.0" encoding="UTF-8"?>
<MSIGDB VERSION="2024.1" SPECIES="Homo sapiens">
  <GENESET NAME="HALLMARK_APOPTOSIS"
           SYSTEMATIC_NAME="M5930"
           STANDARD_NAME="HALLMARK_APOPTOSIS"
           COLLECTION="H"
           ORGANISM="Homo sapiens"
           DESCRIPTION_BRIEF="Genes mediating programmed cell death..."
           DESCRIPTION_FULL="..."
           MEMBERS="CASP3,CASP8,BAX,BCL2,..."
           MEMBERS_EZID="836,841,581,596,..."
           SIZE="161"
           CONTRIBUTOR="Arthur Liberzon"
           EXTERNAL_URL="..."
           PMID=""/>
</MSIGDB>
```

## Gene Identifier Formats

| Version | Identifier Type | Example File |
|---------|-----------------|--------------|
| Symbols | HGNC symbols | h.all.v2024.1.Hs.symbols.gmt |
| Entrez | NCBI Gene IDs | h.all.v2024.1.Hs.entrez.gmt |

## GSEA Enrichment Results

```
NAME	GS<br>follow link	GS DETAILS	SIZE	ES	NES	NOM p-val	FDR q-val	FWER p-val	RANK AT MAX	LEADING EDGE
HALLMARK_APOPTOSIS	HALLMARK_APOPTOSIS	Details...	161	0.42	1.85	0.000	0.001	0.002	4532	tags=35%, list=18%, signal=43%
```

### Result Fields

| Field | Description |
|-------|-------------|
| ES | Enrichment Score |
| NES | Normalized Enrichment Score |
| NOM p-val | Nominal p-value |
| FDR q-val | False Discovery Rate |
| FWER p-val | Family-Wise Error Rate |
| RANK AT MAX | Rank at max ES |
| LEADING EDGE | Core enrichment genes |

## Leading Edge Analysis

```json
{
  "leadingEdge": {
    "geneSet": "HALLMARK_APOPTOSIS",
    "tags": 0.35,
    "list": 0.18,
    "signal": 0.43,
    "coreGenes": ["CASP3", "BAX", "BID", "FADD", "CYCS"]
  }
}
```

- **tags**: Percentage of genes contributing to ES
- **list**: Percentage of ranked list before ES peak
- **signal**: tags * (1 - list)

## See Also

- [Download Documentation](./download.md)
- [MSigDB Collections](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
