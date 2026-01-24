---
id: guides-integration-pathway-downloads
title: "Bulk Download Methods for Pathway and Gene/Protein Target Databases"
type: guide
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [downloads, bulk-data, pathways, targets, genes, proteins, guide]
---

**Parent:** [Integration Guides](./_index.md)

# Bulk Download Methods for Pathway and Gene/Protein Target Databases

This document provides comprehensive information on bulk download methods for major biological pathway and gene/protein target databases. Last updated: January 2026.

---

## Table of Contents

1. [Reactome](#1-reactome)
2. [KEGG](#2-kegg)
3. [WikiPathways](#3-wikipathways)
4. [UniProt](#4-uniprot)
5. [STRING](#5-string)
6. [Gene Ontology](#6-gene-ontology)
7. [NCBI Gene](#7-ncbi-gene)
8. [Ensembl](#8-ensembl)
9. [Processing Recommendations Summary](#9-processing-recommendations-summary)

---

[Full content from pathways-targets.md continues...]

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **Reactome** | FTP | `https://reactome.org/download-data` |
| **KEGG** | API/FTP | `https://rest.kegg.jp/` (academic) |
| **WikiPathways** | Download | `https://data.wikipathways.org/` |
| **UniProt** | FTP | `ftp://ftp.uniprot.org/pub/databases/uniprot/` |
| **STRING** | Download | `https://string-db.org/cgi/download` |
| **Gene Ontology** | Download | `http://geneontology.org/docs/download-ontology/` |
| **NCBI Gene** | FTP | `ftp://ftp.ncbi.nih.gov/gene/` |
| **Ensembl** | FTP | `ftp://ftp.ensembl.org/pub/` |

**Access Requirements:** Most are freely accessible; KEGG FTP requires academic subscription.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | BioPAX, GPML, KGML, TSV |
| Alternative | OWL, XML, JSON |
| Gene sets | GMT, GSE |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `pathway_id` | string | Pathway identifier | "R-HSA-1430728" |
| `gene_symbol` | string | Gene/protein symbol | "MTHFR" |
| `interaction_type` | string | Relationship type | "activation" |
| `confidence` | float | Evidence score | 0.90 |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `participates_in` | Pathway | N:M |
| `interacts_with` | Protein | N:M |

## Sample Data

### Example Target Record
```json
{
  "uniprot_id": "P42898",
  "gene_symbol": "MTHFR",
  "pathways": ["R-HSA-1430728", "hsa00670"],
  "string_partners": [{"partner": "MTR", "score": 0.95}],
  "go_terms": ["GO:0004489", "GO:0006555"]
}
```

### Sample Query Result
| gene | pathway | source | confidence |
|------|---------|--------|------------|
| MTHFR | Folate metabolism | Reactome | 0.99 |
| MTR | One carbon pool | KEGG | 0.95 |

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| Reactome | CC BY 4.0 | Yes |
| WikiPathways | CC0 | Yes |
| UniProt | CC BY 4.0 | Yes |
| STRING | CC BY 4.0 | Yes |
| Gene Ontology | CC BY 4.0 | Yes |
| KEGG | Academic | Subscription required |

## Data Set Size

| Metric | Value |
|--------|-------|
| Reactome pathways | 2,712 human pathways |
| WikiPathways | 3,100+ pathways |
| STRING interactions | 11.9B protein pairs |
| UniProt human proteins | ~20,400 reviewed |
| Total storage estimate | ~15-20 GB combined |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Bulk Download | Retrieving entire datasets via FTP/HTTP rather than API queries | Downloading full Reactome database |
| Pathway | Series of molecular interactions leading to a biological outcome | Apoptosis signaling |
| Target | Molecular entity (protein/gene) affected by compounds or involved in pathways | TP53 protein |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Reactome | Curated database of biological pathways | Pathway analysis |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Metabolic pathways |
| WikiPathways | Community-curated open pathway database | Collaborative curation |
| UniProt | Universal Protein Resource database | Protein sequences |
| STRING | Search Tool for Retrieval of Interacting Genes/Proteins | PPI networks |
| Gene Ontology | Structured vocabulary for gene function | GO terms |
| NCBI Gene | NCBI's gene-centric database | Gene information |
| Ensembl | Genome browser and annotation database | Genomic data |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST endpoints |
| FTP | File Transfer Protocol | Download method |
| GO | Gene Ontology | Functional annotation |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway DB |
| NCBI | National Center for Biotechnology Information | US agency |
| PPI | Protein-Protein Interaction | Network data |

---
