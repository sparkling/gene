---
id: schemas-string-schema
title: "STRING Database Schema"
category: schemas
parent: README.md
last_updated: 2026-01-22
status: migrated
tags: [schema, string, protein-interaction, ppi, network, functional-enrichment]
---

**Parent:** [Schema Documentation](./README.md)

# STRING Database Schema

**Document ID:** STRING-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** STRING v12.0

---

## TL;DR

STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) is a comprehensive database of protein-protein interactions covering 14,094 organisms with 59.3 million proteins and 12.5 billion interactions. Interactions are scored from multiple evidence channels (experimental, text-mining, genomic context, co-expression) and combined into a single confidence score (0-1000). The REST API provides network visualization, functional enrichment, and interaction partner queries.

---

## Database Statistics (v12.0)

| Metric | Count |
|--------|-------|
| **Organisms** | 14,094 |
| **Proteins** | 59,309,604 |
| **Interactions** | 12,535,845,684 |
| **Human Proteins** | ~20,000 |
| **Human Interactions** | ~11.8 million |

---

## License

**License:** Creative Commons Attribution 4.0 (CC BY 4.0)
**Attribution Required:** Yes
**Commercial Use:** Allowed
**URL:** https://string-db.org/cgi/access

---

## API Overview

### Base URL
```
https://string-db.org/api/
```

### Response Formats
- `json` - JSON format
- `tsv` - Tab-separated values
- `tsv-no-header` - TSV without header row
- `xml` - XML format
- `psi-mi` - PSI-MI XML format
- `psi-mi-tab` - PSI-MITAB format

### Request Format
```
https://string-db.org/api/{format}/{method}?{parameters}
```

---

## Core API Endpoints

### 1. Network Retrieval

#### GET /network
Returns the full interaction network for given proteins.

**Parameters:**

| Parameter | Required | Type | Description |
|-----------|----------|------|-------------|
| `identifiers` | Yes | String | Protein identifiers (newline or %0d separated) |
| `species` | No | Integer | NCBI taxonomy ID (e.g., 9606 for human) |
| `required_score` | No | Integer | Minimum combined score (0-1000, default: 400) |
| `network_type` | No | String | `functional` or `physical` (default: functional) |
| `add_nodes` | No | Integer | Number of additional interactors to add (0-100) |
| `show_query_node_labels` | No | Integer | Show labels for query proteins (0/1) |

**Sample Request:**
```
https://string-db.org/api/json/network?identifiers=TP53%0dBRCA1&species=9606&required_score=400
```

**Sample Response:**
```json
[
  {
    "stringId_A": "9606.ENSP00000269305",
    "stringId_B": "9606.ENSP00000418960",
    "preferredName_A": "TP53",
    "preferredName_B": "BRCA1",
    "ncbiTaxonId": 9606,
    "score": 0.999,
    "nscore": 0.000,
    "fscore": 0.900,
    "pscore": 0.000,
    "ascore": 0.000,
    "escore": 0.994,
    "dscore": 0.900,
    "tscore": 0.969
  }
]
```

---

#### GET /interaction_partners
Returns interaction partners for given proteins.

**Parameters:**

| Parameter | Required | Type | Description |
|-----------|----------|------|-------------|
| `identifiers` | Yes | String | Query protein identifier(s) |
| `species` | No | Integer | NCBI taxonomy ID |
| `limit` | No | Integer | Max partners per protein (default: 10) |
| `required_score` | No | Integer | Minimum score threshold (0-1000) |
| `network_type` | No | String | `functional` or `physical` |

**Sample Request:**
```
https://string-db.org/api/json/interaction_partners?identifiers=TP53&species=9606&limit=5
```

**Sample Response:**
```json
[
  {
    "stringId_A": "9606.ENSP00000269305",
    "stringId_B": "9606.ENSP00000357274",
    "preferredName_A": "TP53",
    "preferredName_B": "MDM2",
    "ncbiTaxonId": 9606,
    "score": 0.999,
    "nscore": 0.000,
    "fscore": 0.999,
    "pscore": 0.000,
    "ascore": 0.000,
    "escore": 0.998,
    "dscore": 0.900,
    "tscore": 0.980
  }
]
```

---

### 2. Functional Enrichment

#### GET /enrichment
Performs functional enrichment analysis on a protein set.

**Parameters:**

| Parameter | Required | Type | Description |
|-----------|----------|------|-------------|
| `identifiers` | Yes | String | Protein identifiers |
| `species` | No | Integer | NCBI taxonomy ID |
| `caller_identity` | No | String | Your application ID |

**Sample Request:**
```
https://string-db.org/api/json/enrichment?identifiers=TP53%0dBRCA1%0dBRCA2%0dATM&species=9606
```

**Sample Response:**
```json
[
  {
    "category": "Process",
    "term": "GO:0006281",
    "description": "DNA repair",
    "number_of_genes": 4,
    "number_of_genes_in_background": 484,
    "p_value": 1.52e-8,
    "fdr": 3.21e-6,
    "preferredNames": "TP53,BRCA1,BRCA2,ATM"
  },
  {
    "category": "KEGG",
    "term": "hsa03440",
    "description": "Homologous recombination",
    "number_of_genes": 3,
    "number_of_genes_in_background": 41,
    "p_value": 4.67e-7,
    "fdr": 2.11e-5,
    "preferredNames": "BRCA1,BRCA2,ATM"
  }
]
```

**Enrichment Categories:**
- `Process` - GO Biological Process
- `Function` - GO Molecular Function
- `Component` - GO Cellular Component
- `KEGG` - KEGG Pathways
- `Pfam` - Protein families
- `InterPro` - Protein domains
- `SMART` - Domain annotations
- `Reactome` - Reactome pathways
- `WikiPathways` - WikiPathways
- `HPO` - Human Phenotype Ontology
- `DisGeNET` - Disease-gene associations
- `NetworkNeighborAL` - Network-based annotations

---

### 3. Protein Mapping

#### GET /get_string_ids
Maps identifiers to STRING internal IDs.

**Parameters:**

| Parameter | Required | Type | Description |
|-----------|----------|------|-------------|
| `identifiers` | Yes | String | Query identifiers |
| `species` | No | Integer | NCBI taxonomy ID |
| `limit` | No | Integer | Max matches per identifier |
| `echo_query` | No | Integer | Include query in response (0/1) |
| `caller_identity` | No | String | Application identifier |

**Sample Request:**
```
https://string-db.org/api/json/get_string_ids?identifiers=p53&species=9606&limit=1
```

**Sample Response:**
```json
[
  {
    "queryIndex": 0,
    "queryItem": "p53",
    "stringId": "9606.ENSP00000269305",
    "ncbiTaxonId": 9606,
    "taxonName": "Homo sapiens",
    "preferredName": "TP53",
    "annotation": "Tumor protein p53"
  }
]
```

---

### 4. Visualization

#### GET /network_image
Returns a PNG image of the interaction network.

**Parameters:**

| Parameter | Required | Type | Description |
|-----------|----------|------|-------------|
| `identifiers` | Yes | String | Protein identifiers |
| `species` | No | Integer | NCBI taxonomy ID |
| `required_score` | No | Integer | Minimum score (0-1000) |
| `network_flavor` | No | String | `evidence`, `confidence`, `actions` |
| `add_color_nodes` | No | Integer | Color nodes (0/1) |
| `add_white_nodes` | No | Integer | Add white interactors |

**Sample Request:**
```
https://string-db.org/api/image/network?identifiers=TP53%0dMDM2%0dCDKN1A&species=9606
```

---

### 5. Homology Mapping

#### GET /homology
Maps proteins to homologs in another species.

**Parameters:**

| Parameter | Required | Type | Description |
|-----------|----------|------|-------------|
| `identifiers` | Yes | String | Source protein identifiers |
| `species` | Yes | Integer | Source species taxonomy ID |
| `species_b` | Yes | Integer | Target species taxonomy ID |

**Sample Response:**
```json
[
  {
    "stringId_A": "9606.ENSP00000269305",
    "stringId_B": "10090.ENSMUSP00000104298",
    "ncbiTaxonId_A": 9606,
    "ncbiTaxonId_B": 10090,
    "preferredName_A": "TP53",
    "preferredName_B": "Trp53",
    "bitscore": 356.0
  }
]
```

---

## Score Channels (Field Dictionary)

STRING combines evidence from multiple sources into channel-specific scores and a combined score.

### Score Fields

| Field | Type | Range | Description |
|-------|------|-------|-------------|
| `score` | Float | 0-1 | Combined confidence score |
| `nscore` | Float | 0-1 | **Neighborhood score** - Gene neighborhood in genome |
| `fscore` | Float | 0-1 | **Fusion score** - Gene fusion events |
| `pscore` | Float | 0-1 | **Phylogenetic score** - Co-occurrence across genomes |
| `ascore` | Float | 0-1 | **Co-expression score** - Co-expression patterns |
| `escore` | Float | 0-1 | **Experimental score** - Laboratory experiments |
| `dscore` | Float | 0-1 | **Database score** - Curated pathway databases |
| `tscore` | Float | 0-1 | **Text-mining score** - Literature co-mentions |

### Score Computation

The combined score uses a probabilistic integration:
```
combined = 1 - ∏(1 - score_i)
```
where score_i is each individual channel score after prior correction.

### Score Confidence Levels

| Combined Score | Confidence | Interpretation |
|---------------|------------|----------------|
| 0.900 - 1.000 | Highest | Very high confidence |
| 0.700 - 0.899 | High | High confidence |
| 0.400 - 0.699 | Medium | Medium confidence |
| 0.150 - 0.399 | Low | Low confidence |
| 0.000 - 0.149 | Lowest | Very low confidence |

---

## Interaction Data Schema

### Network Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `stringId_A` | String | STRING ID for protein A (taxid.ENSP#####) |
| `stringId_B` | String | STRING ID for protein B |
| `preferredName_A` | String | Gene symbol for protein A |
| `preferredName_B` | String | Gene symbol for protein B |
| `ncbiTaxonId` | Integer | NCBI taxonomy ID |
| `score` | Float | Combined confidence score |
| `nscore` | Float | Neighborhood score |
| `fscore` | Float | Fusion score |
| `pscore` | Float | Phylogenetic profile score |
| `ascore` | Float | Co-expression score |
| `escore` | Float | Experimental score |
| `dscore` | Float | Database score |
| `tscore` | Float | Text-mining score |

### Enrichment Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `category` | String | Annotation category (GO, KEGG, etc.) |
| `term` | String | Term identifier |
| `description` | String | Term description |
| `number_of_genes` | Integer | Genes in input matching term |
| `number_of_genes_in_background` | Integer | Total genes with term |
| `p_value` | Float | Hypergeometric p-value |
| `fdr` | Float | False discovery rate (Benjamini-Hochberg) |
| `preferredNames` | String | Comma-separated gene names |

### Protein Mapping Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `queryIndex` | Integer | Index of query identifier |
| `queryItem` | String | Original query string |
| `stringId` | String | STRING protein identifier |
| `ncbiTaxonId` | Integer | NCBI taxonomy ID |
| `taxonName` | String | Species name |
| `preferredName` | String | Preferred gene/protein name |
| `annotation` | String | Protein description |

---

## Bulk Data Downloads

### Available Files

| File | Description | Format |
|------|-------------|--------|
| `protein.links.v12.0.txt.gz` | All interactions | TSV |
| `protein.links.detailed.v12.0.txt.gz` | Interactions with channel scores | TSV |
| `protein.links.full.v12.0.txt.gz` | Full interaction data | TSV |
| `protein.info.v12.0.txt.gz` | Protein annotations | TSV |
| `protein.aliases.v12.0.txt.gz` | Protein name aliases | TSV |
| `protein.sequences.v12.0.fa.gz` | Protein sequences | FASTA |
| `species.v12.0.txt` | Species list | TSV |

**Download URL:** https://string-db.org/cgi/download

### Protein Links File Schema

```
protein1    protein2    combined_score
9606.ENSP00000269305    9606.ENSP00000357274    999
9606.ENSP00000269305    9606.ENSP00000418960    999
```

### Detailed Links File Schema

```
protein1    protein2    neighborhood    fusion    cooccurence    coexpression    experimental    database    textmining    combined_score
9606.ENSP00000269305    9606.ENSP00000357274    0    0    0    0    998    999    980    999
```

---

## Cross-References

### Supported Identifier Types

| Database | ID Format | Example |
|----------|-----------|---------|
| STRING Internal | taxid.ENSP##### | 9606.ENSP00000269305 |
| Ensembl Protein | ENSP########### | ENSP00000269305 |
| Ensembl Gene | ENSG########### | ENSG00000141510 |
| UniProt | [A-Z][0-9]{5} | P04637 |
| Gene Symbol | Text | TP53 |
| RefSeq | NP_###### | NP_000537 |
| NCBI Gene ID | Integer | 7157 |

### Species Taxonomy IDs (Common)

| Species | NCBI Tax ID |
|---------|-------------|
| Homo sapiens | 9606 |
| Mus musculus | 10090 |
| Rattus norvegicus | 10116 |
| Drosophila melanogaster | 7227 |
| Caenorhabditis elegans | 6239 |
| Saccharomyces cerevisiae | 4932 |
| Escherichia coli K12 | 83333 |
| Arabidopsis thaliana | 3702 |
| Danio rerio | 7955 |

---

## Sample Data

### Sample Network Query (Human TP53-BRCA1)

**Request:**
```bash
curl "https://string-db.org/api/json/network?identifiers=TP53%0dBRCA1&species=9606"
```

**Response:**
```json
[
  {
    "stringId_A": "9606.ENSP00000269305",
    "stringId_B": "9606.ENSP00000418960",
    "preferredName_A": "TP53",
    "preferredName_B": "BRCA1",
    "ncbiTaxonId": 9606,
    "score": 0.999,
    "nscore": 0.000,
    "fscore": 0.000,
    "pscore": 0.000,
    "ascore": 0.076,
    "escore": 0.994,
    "dscore": 0.900,
    "tscore": 0.969
  }
]
```

### Sample Enrichment Query

**Request:**
```bash
curl "https://string-db.org/api/json/enrichment?identifiers=TP53%0dBRCA1%0dBRCA2&species=9606"
```

**Response (abbreviated):**
```json
[
  {
    "category": "Process",
    "term": "GO:0006281",
    "description": "DNA repair",
    "number_of_genes": 3,
    "number_of_genes_in_background": 484,
    "p_value": 2.45e-6,
    "fdr": 1.82e-4,
    "preferredNames": "TP53,BRCA1,BRCA2"
  },
  {
    "category": "Process",
    "term": "GO:0006974",
    "description": "cellular response to DNA damage stimulus",
    "number_of_genes": 3,
    "number_of_genes_in_background": 594,
    "p_value": 4.62e-6,
    "fdr": 2.31e-4,
    "preferredNames": "TP53,BRCA1,BRCA2"
  }
]
```

---

## Integration Notes

### Rate Limits

- API calls are limited; add `caller_identity` parameter for tracking
- For bulk queries, use downloaded flat files
- Batch identifiers (up to 2000) in single request using newline separators

### Data Quality

- **Experimental scores (escore)** are most reliable for direct interactions
- **Text-mining scores (tscore)** may include false positives from sentence co-occurrence
- **Database scores (dscore)** from curated sources like Reactome and KEGG
- Use `required_score >= 700` for high-confidence interactions

### Network Types

- **Functional**: All association types (default)
- **Physical**: Only physical binding interactions

---

## Data Set Size

| Metric | Value |
|--------|-------|
| **Total Organisms** | 14,094 |
| **Total Proteins** | 59,309,604 |
| **Total Interactions** | 12,535,845,684 |
| **Human Proteins** | ~20,000 |
| **Human Interactions** | ~11.8 million |
| **Human Protein-Protein Links** | 500,000+ high-confidence (score >= 700) |
| **Interaction Pairs** | 68,000,000+ (all evidence types) |
| **Average Proteins per Organism** | ~4,200 |
| **Average Interactions per Organism** | ~890,000 |
| **Text-Mining Evidence** | 1+ billion co-mentions |
| **Experimental Evidence** | 5+ million interactions |
| **Curated Database Evidence** | 2+ million interactions |
| **Database Downloads** | ~500 GB total (all formats) |
| **Bulk File Size** | protein.links.v12.0.txt.gz: ~3 GB |
| **Detailed Links File** | protein.links.detailed.v12.0.txt.gz: ~8 GB |
| **Update Frequency** | Approximately annual major releases |
| **Last Version** | v12.0 (released 2023) |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 12,500,000,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | TSV (downloadable), JSON (API) |
| Alternative | XML |
| Encoding | UTF-8 |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "9606.ENSP00000000003" |
| `name` | string | Entity name | "Protein A" |
| `type` | string | Record type | "protein" / "interaction" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `interacts_with` | Protein | N:M |
| `has_score` | Evidence | N:M |
| `cross_references` | Database | N:M |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| STRING | HTTP | https://string-db.org/download |
| Protein Links | FTP | https://string-db.org/download/protein.links.full.v12.0.txt.gz |
| Protein Sequences | FTP | https://string-db.org/download/protein.sequences.v12.0.fa.gz |

**Access Requirements:** Open access, no registration required

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| STRING | CC BY 4.0 | Yes |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `stringId` | STRING protein identifier (taxid.ENSP#####) | 9606.ENSP00000269305 |
| `preferredName` | Gene symbol or preferred protein name | TP53 |
| `score` | Combined confidence score (0-1 or 0-1000) | 0.999 or 999 |
| `nscore` | Neighborhood score - gene proximity in genome | 0.000 |
| `fscore` | Fusion score - gene fusion events | 0.000 |
| `pscore` | Phylogenetic profile score - co-occurrence | 0.000 |
| `ascore` | Co-expression score - expression correlation | 0.076 |
| `escore` | Experimental score - laboratory evidence | 0.994 |
| `dscore` | Database score - curated pathway sources | 0.900 |
| `tscore` | Text-mining score - literature co-mentions | 0.969 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Protein-Protein Interaction | Physical or functional association between proteins | Core data model |
| Functional Network | Network of proteins connected by any evidence type | network_type: functional |
| Physical Network | Network limited to physical binding interactions | network_type: physical |
| Evidence Channel | Source type for interaction (experimental, text-mining, etc.) | Score components |
| Confidence Score | Probabilistic measure of interaction reliability | Combined score |
| Functional Enrichment | Statistical over-representation of GO/KEGG terms | /enrichment endpoint |
| Interaction Partner | Protein directly associated with query protein | /interaction_partners |
| Homology Mapping | Finding orthologous proteins across species | /homology endpoint |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| STRING | Search Tool for Retrieval of Interacting Genes/Proteins | Database name |
| PPI | Protein-Protein Interaction | Core data type |
| NCBI | National Center for Biotechnology Information | Taxonomy source |
| ENSP | Ensembl Protein | Identifier format |
| GO | Gene Ontology | Enrichment category |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway enrichment |
| FDR | False Discovery Rate | Enrichment statistics |
| CC BY | Creative Commons Attribution | License type |
| TSV | Tab-Separated Values | Download format |
| PSI-MI | Proteomics Standards Initiative Molecular Interaction | Interaction standard |
| BioGRID | Biological General Repository for Interaction Datasets | Related database |
| HPO | Human Phenotype Ontology | Enrichment category |

---

## References

1. Szklarczyk D, et al. (2023) "The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest." Nucleic Acids Res. 51(D1):D D99-D707.

2. STRING API Documentation: https://string-db.org/cgi/help?subpage=api

3. STRING Downloads: https://string-db.org/cgi/download

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation for STRING v12.0 |
