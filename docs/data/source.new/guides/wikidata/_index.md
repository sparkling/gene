---
id: guides-wikidata
title: "Wikidata Guides"
type: directory
parent: ../_index.md
description: "SPARQL queries and integration guides for Wikidata biomedical data"
last_updated: 2026-01-23
status: active
children:
  - master-reference.md
  - schema.md
  - pharmaceutical-queries.md
  - traditional-medicine.md
  - supplements.md
  - bulk-download.md
  - wikipedia-integration.md
tags: [wikidata, sparql, knowledge-graph, cc0, guides]
---

**Parent:** [Guides](../_index.md)

# Wikidata Guides

Comprehensive guides for extracting biomedical data from Wikidata using SPARQL queries.

## Overview

Wikidata is a free, structured knowledge base containing:
- 59,000+ human genes with Entrez/Ensembl IDs
- 27,000+ proteins with UniProt mappings
- 200,000+ diseases with OMIM/MONDO/ICD codes
- 45,000+ medications with DrugBank/ChEMBL links
- 10,000+ medicinal plants with compound data

**License:** CC0 (Public Domain) - No attribution required

## Guides in This Section

| Guide | Description |
|-------|-------------|
| [master-reference.md](./master-reference.md) | Complete Q-ID and P-code reference for biomedical extraction |
| [schema.md](./schema.md) | Wikidata biomedical schema and property documentation |
| [pharmaceutical-queries.md](./pharmaceutical-queries.md) | Drug, target, and pharmacological data queries |
| [traditional-medicine.md](./traditional-medicine.md) | TCM, Ayurveda, Kampo, and medicinal plant queries |
| [supplements.md](./supplements.md) | Dietary supplements and nutraceutical queries |
| [bulk-download.md](./bulk-download.md) | Wikidata, Wikipedia, and DBpedia bulk download methods |
| [wikipedia-integration.md](./wikipedia-integration.md) | Wikipedia API and Wikidata cross-linking |

## Quick Start

### SPARQL Endpoint

```
Query Interface: https://query.wikidata.org/
SPARQL Endpoint: https://query.wikidata.org/sparql
```

### Example Query: Human Genes with Identifiers

```sparql
SELECT ?gene ?geneLabel ?entrezId ?ensembl ?uniprot WHERE {
  ?gene wdt:P31 wd:Q7187 .           # instance of: gene
  ?gene wdt:P703 wd:Q15978631 .      # found in taxon: Homo sapiens
  ?gene wdt:P351 ?entrezId .         # Entrez Gene ID
  OPTIONAL { ?gene wdt:P594 ?ensembl . }
  OPTIONAL { ?gene wdt:P352 ?uniprot . }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

## Key Properties Reference

### Gene Properties
| Property | Name | Example |
|----------|------|---------|
| P351 | Entrez Gene ID | "7157" |
| P352 | UniProt ID | "P04637" |
| P353 | HGNC Symbol | "TP53" |
| P594 | Ensembl Gene ID | "ENSG00000141510" |

### Disease Properties
| Property | Name | Example |
|----------|------|---------|
| P699 | Disease Ontology ID | "DOID:162" |
| P492 | OMIM ID | "191170" |
| P5270 | MONDO ID | "MONDO:0005148" |

### Compound Properties
| Property | Name | Example |
|----------|------|---------|
| P662 | PubChem CID | "2244" |
| P592 | ChEMBL ID | "CHEMBL25" |
| P715 | DrugBank ID | "DB00945" |
| P231 | CAS Number | "50-78-2" |

## Related Resources

- [Wikidata Resource Schema](../../resource/08.literature.knowledge/8.2.knowledge.bases/wikidata/schema.md)
- [Integration Guides](../integration/_index.md)

---

*Last Updated: January 2026*
