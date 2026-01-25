---
id: download-wikipathways
title: "WikiPathways Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# WikiPathways - Download Documentation

## Overview

WikiPathways provides pathway data via REST API, SPARQL endpoint, and bulk downloads. All data is CC0 (public domain), allowing unrestricted use.

## REST API

### Base URL

```
https://webservice.wikipathways.org
```

### Get Pathway

```bash
# Get pathway info
curl "https://webservice.wikipathways.org/getPathwayInfo?pwId=WP254&format=json"

# Get pathway GPML
curl "https://webservice.wikipathways.org/getPathway?pwId=WP254&format=json"

# Get specific revision
curl "https://webservice.wikipathways.org/getPathway?pwId=WP254&revision=129662&format=json"
```

### Search Pathways

```bash
# Search by text
curl "https://webservice.wikipathways.org/findPathwaysByText?query=apoptosis&species=Homo%20sapiens&format=json"

# Search by gene/protein
curl "https://webservice.wikipathways.org/findPathwaysByXref?ids=ENSG00000141510&codes=En&format=json"

# Search by literature
curl "https://webservice.wikipathways.org/findPathwaysByLiterature?query=12345678&format=json"
```

### List Pathways

```bash
# List all human pathways
curl "https://webservice.wikipathways.org/listPathways?organism=Homo%20sapiens&format=json"

# List all organisms
curl "https://webservice.wikipathways.org/listOrganisms?format=json"

# List curated pathways
curl "https://webservice.wikipathways.org/listPathways?organism=Homo%20sapiens&tag=Curation:AnalysisCollection&format=json"
```

### Get Pathway Components

```bash
# Get genes in pathway
curl "https://webservice.wikipathways.org/getXrefList?pwId=WP254&code=En&format=json"

# Get all cross-references
curl "https://webservice.wikipathways.org/getXrefList?pwId=WP254&format=json"
```

### Export Formats

```bash
# Get as PNG image
curl "https://webservice.wikipathways.org/getColoredPathway?pwId=WP254&graphId=&color=&fileType=png" > pathway.png

# Get as SVG
curl "https://webservice.wikipathways.org/getColoredPathway?pwId=WP254&fileType=svg" > pathway.svg

# Get as PDF
curl "https://webservice.wikipathways.org/getColoredPathway?pwId=WP254&fileType=pdf" > pathway.pdf
```

## Bulk Downloads

### Download Portal

```
https://data.wikipathways.org/
```

### GMT Files (For Enrichment Analysis)

```bash
# Human GMT (gene symbols)
wget https://data.wikipathways.org/current/gmt/wikipathways-20240110-gmt-Homo_sapiens.gmt

# Mouse GMT
wget https://data.wikipathways.org/current/gmt/wikipathways-20240110-gmt-Mus_musculus.gmt
```

### GPML Archives

```bash
# Human GPML files
wget https://data.wikipathways.org/current/gpml/wikipathways-20240110-gpml-Homo_sapiens.zip

# All organisms
wget https://data.wikipathways.org/current/gpml/wikipathways-20240110-gpml-all.zip
```

### RDF/Turtle

```bash
# RDF data dump
wget https://data.wikipathways.org/current/rdf/wikipathways-20240110-rdf-wp.zip

# Void descriptions
wget https://data.wikipathways.org/current/rdf/wikipathways-20240110-rdf-void.ttl
```

### Other Formats

```bash
# BioPAX Level 3
wget https://data.wikipathways.org/current/biopax3/wikipathways-20240110-biopax3-Homo_sapiens.zip

# Gene lists
wget https://data.wikipathways.org/current/gene_lists/wikipathways-20240110-genes-Homo_sapiens.txt
```

## SPARQL Endpoint

### Endpoint URL

```
https://sparql.wikipathways.org/sparql
```

### Example Queries

```sparql
# List all human pathways
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT DISTINCT ?pathway ?title WHERE {
  ?pathway a wp:Pathway ;
           wp:organismName "Homo sapiens" ;
           dcterms:title ?title .
}
ORDER BY ?title
```

```sparql
# Find pathways containing TP53
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT DISTINCT ?pathway ?title ?gene WHERE {
  ?gene a wp:GeneProduct ;
        rdfs:label "TP53"^^xsd:string ;
        dcterms:isPartOf ?pathway .
  ?pathway a wp:Pathway ;
           dcterms:title ?title .
}
```

```sparql
# Get all genes in a pathway
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT DISTINCT ?gene ?label WHERE {
  ?gene a wp:GeneProduct ;
        dcterms:isPartOf <https://identifiers.org/wikipathways/WP254> ;
        rdfs:label ?label .
}
```

## Python Examples

### API Access

```python
import requests

class WikiPathwaysClient:
    def __init__(self):
        self.base_url = "https://webservice.wikipathways.org"

    def get_pathway(self, pathway_id):
        """Get pathway information."""
        response = requests.get(
            f"{self.base_url}/getPathwayInfo",
            params={"pwId": pathway_id, "format": "json"}
        )
        return response.json()

    def search_pathways(self, query, organism="Homo sapiens"):
        """Search pathways by text."""
        response = requests.get(
            f"{self.base_url}/findPathwaysByText",
            params={
                "query": query,
                "species": organism,
                "format": "json"
            }
        )
        return response.json()

    def get_genes(self, pathway_id, code="En"):
        """Get Ensembl gene IDs in pathway."""
        response = requests.get(
            f"{self.base_url}/getXrefList",
            params={
                "pwId": pathway_id,
                "code": code,
                "format": "json"
            }
        )
        return response.json()

    def list_pathways(self, organism="Homo sapiens"):
        """List all pathways for organism."""
        response = requests.get(
            f"{self.base_url}/listPathways",
            params={
                "organism": organism,
                "format": "json"
            }
        )
        return response.json()

# Example usage
client = WikiPathwaysClient()
pathways = client.search_pathways("apoptosis")
genes = client.get_genes("WP254")
```

### Load GMT File

```python
def load_gmt(filepath):
    """Load WikiPathways GMT file."""
    pathways = {}

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            url = parts[1]
            genes = parts[2:]

            # Extract pathway ID from name
            wp_id = name.split('%')[0] if '%' in name else name

            pathways[wp_id] = {
                'name': name,
                'url': url,
                'genes': set(genes)
            }

    return pathways

# Example
wp_human = load_gmt('wikipathways-20240110-gmt-Homo_sapiens.gmt')
print(f"Loaded {len(wp_human)} pathways")
```

### SPARQL Query

```python
from SPARQLWrapper import SPARQLWrapper, JSON

def query_wikipathways(sparql_query):
    """Execute SPARQL query against WikiPathways."""
    endpoint = SPARQLWrapper("https://sparql.wikipathways.org/sparql")
    endpoint.setQuery(sparql_query)
    endpoint.setReturnFormat(JSON)

    results = endpoint.query().convert()
    return results["results"]["bindings"]

# Example: Find pathways with a gene
query = """
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT DISTINCT ?pathway ?title WHERE {
  ?gene a wp:GeneProduct ;
        rdfs:label "EGFR"^^xsd:string ;
        dcterms:isPartOf ?pathway .
  ?pathway a wp:Pathway ;
           dcterms:title ?title ;
           wp:organismName "Homo sapiens" .
}
"""

results = query_wikipathways(query)
for r in results:
    print(f"{r['pathway']['value']}: {r['title']['value']}")
```

## Data Identifier Codes

| Code | Database |
|------|----------|
| En | Ensembl Gene |
| L | Entrez Gene |
| H | HGNC |
| S | UniProt |
| Ce | ChEBI |
| Cs | CAS |
| Wd | Wikidata |
| Re | Reactome |
| Wp | WikiPathways |

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| REST API | No strict limit |
| SPARQL | Timeout at 60 seconds |
| Bulk download | Unlimited |

---

## Dataset Versions

### Current Release: January 2024

| Property | Value |
|----------|-------|
| Version | 20240110 |
| Release Date | 2024-01-10 |
| Total Size | ~500 MB (all formats) |
| Human Pathways | ~800 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| wikipathways-*-gmt-Homo_sapiens.gmt | ~2 MB | ~800 | GSEA gene sets |
| wikipathways-*-gpml-Homo_sapiens.zip | ~50 MB | ~800 | GPML pathway files |
| wikipathways-*-rdf-wp.zip | ~200 MB | All | RDF/Turtle data |
| wikipathways-*-biopax3-Homo_sapiens.zip | ~100 MB | ~800 | BioPAX Level 3 |

### Previous Versions

| Version | Release | Status |
|---------|---------|--------|
| 20231210 | 2023-12-10 | Archived |
| 20231110 | 2023-11-10 | Archived |
| 20231010 | 2023-10-10 | Archived |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| REST Base URL | `https://webservice.wikipathways.org` |
| SPARQL Endpoint | `https://sparql.wikipathways.org/sparql` |
| Authentication | None required |
| Rate Limit | No strict limit (60s SPARQL timeout) |
| Response Format | JSON, XML, PNG, SVG, PDF |

### REST API Endpoints

| Operation | Endpoint | Example |
|-----------|----------|---------|
| Get Pathway | `/getPathway` | `?pwId=WP254&format=json` |
| Search | `/findPathwaysByText` | `?query=apoptosis&species=Homo%20sapiens` |
| List Pathways | `/listPathways` | `?organism=Homo%20sapiens&format=json` |
| Get Genes | `/getXrefList` | `?pwId=WP254&code=En&format=json` |
| Export Image | `/getColoredPathway` | `?pwId=WP254&fileType=png` |

---

## Update Frequency

| Type | Frequency |
|------|-----------|
| Monthly release | First of month |
| Continuous edits | Real-time (API) |
| RDF dump | Monthly |

## See Also

- [Schema Documentation](./schema.md)
- [WikiPathways API Documentation](https://www.wikipathways.org/index.php/Help:WikiPathways_Webservice)
