---
id: download-gene-ontology
title: "Gene Ontology Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# Gene Ontology - Download Documentation

## Overview

Gene Ontology provides the ontology files and gene annotations via REST API, SPARQL endpoint, and bulk downloads. All data is available under CC BY 4.0.

## Download Portal

### Primary URLs

```
https://current.geneontology.org/
http://geneontology.org/docs/download-ontology/
http://geneontology.org/docs/download-go-annotations/
```

## Ontology Files

### OBO Format (Recommended)

```bash
# Standard ontology (most common)
wget http://purl.obolibrary.org/obo/go/go-basic.obo

# Full ontology with all relationships
wget http://purl.obolibrary.org/obo/go.obo

# GO Plus (with axioms and cross-ontology links)
wget http://purl.obolibrary.org/obo/go/extensions/go-plus.owl
```

### OWL Format

```bash
# OWL format
wget http://purl.obolibrary.org/obo/go.owl

# GO-CAM compatible
wget http://purl.obolibrary.org/obo/go/extensions/go-lego.owl
```

### JSON Format

```bash
wget http://purl.obolibrary.org/obo/go.json
```

## Gene Annotations (GAF Files)

### Human Annotations

```bash
# Human (UniProt-GOA)
wget http://current.geneontology.org/annotations/goa_human.gaf.gz

# Human with IEA (electronic annotations)
wget http://current.geneontology.org/annotations/goa_human_isoform.gaf.gz
```

### Model Organisms

```bash
# Mouse (MGI)
wget http://current.geneontology.org/annotations/mgi.gaf.gz

# Rat (RGD)
wget http://current.geneontology.org/annotations/rgd.gaf.gz

# Zebrafish (ZFIN)
wget http://current.geneontology.org/annotations/zfin.gaf.gz

# Drosophila (FlyBase)
wget http://current.geneontology.org/annotations/fb.gaf.gz

# C. elegans (WormBase)
wget http://current.geneontology.org/annotations/wb.gaf.gz

# Yeast (SGD)
wget http://current.geneontology.org/annotations/sgd.gaf.gz

# Arabidopsis (TAIR)
wget http://current.geneontology.org/annotations/tair.gaf.gz

# E. coli (EcoCyc)
wget http://current.geneontology.org/annotations/ecocyc.gaf.gz
```

### All Species

```bash
# UniProt-GOA (all species)
wget http://current.geneontology.org/annotations/goa_uniprot_all.gaf.gz
```

## REST API

### Base URL

```
https://api.geneontology.org/api
```

### Get Term Details

```bash
# Get GO term
curl "https://api.geneontology.org/api/ontology/term/GO:0008150"

# Get term with relations
curl "https://api.geneontology.org/api/ontology/term/GO:0008150/graph"
```

### Search Terms

```bash
# Search by name
curl "https://api.geneontology.org/api/search/entity?q=apoptosis&category=ontology_class"

# Search with filters
curl "https://api.geneontology.org/api/search/entity?q=apoptosis&category=ontology_class&rows=100"
```

### Get Gene Annotations

```bash
# Get annotations for gene (by HGNC)
curl "https://api.geneontology.org/api/bioentity/gene/HGNC:11998/function"

# Get annotations for UniProt protein
curl "https://api.geneontology.org/api/bioentity/gene/UniProtKB:P04637/function"
```

### Get Genes for Term

```bash
# Get genes annotated to term
curl "https://api.geneontology.org/api/bioentity/function/GO:0006915/genes?rows=1000"
```

## AmiGO 2 Browser

### Web Interface

```
https://amigo.geneontology.org/amigo
```

### Direct Links

```bash
# Term page
https://amigo.geneontology.org/amigo/term/GO:0006915

# Gene page
https://amigo.geneontology.org/amigo/gene_product/UniProtKB:P04637

# Search
https://amigo.geneontology.org/amigo/search/bioentity?q=TP53
```

## QuickGO (EBI)

### Web Interface

```
https://www.ebi.ac.uk/QuickGO/
```

### REST API

```bash
# Get term info
curl "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0006915"

# Get annotations
curl "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=P04637"

# Download annotations as TSV
curl "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId=P04637" -H "Accept: text/tsv"
```

## SPARQL Endpoint

### URL

```
http://rdf.geneontology.org/sparql
```

### Example Queries

```sparql
# Get all children of a term
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX obo: <http://purl.obolibrary.org/obo/>

SELECT ?child ?label WHERE {
  ?child rdfs:subClassOf obo:GO_0006915 ;
         rdfs:label ?label .
}
```

```sparql
# Get annotations for a gene
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX obo: <http://purl.obolibrary.org/obo/>

SELECT ?term ?label WHERE {
  ?annotation rdf:type obo:GO_0003674 ;
              obo:RO_0002333 <http://identifiers.org/uniprot/P04637> ;
              obo:RO_0002331 ?term .
  ?term rdfs:label ?label .
}
```

## Python Examples

### Parse OBO File

```python
from goatools import obo_parser

# Load ontology
go = obo_parser.GODag('go-basic.obo')

# Get term
term = go['GO:0006915']
print(f"{term.id}: {term.name}")
print(f"Namespace: {term.namespace}")
print(f"Definition: {term.defn}")

# Get ancestors
ancestors = term.get_all_parents()
print(f"Ancestors: {len(ancestors)}")

# Get children
children = term.children
print(f"Children: {len(children)}")
```

### Parse GAF File

```python
from goatools.associations import read_gaf

# Load annotations
gene2go = read_gaf('goa_human.gaf.gz')

# Get GO terms for gene
if 'P04637' in gene2go:
    terms = gene2go['P04637']
    print(f"P04637 annotated to {len(terms)} terms")
```

### Enrichment Analysis

```python
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag

# Load ontology and annotations
go = GODag('go-basic.obo')
gene2go = read_gaf('goa_human.gaf.gz')

# Background gene set
background = set(gene2go.keys())

# Study gene set
study_genes = {'P04637', 'P38398', 'Q13315', ...}

# Run enrichment
goeaobj = GOEnrichmentStudy(
    background,
    gene2go,
    go,
    propagate_counts=True,
    alpha=0.05,
    methods=['fdr_bh']
)

results = goeaobj.run_study(study_genes)

# Filter significant
significant = [r for r in results if r.p_fdr_bh < 0.05]
```

### API Access

```python
import requests

class GOClient:
    def __init__(self):
        self.base_url = "https://api.geneontology.org/api"

    def get_term(self, go_id):
        """Get GO term details."""
        response = requests.get(f"{self.base_url}/ontology/term/{go_id}")
        return response.json()

    def search_terms(self, query, limit=100):
        """Search for GO terms."""
        response = requests.get(
            f"{self.base_url}/search/entity",
            params={
                "q": query,
                "category": "ontology_class",
                "rows": limit
            }
        )
        return response.json()

    def get_annotations(self, gene_id):
        """Get GO annotations for gene."""
        response = requests.get(
            f"{self.base_url}/bioentity/gene/{gene_id}/function"
        )
        return response.json()

# Example
client = GOClient()
term = client.get_term("GO:0006915")
print(f"Term: {term['label']}")
```

## GO Slims

Pre-defined subsets for high-level analysis:

```bash
# Generic slim
wget http://current.geneontology.org/ontology/subsets/goslim_generic.obo

# ChEMBL drug discovery slim
wget http://current.geneontology.org/ontology/subsets/goslim_chembl.obo

# Plant slim
wget http://current.geneontology.org/ontology/subsets/goslim_plant.obo
```

---

## Dataset Versions

### Current Release: 2024-01-17

| Property | Value |
|----------|-------|
| Version | 2024-01-17 |
| Release Date | 2024-01-17 |
| GO Terms | ~45,000 |
| Annotations (Human) | ~700,000 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| go-basic.obo | ~35 MB | ~45K terms | Standard ontology |
| go.owl | ~150 MB | ~45K terms | OWL format |
| goa_human.gaf.gz | ~20 MB | ~700K | Human annotations |
| goa_uniprot_all.gaf.gz | ~5 GB | ~400M | All species |

### Previous Versions

| Version | Release | GO Terms | Status |
|---------|---------|----------|--------|
| 2023-12-15 | 2023-12-15 | ~44,800 | Archived |
| 2023-11-15 | 2023-11-15 | ~44,700 | Archived |
| 2023-10-15 | 2023-10-15 | ~44,600 | Archived |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| GO API | `https://api.geneontology.org/api` |
| QuickGO | `https://www.ebi.ac.uk/QuickGO/services` |
| SPARQL | `http://rdf.geneontology.org/sparql` |
| Authentication | None required |
| Rate Limit | No strict limit |

### GO API Endpoints

| Operation | Endpoint | Example |
|-----------|----------|---------|
| Get Term | `/ontology/term/{id}` | `/ontology/term/GO:0008150` |
| Search | `/search/entity` | `?q=apoptosis&category=ontology_class` |
| Gene Annotations | `/bioentity/gene/{id}/function` | `/bioentity/gene/UniProtKB:P04637/function` |
| Genes for Term | `/bioentity/function/{id}/genes` | `/bioentity/function/GO:0006915/genes` |

### QuickGO Endpoints

| Operation | Endpoint | Example |
|-----------|----------|---------|
| Term Info | `/ontology/go/terms/{id}` | `/ontology/go/terms/GO:0006915` |
| Annotations | `/annotation/search` | `?geneProductId=P04637` |
| Download TSV | `/annotation/downloadSearch` | `?geneProductId=P04637` |

---

## Update Frequency

| Resource | Frequency |
|----------|-----------|
| Ontology | Monthly |
| Annotations | Monthly |
| API | Real-time |

## See Also

- [Schema Documentation](./schema.md)
- [GO Documentation](http://geneontology.org/docs/)
