---
id: download-wikipedia
title: "Wikipedia Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# Wikipedia - Download Documentation

## Overview

Wikipedia provides multiple access methods: REST API for content retrieval, MediaWiki API for detailed queries, and database dumps for bulk access.

## REST API (Recommended for Content)

### Base URL

```
https://en.wikipedia.org/api/rest_v1
```

### Get Page Summary

```bash
# Get summary with extract
curl "https://en.wikipedia.org/api/rest_v1/page/summary/TP53"

# Get specific language
curl "https://de.wikipedia.org/api/rest_v1/page/summary/P53"
```

### Get Full HTML

```bash
# Get parsed HTML
curl "https://en.wikipedia.org/api/rest_v1/page/html/TP53"

# Get specific revision
curl "https://en.wikipedia.org/api/rest_v1/page/html/TP53/1234567890"
```

### Get Mobile-Optimized Content

```bash
# Get mobile sections
curl "https://en.wikipedia.org/api/rest_v1/page/mobile-sections/TP53"

# Get lead section only
curl "https://en.wikipedia.org/api/rest_v1/page/mobile-sections-lead/TP53"
```

### Get Related Pages

```bash
# Get related articles
curl "https://en.wikipedia.org/api/rest_v1/page/related/TP53"
```

## MediaWiki API

### Base URL

```
https://en.wikipedia.org/w/api.php
```

### Query Page Content

```bash
# Get page content (wikitext)
curl "https://en.wikipedia.org/w/api.php?action=query&titles=TP53&prop=revisions&rvprop=content&rvslots=main&format=json"

# Get parsed HTML
curl "https://en.wikipedia.org/w/api.php?action=parse&page=TP53&format=json"

# Get plain text extract
curl "https://en.wikipedia.org/w/api.php?action=query&titles=TP53&prop=extracts&exintro&explaintext&format=json"
```

### Get Page Metadata

```bash
# Get categories
curl "https://en.wikipedia.org/w/api.php?action=query&titles=TP53&prop=categories&format=json"

# Get internal links
curl "https://en.wikipedia.org/w/api.php?action=query&titles=TP53&prop=links&pllimit=500&format=json"

# Get external links
curl "https://en.wikipedia.org/w/api.php?action=query&titles=TP53&prop=extlinks&format=json"

# Get Wikidata ID
curl "https://en.wikipedia.org/w/api.php?action=query&titles=TP53&prop=pageprops&format=json"
```

### Search

```bash
# Full-text search
curl "https://en.wikipedia.org/w/api.php?action=query&list=search&srsearch=tumor+suppressor&format=json"

# Prefix search (autocomplete)
curl "https://en.wikipedia.org/w/api.php?action=opensearch&search=TP5&format=json"
```

### Bulk Queries

```bash
# Multiple pages
curl "https://en.wikipedia.org/w/api.php?action=query&titles=TP53|BRCA1|EGFR&prop=extracts&exintro&explaintext&format=json"

# Generator (pages in category)
curl "https://en.wikipedia.org/w/api.php?action=query&generator=categorymembers&gcmtitle=Category:Tumor_suppressor_genes&prop=extracts&exintro&format=json"
```

## Database Dumps

### Dump URLs

```
https://dumps.wikimedia.org/
https://dumps.wikimedia.org/enwiki/
```

### Available Dumps

| File | Content | Size |
|------|---------|------|
| pages-articles.xml.bz2 | Current article text | ~20 GB |
| pages-meta-current.xml.bz2 | Current + metadata | ~25 GB |
| pages-meta-history.xml.bz2 | Full history | ~1 TB |
| categorylinks.sql.gz | Category relationships | ~1 GB |
| pagelinks.sql.gz | Internal links | ~6 GB |
| externallinks.sql.gz | External links | ~3 GB |

### Download Commands

```bash
# Download latest article dump
wget https://dumps.wikimedia.org/enwiki/latest/enwiki-latest-pages-articles.xml.bz2

# Download category links
wget https://dumps.wikimedia.org/enwiki/latest/enwiki-latest-categorylinks.sql.gz

# Download page table
wget https://dumps.wikimedia.org/enwiki/latest/enwiki-latest-page.sql.gz
```

### Process Dump

```python
import bz2
import xml.etree.ElementTree as ET

def process_dump(filepath):
    """Process Wikipedia XML dump."""
    with bz2.open(filepath, 'rt') as f:
        for event, elem in ET.iterparse(f, events=['end']):
            if elem.tag.endswith('page'):
                title = elem.find('.//{*}title')
                text = elem.find('.//{*}revision/{*}text')

                if title is not None and text is not None:
                    yield {
                        'title': title.text,
                        'text': text.text
                    }

                elem.clear()
```

## DBpedia

### SPARQL Endpoint

```
https://dbpedia.org/sparql
```

### Query Examples

```sparql
# Get gene information
SELECT ?gene ?label ?entrez ?uniprot
WHERE {
  ?gene a dbo:Gene ;
        rdfs:label ?label ;
        dbo:entrezgene ?entrez ;
        dbo:uniprot ?uniprot .
  FILTER (lang(?label) = 'en')
}
LIMIT 100
```

### DBpedia Downloads

```bash
# Download DBpedia datasets
wget https://downloads.dbpedia.org/repo/dbpedia/generic/infobox-properties/2023.09.01/infobox-properties_lang=en.ttl.bz2
```

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| REST API | No hard limit (be reasonable) |
| MediaWiki API | maxlag parameter recommended |
| Dumps | Unlimited |

### User-Agent Requirement

All API requests should include a User-Agent header:

```bash
curl -H "User-Agent: MyApp/1.0 (contact@example.com)" \
  "https://en.wikipedia.org/api/rest_v1/page/summary/TP53"
```

## Python Examples

### REST API Access

```python
import requests

def get_wikipedia_extract(title, lang='en'):
    """Get article summary from Wikipedia."""
    url = f"https://{lang}.wikipedia.org/api/rest_v1/page/summary/{title}"
    headers = {"User-Agent": "MyApp/1.0 (contact@example.com)"}

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()
        return {
            'title': data['title'],
            'extract': data.get('extract', ''),
            'wikidata_id': data.get('wikibase_item'),
            'url': data['content_urls']['desktop']['page']
        }
    return None
```

### Batch Extraction

```python
import requests
import time

def get_gene_articles(gene_list):
    """Get Wikipedia extracts for multiple genes."""
    results = []

    for gene in gene_list:
        data = get_wikipedia_extract(gene)
        if data:
            results.append(data)
        time.sleep(0.1)  # Be polite

    return results

# Example
genes = ['TP53', 'BRCA1', 'EGFR', 'KRAS']
articles = get_gene_articles(genes)
```

## Other Language Wikis

| Language | API Base |
|----------|----------|
| English | en.wikipedia.org |
| German | de.wikipedia.org |
| French | fr.wikipedia.org |
| Spanish | es.wikipedia.org |
| Japanese | ja.wikipedia.org |
| Chinese | zh.wikipedia.org |

## See Also

- [Schema Documentation](./schema.md)
- [MediaWiki API Help](https://www.mediawiki.org/wiki/API:Main_page)
- [Wikimedia Downloads](https://dumps.wikimedia.org/)
