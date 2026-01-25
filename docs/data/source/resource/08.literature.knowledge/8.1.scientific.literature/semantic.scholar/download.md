---
id: download-semantic-scholar
title: "Semantic Scholar Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# Semantic Scholar - Download Documentation

## Overview

Semantic Scholar provides a REST API with free tier access and bulk dataset downloads for research purposes.

## REST API

### Base URL

```
https://api.semanticscholar.org
```

### Paper Endpoints

```bash
# Get paper by S2 ID
curl "https://api.semanticscholar.org/graph/v1/paper/649def34f8be52c8b66281af98ae884c09aef38b"

# Get paper by DOI
curl "https://api.semanticscholar.org/graph/v1/paper/DOI:10.1038/s41586-019-1666-5"

# Get paper by ArXiv ID
curl "https://api.semanticscholar.org/graph/v1/paper/ARXIV:1910.10683"

# Get paper by PMID
curl "https://api.semanticscholar.org/graph/v1/paper/PMID:31534227"

# Get paper by CorpusId
curl "https://api.semanticscholar.org/graph/v1/paper/CorpusId:123456789"
```

### Specify Fields

```bash
# Select specific fields
curl "https://api.semanticscholar.org/graph/v1/paper/649def34...?fields=title,authors,year,citationCount,abstract"

# Get citations with context
curl "https://api.semanticscholar.org/graph/v1/paper/649def34.../citations?fields=title,authors,year,contexts,intents"

# Get references
curl "https://api.semanticscholar.org/graph/v1/paper/649def34.../references?fields=title,authors,year,isInfluential"
```

### Paper Search

```bash
# Keyword search
curl "https://api.semanticscholar.org/graph/v1/paper/search?query=CRISPR+gene+editing&limit=100"

# Search with year filter
curl "https://api.semanticscholar.org/graph/v1/paper/search?query=deep+learning&year=2020-2024"

# Search with fields of study
curl "https://api.semanticscholar.org/graph/v1/paper/search?query=transformer&fieldsOfStudy=Computer+Science"

# Search with venue
curl "https://api.semanticscholar.org/graph/v1/paper/search?query=neural+networks&venue=Nature"

# Bulk paper lookup (up to 500)
curl -X POST "https://api.semanticscholar.org/graph/v1/paper/batch" \
  -H "Content-Type: application/json" \
  -d '{"ids": ["649def34...", "abc123..."]}'
```

### Author Endpoints

```bash
# Get author by ID
curl "https://api.semanticscholar.org/graph/v1/author/1741101?fields=name,affiliations,paperCount,citationCount,hIndex"

# Get author's papers
curl "https://api.semanticscholar.org/graph/v1/author/1741101/papers?fields=title,year,citationCount&limit=100"

# Author search
curl "https://api.semanticscholar.org/graph/v1/author/search?query=John+Doe"
```

### Recommendations

```bash
# Get recommended papers based on a paper
curl "https://api.semanticscholar.org/recommendations/v1/papers/forpaper/649def34...?fields=title,authors,year"

# Get recommendations for multiple positive examples
curl -X POST "https://api.semanticscholar.org/recommendations/v1/papers/" \
  -H "Content-Type: application/json" \
  -d '{"positivePaperIds": ["649def34...", "abc123..."]}'
```

## Available Fields

### Paper Fields

| Field | Description |
|-------|-------------|
| paperId | S2 paper ID |
| corpusId | Numeric corpus ID |
| externalIds | DOI, ArXiv, PMID, etc. |
| url | S2 URL |
| title | Paper title |
| abstract | Abstract text |
| venue | Publication venue |
| publicationVenue | Venue details |
| year | Publication year |
| referenceCount | Number of references |
| citationCount | Citation count |
| influentialCitationCount | Influential citations |
| isOpenAccess | OA status |
| openAccessPdf | OA PDF URL |
| fieldsOfStudy | Research fields |
| s2FieldsOfStudy | S2-assigned fields |
| publicationTypes | Article types |
| publicationDate | Full date |
| journal | Journal details |
| authors | Author list |
| citations | Citing papers |
| references | Referenced papers |
| tldr | AI summary |
| embedding | SPECTER vector |

### Author Fields

| Field | Description |
|-------|-------------|
| authorId | S2 author ID |
| externalIds | ORCID, DBLP |
| name | Display name |
| aliases | Alternative names |
| affiliations | Institutions |
| homepage | Personal website |
| paperCount | Total papers |
| citationCount | Total citations |
| hIndex | h-index |

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| Anonymous | 100 requests/5 minutes |
| API Key | 1 request/second |

### Getting an API Key

Request at: https://www.semanticscholar.org/product/api

```bash
# Use API key
curl -H "x-api-key: YOUR_KEY" "https://api.semanticscholar.org/graph/v1/paper/search?query=example"
```

## Bulk Dataset Downloads

### Dataset API

```bash
# Get latest release info
curl "https://api.semanticscholar.org/datasets/v1/release/latest"

# List available datasets
curl "https://api.semanticscholar.org/datasets/v1/release/2024-01-02"

# Get download links for papers dataset
curl -H "x-api-key: YOUR_KEY" "https://api.semanticscholar.org/datasets/v1/release/2024-01-02/dataset/papers"
```

### Available Datasets

| Dataset | Description | Size |
|---------|-------------|------|
| papers | All paper metadata | ~100 GB |
| authors | Author information | ~10 GB |
| abstracts | Paper abstracts | ~30 GB |
| citations | Citation links | ~50 GB |
| tldrs | AI summaries | ~5 GB |
| embeddings | SPECTER vectors | ~200 GB |

### Download Example

```bash
# Get download manifest
manifest=$(curl -H "x-api-key: YOUR_KEY" \
  "https://api.semanticscholar.org/datasets/v1/release/2024-01-02/dataset/papers")

# Download each file
for url in $(echo $manifest | jq -r '.files[].url'); do
  wget "$url"
done
```

## Python Examples

### API Access

```python
import requests
import time

class SemanticScholarClient:
    def __init__(self, api_key=None):
        self.base_url = "https://api.semanticscholar.org/graph/v1"
        self.headers = {"x-api-key": api_key} if api_key else {}
        self.delay = 1.0 if api_key else 3.0

    def get_paper(self, paper_id, fields=None):
        """Get paper by ID."""
        url = f"{self.base_url}/paper/{paper_id}"
        params = {"fields": ",".join(fields)} if fields else {}

        response = requests.get(url, headers=self.headers, params=params)
        time.sleep(self.delay)
        return response.json()

    def search_papers(self, query, limit=100, fields=None):
        """Search for papers."""
        url = f"{self.base_url}/paper/search"
        params = {"query": query, "limit": limit}
        if fields:
            params["fields"] = ",".join(fields)

        response = requests.get(url, headers=self.headers, params=params)
        time.sleep(self.delay)
        return response.json()

    def get_citations(self, paper_id, limit=1000):
        """Get all citations for a paper."""
        url = f"{self.base_url}/paper/{paper_id}/citations"
        citations = []
        offset = 0

        while True:
            params = {
                "fields": "title,authors,year,citationCount,contexts,intents",
                "limit": 100,
                "offset": offset
            }
            response = requests.get(url, headers=self.headers, params=params)
            data = response.json()

            if not data.get("data"):
                break

            citations.extend(data["data"])
            offset += 100

            if len(citations) >= limit or offset >= data.get("total", 0):
                break

            time.sleep(self.delay)

        return citations[:limit]
```

### Process Bulk Data

```python
import gzip
import json

def process_papers_file(filepath):
    """Process a gzipped JSONL papers file."""
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            paper = json.loads(line)
            yield paper

# Example: Find highly cited papers
for paper in process_papers_file('papers-part0.jsonl.gz'):
    if paper.get('citationCount', 0) > 1000:
        print(f"{paper['title']}: {paper['citationCount']} citations")
```

## Update Frequency

| Data Type | Update Frequency |
|-----------|------------------|
| API | Weekly |
| Datasets | Weekly releases |
| Citation counts | Weekly |
| TLDRs | Continuous |

## See Also

- [Schema Documentation](./schema.md)
- [S2 API Documentation](https://api.semanticscholar.org/api-docs/)
