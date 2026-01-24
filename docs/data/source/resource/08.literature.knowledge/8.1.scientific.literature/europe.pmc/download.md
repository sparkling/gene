---
id: download-europe-pmc
title: "Europe PMC Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# Europe PMC - Download Documentation

## Overview

Europe PMC provides multiple access methods: REST API for programmatic queries, OAI-PMH for metadata harvesting, and FTP for bulk downloads.

## REST API

### Base URL

```
https://www.ebi.ac.uk/europepmc/webservices/rest
```

### Search Endpoint

```bash
# Basic search
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TP53&format=json"

# With pagination
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TP53&pageSize=100&cursorMark=*&format=json"

# Filtered by source
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TP53+SRC:MED&format=json"

# Open access only
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TP53+OPEN_ACCESS:Y&format=json"

# By date range
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TP53+FIRST_PDATE:[2020-01-01+TO+2024-12-31]&format=json"
```

### Article Retrieval

```bash
# Get article by PMID
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:12345678+SRC:MED&format=json"

# Get full-text XML
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1234567/fullTextXML"

# Get references
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/MED/12345678/references?format=json"

# Get citations
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/MED/12345678/citations?format=json"
```

## Annotations API

```bash
# Get text-mined annotations
curl "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=PMC:PMC1234567&format=JSON"

# Filter by annotation type
curl "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=PMC:PMC1234567&type=Gene_Proteins&format=JSON"
```

## OAI-PMH Harvesting

### Base URL

```
https://www.ebi.ac.uk/europepmc/oai.cgi
```

### Commands

```bash
# Identify repository
curl "https://www.ebi.ac.uk/europepmc/oai.cgi?verb=Identify"

# List metadata formats
curl "https://www.ebi.ac.uk/europepmc/oai.cgi?verb=ListMetadataFormats"

# List records (incremental harvesting)
curl "https://www.ebi.ac.uk/europepmc/oai.cgi?verb=ListRecords&metadataPrefix=pmc&from=2024-01-01"

# Get specific record
curl "https://www.ebi.ac.uk/europepmc/oai.cgi?verb=GetRecord&identifier=PMC1234567&metadataPrefix=pmc"
```

## FTP Bulk Downloads

### Base URL

```
https://europepmc.org/ftp/
```

### Available Files

| Directory | Content | Update Frequency |
|-----------|---------|------------------|
| /oa | Open Access full text | Weekly |
| /manuscripts | Author manuscripts | Weekly |
| /suppl_data | Supplementary files | As available |

### Download Commands

```bash
# Download OA package list
curl -O "https://europepmc.org/ftp/oa/oa_comm_xml.PMC000xxxxxx.baseline.YYYY-MM-DD.filelist.txt"

# Download specific package
curl -O "https://europepmc.org/ftp/oa/oa_comm_xml.PMC000xxxxxx.baseline.YYYY-MM-DD.tar.gz"
```

## Query Syntax

### Field Prefixes

| Prefix | Field | Example |
|--------|-------|---------|
| TITLE: | Article title | `TITLE:cancer` |
| ABSTRACT: | Abstract text | `ABSTRACT:mutation` |
| AUTH: | Author name | `AUTH:Smith` |
| JOURNAL: | Journal name | `JOURNAL:"Nature"` |
| PMID: | PubMed ID | `PMID:12345678` |
| DOI: | DOI | `DOI:10.1000/example` |
| PUB_YEAR: | Publication year | `PUB_YEAR:2024` |
| SRC: | Source type | `SRC:MED` |

### Boolean Operators

```
TP53 AND cancer
TP53 OR p53
TP53 NOT "cell line"
(TP53 OR BRCA1) AND mutation
```

## Rate Limits

| Access Type | Rate Limit |
|-------------|------------|
| Anonymous | 10 requests/second |
| With email parameter | 25 requests/second |
| Registered API key | 50 requests/second |

## Python Example

```python
import requests
import time

def search_europepmc(query, max_results=1000):
    """Search Europe PMC with pagination."""
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    results = []
    cursor = "*"

    while len(results) < max_results:
        params = {
            "query": query,
            "format": "json",
            "pageSize": 100,
            "cursorMark": cursor,
            "resultType": "core"
        }

        response = requests.get(base_url, params=params)
        data = response.json()

        if "resultList" not in data or not data["resultList"]["result"]:
            break

        results.extend(data["resultList"]["result"])
        cursor = data.get("nextCursorMark", cursor)

        if cursor == data.get("cursorMark"):
            break

        time.sleep(0.1)  # Rate limiting

    return results[:max_results]
```

## Output Formats

| Format | Parameter | Use Case |
|--------|-----------|----------|
| JSON | `format=json` | Programmatic access |
| XML | `format=xml` | Standard exchange |
| Dublin Core | `format=dc` | Metadata harvesting |

## See Also

- [Schema Documentation](./schema.md)
- [Europe PMC API Documentation](https://europepmc.org/RestfulWebService)
