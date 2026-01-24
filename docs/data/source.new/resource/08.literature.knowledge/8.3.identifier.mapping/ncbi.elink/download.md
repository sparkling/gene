---
id: download-ncbi-elink
title: "NCBI E-Link Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# NCBI E-Link - Download Documentation

## Overview

NCBI E-Link provides programmatic access to discover links between records across 40+ NCBI databases. It is part of the Entrez E-utilities suite.

## Base URL

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi
```

## Basic Link Discovery

### Find Links Between Databases

```bash
# PubMed to Gene links
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gene&id=12345678&retmode=json"

# Gene to PubMed links
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=pubmed&id=7157&retmode=json"

# Protein to Structure links
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=structure&id=P04637&retmode=json"
```

### Multiple Input IDs

```bash
# Batch query (comma-separated)
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gene&id=12345678,23456789,34567890&retmode=json"

# One LinkSet per ID
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gene&id=12345678&id=23456789&retmode=json"
```

## Command Modes

### cmd=neighbor (Related Records)

```bash
# Find related PubMed articles
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pubmed&id=12345678&cmd=neighbor&retmode=json"

# With relevance scores
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pubmed&id=12345678&cmd=neighbor_score&retmode=json"
```

### cmd=acheck (Available Link Types)

```bash
# Check what links are available
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12345678&cmd=acheck&retmode=json"
```

### cmd=llinks (External Links)

```bash
# Get LinkOut URLs
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12345678&cmd=llinks&retmode=json"

# Get external links for genes
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&id=7157&cmd=llinks&retmode=json"
```

### cmd=llinkslib (Library Links)

```bash
# Get library-specific links
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12345678&cmd=llinkslib&retmode=json"
```

### cmd=prlinks (Primary Links)

```bash
# Get primary provider links
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=12345678&cmd=prlinks&retmode=json"
```

## Link Names

### Specific Link Type

```bash
# Only pubmed_gene_rif links (Gene Reference Into Function)
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gene&linkname=pubmed_gene_rif&id=12345678&retmode=json"

# Only pubmed_pubmed_citedin (citing articles)
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pubmed&linkname=pubmed_pubmed_citedin&id=12345678&retmode=json"
```

### List Available Link Names

Use the EInfo endpoint:

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=pubmed&retmode=json"
```

## Using History Server

### Store Results

```bash
# Store link results on history server
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gene&id=12345678&cmd=neighbor_history&retmode=json"
```

### Use with ESearch

```bash
# Search and then link
# Step 1: Search
search_result=$(curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=TP53+cancer&usehistory=y&retmode=json")

# Extract WebEnv and QueryKey
webenv=$(echo $search_result | jq -r '.esearchresult.webenv')
query_key=$(echo $search_result | jq -r '.esearchresult.querykey')

# Step 2: Link using history
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gene&query_key=${query_key}&WebEnv=${webenv}&retmode=json"
```

## Filter Parameters

### Date Range

```bash
# Links with date filter
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pubmed&id=12345678&cmd=neighbor&datetype=pdat&mindate=2020&maxdate=2024&retmode=json"
```

### Taxonomy Filter

```bash
# Human-specific links
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id=7157&linkname=gene_protein&term=txid9606[ORGN]&retmode=json"
```

## Rate Limits

| Access Type | Rate Limit |
|-------------|------------|
| Without API key | 3 requests/second |
| With API key | 10 requests/second |

### Using API Key

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gene&id=12345678&api_key=YOUR_API_KEY&retmode=json"
```

Register for API key: https://www.ncbi.nlm.nih.gov/account/

## Python Examples

### Basic Link Discovery

```python
from Bio import Entrez
import time

Entrez.email = "your@email.com"
Entrez.api_key = "your_api_key"

def get_gene_links(pmid):
    """Get genes linked to a PubMed article."""
    handle = Entrez.elink(
        dbfrom="pubmed",
        db="gene",
        id=pmid
    )
    result = Entrez.read(handle)
    handle.close()

    genes = []
    for linkset in result:
        if "LinkSetDb" in linkset:
            for db in linkset["LinkSetDb"]:
                if db["LinkName"] == "pubmed_gene":
                    genes.extend([link["Id"] for link in db["Link"]])

    return genes
```

### Batch Processing

```python
def batch_link_discovery(pmids, target_db="gene", batch_size=100):
    """Get links for multiple PMIDs."""
    all_links = {}

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]

        handle = Entrez.elink(
            dbfrom="pubmed",
            db=target_db,
            id=batch
        )
        result = Entrez.read(handle)
        handle.close()

        for linkset in result:
            from_id = linkset["IdList"][0]
            links = []
            if "LinkSetDb" in linkset:
                for db in linkset["LinkSetDb"]:
                    links.extend([link["Id"] for link in db["Link"]])
            all_links[from_id] = links

        time.sleep(0.1)

    return all_links
```

### Find Related Articles

```python
def find_related_articles(pmid, max_results=20):
    """Find related PubMed articles with scores."""
    handle = Entrez.elink(
        dbfrom="pubmed",
        db="pubmed",
        id=pmid,
        cmd="neighbor_score"
    )
    result = Entrez.read(handle)
    handle.close()

    related = []
    for linkset in result:
        if "LinkSetDb" in linkset:
            for db in linkset["LinkSetDb"]:
                if db["LinkName"] == "pubmed_pubmed":
                    for link in db["Link"][:max_results]:
                        related.append({
                            "pmid": link["Id"],
                            "score": int(link.get("Score", 0))
                        })

    return sorted(related, key=lambda x: x["score"], reverse=True)
```

## EDirect Command Line

```bash
# Install EDirect
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Use elink command
esearch -db pubmed -query "TP53 cancer" | elink -target gene | efetch -format docsum

# Chain multiple links
esearch -db gene -query "7157" | elink -target protein | elink -target structure | efetch -format docsum
```

## See Also

- [Schema Documentation](./schema.md)
- [E-utilities Help](https://www.ncbi.nlm.nih.gov/books/NBK25500/)
