---
id: download-pubmed-central
title: "PubMed Central (PMC) Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# PubMed Central (PMC) - Download Documentation

## Overview

PubMed Central provides access via E-utilities API, OAI-PMH harvesting, and FTP bulk downloads. The Open Access Subset allows commercial text mining.

## E-Utilities API

### Base URL

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
```

### Search PMC

```bash
# Search by keyword
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=TP53+cancer&retmode=json"

# Search by date
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=cancer&mindate=2024/01/01&maxdate=2024/12/31&datetype=pdat&retmode=json"

# Search Open Access subset
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=open+access[filter]+AND+cancer&retmode=json"
```

### Fetch Article XML

```bash
# Fetch by PMCID
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMC1234567&rettype=xml"

# Fetch multiple
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMC1234567,PMC1234568&rettype=xml"
```

### ID Conversion

```bash
# Convert PMID to PMCID
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678&format=json"

# Convert DOI to PMCID
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=10.1038/ncomms12345&format=json"

# Batch conversion (up to 200)
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678,12345679,12345680&format=json"
```

## FTP Bulk Downloads

### FTP Server

```
ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/
```

### Available Directories

| Directory | Content | Access |
|-----------|---------|--------|
| /oa_bulk/ | Open Access articles | Commercial OK |
| /oa_pdf/ | OA PDF files | Commercial OK |
| /oa_comm/ | Commercial use subset | Commercial OK |
| /oa_noncomm/ | Non-commercial OA | NC only |
| /manuscript/ | Author manuscripts | Public |

### Download OA Subset

```bash
# Get file list
wget ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.txt

# Download specific package
wget ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/xml/PMC000xxxxxx.baseline.2024-01-01.tar.gz

# Mirror entire OA subset (large!)
wget -r -np ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/
```

### File List Format

```
oa_file_list.txt columns:
1. File path
2. Journal
3. PMCID
4. License
5. Last updated

Example:
oa_package/08/11/PMC81234.tar.gz<TAB>J Biol Chem<TAB>PMC81234<TAB>CC BY<TAB>2024-01-15
```

## OAI-PMH Harvesting

### Base URL

```
https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi
```

### Commands

```bash
# Identify
curl "https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=Identify"

# List sets
curl "https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=ListSets"

# Harvest records
curl "https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=ListRecords&metadataPrefix=pmc&from=2024-01-01"

# Get specific record
curl "https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:1234567&metadataPrefix=pmc"
```

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| Without API key | 3 requests/second |
| With API key | 10 requests/second |
| Bulk FTP | Unlimited |

### Getting an API Key

Register at: https://www.ncbi.nlm.nih.gov/account/

```bash
# Use API key
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=cancer&api_key=YOUR_KEY"
```

## Python Examples

### Search and Download

```python
from Bio import Entrez
import time

Entrez.email = "your@email.com"
Entrez.api_key = "your_api_key"

def search_pmc(query, max_results=1000):
    """Search PMC and retrieve results."""
    # Search
    handle = Entrez.esearch(db="pmc", term=query, retmax=max_results)
    results = Entrez.read(handle)
    pmcids = results["IdList"]

    # Fetch in batches
    articles = []
    for i in range(0, len(pmcids), 100):
        batch = pmcids[i:i+100]
        handle = Entrez.efetch(db="pmc", id=batch, rettype="xml")
        articles.append(handle.read())
        time.sleep(0.1)

    return articles
```

### Process Downloaded XML

```python
from lxml import etree

def extract_abstract(xml_file):
    """Extract abstract from PMC XML."""
    tree = etree.parse(xml_file)

    # Extract abstract paragraphs
    abstracts = tree.xpath('//abstract//p/text()')
    return ' '.join(abstracts)

def extract_body_text(xml_file):
    """Extract full body text."""
    tree = etree.parse(xml_file)

    # Get all paragraph text
    paragraphs = tree.xpath('//body//p')
    text = []
    for p in paragraphs:
        text.append(etree.tostring(p, method='text', encoding='unicode'))

    return '\n\n'.join(text)
```

## Data Size Estimates

| Content | Approximate Size |
|---------|------------------|
| OA Subset (all) | ~500 GB compressed |
| Commercial use subset | ~300 GB compressed |
| Author manuscripts | ~100 GB compressed |
| PDF subset | ~1 TB |

## Update Frequency

| Method | Frequency |
|--------|-----------|
| FTP baseline | Monthly |
| FTP updates | Daily |
| API | Real-time |

## License Information

| Subset | License | Commercial Use |
|--------|---------|----------------|
| oa_comm | CC BY, CC0 | Yes |
| oa_noncomm | CC BY-NC | No |
| manuscript | Public access | Research only |

## See Also

- [Schema Documentation](./schema.md)
- [PMC Open Access Subset](https://www.ncbi.nlm.nih.gov/pmc/tools/openftlist/)
