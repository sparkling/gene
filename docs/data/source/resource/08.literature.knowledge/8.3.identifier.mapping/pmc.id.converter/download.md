---
id: download-pmc-id-converter
title: "PMC ID Converter Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# PMC ID Converter - Download Documentation

## Overview

The PMC ID Converter service provides batch conversion between PMID, PMCID, DOI, and Manuscript IDs. It supports up to 200 IDs per request.

## API Endpoint

```
https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/
```

## Basic Conversions

### PMID to PMCID/DOI

```bash
# Single PMID
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678&format=json"

# Multiple PMIDs
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678,23456789,34567890&format=json"
```

### PMCID to PMID/DOI

```bash
# With PMC prefix
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=PMC1234567&format=json"

# Without prefix (specify idtype)
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=1234567&idtype=pmcid&format=json"
```

### DOI to PMID/PMCID

```bash
# Single DOI
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=10.1038/example&idtype=doi&format=json"

# Multiple DOIs
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=10.1038/example1,10.1038/example2&idtype=doi&format=json"
```

### Manuscript ID Conversion

```bash
# NIHMS ID
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=NIHMS123456&format=json"

# With MID prefix
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=MID:NIHMS123456&format=json"
```

## Parameters

| Parameter | Required | Values | Description |
|-----------|----------|--------|-------------|
| ids | Yes | string | Comma-separated identifiers (max 200) |
| idtype | No | pmid, pmcid, doi, mid | Input type (auto-detected if omitted) |
| format | No | json, xml, csv | Response format (default: xml) |
| versions | No | yes, no | Include version history |
| tool | No | string | Application name |
| email | No | string | Contact email |

## Output Formats

### JSON Format

```bash
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678&format=json"
```

### XML Format

```bash
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678&format=xml"
```

### CSV Format

```bash
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678&format=csv"
```

## Include Version History

```bash
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=PMC1234567&versions=yes&format=json"
```

## Batch Processing

### Maximum Batch Size

Up to 200 IDs per request.

```bash
# Large batch
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678,12345679,12345680,...&format=json"
```

### Mixed ID Types

The service auto-detects ID types:

```bash
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678,PMC1234567,10.1038/example&format=json"
```

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| General | No strict limit |
| Recommended | 3 requests/second |
| With tool/email | Preferred for tracking |

### Best Practices

```bash
# Include identification
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678&format=json&tool=MyApp&email=user@example.com"
```

## Python Examples

### Basic Conversion

```python
import requests
import time

def convert_pmids_to_pmcids(pmids):
    """Convert list of PMIDs to PMCIDs."""
    base_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"

    results = {}
    batch_size = 200

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]

        params = {
            "ids": ",".join(map(str, batch)),
            "format": "json",
            "tool": "MyApp",
            "email": "user@example.com"
        }

        response = requests.get(base_url, params=params)
        data = response.json()

        if "records" in data:
            for record in data["records"]:
                pmid = record.get("pmid")
                pmcid = record.get("pmcid")
                if pmid and pmcid:
                    results[pmid] = pmcid

        time.sleep(0.35)  # Rate limiting

    return results
```

### Bidirectional Conversion

```python
def convert_ids(ids, target_type="pmcid"):
    """Convert IDs to target type."""
    base_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"

    params = {
        "ids": ",".join(ids),
        "format": "json",
        "tool": "MyApp",
        "email": "user@example.com"
    }

    response = requests.get(base_url, params=params)
    data = response.json()

    results = {}
    for record in data.get("records", []):
        if "errmsg" not in record:
            requested = record.get("requested-id", record.get("pmid", ""))
            target = record.get(target_type)
            if target:
                results[requested] = target

    return results

# Examples
pmid_to_pmcid = convert_ids(["12345678", "23456789"], "pmcid")
pmcid_to_doi = convert_ids(["PMC1234567"], "doi")
doi_to_pmid = convert_ids(["10.1038/example"], "pmid")
```

### Handle Errors

```python
def safe_convert(ids):
    """Convert IDs with error handling."""
    base_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"

    params = {
        "ids": ",".join(ids),
        "format": "json"
    }

    response = requests.get(base_url, params=params)
    data = response.json()

    results = {"success": {}, "errors": {}}

    for record in data.get("records", []):
        requested_id = record.get("requested-id", "")

        if "errmsg" in record:
            results["errors"][requested_id] = record["errmsg"]
        else:
            results["success"][requested_id] = {
                "pmid": record.get("pmid"),
                "pmcid": record.get("pmcid"),
                "doi": record.get("doi")
            }

    return results
```

## Command Line Usage

### Using curl and jq

```bash
# Convert and extract PMCIDs
curl -s "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678,23456789&format=json" | jq '.records[] | {pmid, pmcid}'

# Get DOIs only
curl -s "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=PMC1234567&format=json" | jq -r '.records[].doi'
```

### Batch from File

```bash
# IDs in file (one per line)
ids=$(cat pmids.txt | tr '\n' ',' | sed 's/,$//')
curl "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=${ids}&format=csv" > results.csv
```

## Integration with PMC Downloads

After converting to PMCIDs, use the IDs to download full text:

```bash
# 1. Convert PMID to PMCID
pmcid=$(curl -s "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=12345678&format=json" | jq -r '.records[0].pmcid')

# 2. Download full text XML
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=${pmcid}&rettype=xml" > article.xml
```

## See Also

- [Schema Documentation](./schema.md)
- [PubMed Central](../../8.1.scientific.literature/pubmed.central/README.md)
