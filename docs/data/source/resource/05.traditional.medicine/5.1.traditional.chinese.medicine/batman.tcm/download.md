---
id: download-batman-tcm
title: "BATMAN-TCM 2.0 Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# BATMAN-TCM 2.0 Download Instructions

## Quick Start

```bash
# REST API for programmatic access
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=herb&query=Ginseng&output=json"

# Web interface for bulk downloads
# URL: http://bionet.ncpsb.org.cn/batman-tcm/
```

## Prerequisites

- Python 3.8+ with requests library
- cURL or similar HTTP client
- Sufficient storage for bulk data (~500MB)

## Download Methods

### Method 1: REST API (Recommended)

BATMAN-TCM 2.0 provides full REST API access:

```bash
# Query herb targets
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=herb&query=Ginseng&output=json"

# Query compound targets
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=compound&query=ginsenoside_Rg1&output=json"

# Reverse query: find herbs for gene set
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=genes&genes=TP53,EGFR,MTOR&output=json"
```

### Method 2: Bulk Download

Tab-delimited files available via web interface:
- **URL**: http://bionet.ncpsb.org.cn/batman-tcm/
- **Download Section**: Bulk TSV files

## File Inventory

| File | Records | Description |
|------|---------|-------------|
| formulas.tsv | 54,832 | TCM formula data |
| herbs.tsv | 8,404 | Herb metadata |
| ingredients.tsv | 39,171 | Compound information |
| known_tti.tsv | 17,068 | Known target interactions |
| predicted_tti.tsv | 2,319,272 | Predicted interactions |
| targets.tsv | ~15,000 | Protein annotations |

## API Parameters

| Parameter | Values | Description |
|-----------|--------|-------------|
| query_type | herb, compound, formula, genes | Query mode |
| query | String | Search term |
| output | json, html | Response format |

## Python Download Script

```python
import requests
import time
import json

BASE_URL = "http://bionet.ncpsb.org.cn/batman-tcm/api"

def query_with_retry(params, retries=3, delay=5):
    """Query API with retry logic for server issues."""
    for attempt in range(retries):
        try:
            response = requests.get(BASE_URL, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except (requests.Timeout, requests.ConnectionError) as e:
            if attempt < retries - 1:
                print(f"Retry {attempt + 1}/{retries} after error: {e}")
                time.sleep(delay * (attempt + 1))
            else:
                raise

# Query herb targets
params = {"query_type": "herb", "query": "Ginseng", "output": "json"}
result = query_with_retry(params)
print(f"Found {len(result.get('results', []))} targets")

# Add delay between requests (server may rate limit)
time.sleep(3)
```

## Verification

After download, verify against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Formulas | 54,832 |
| Herbs | 8,404 |
| Ingredients | 39,171 |
| Known TTIs | 17,068 |
| Predicted TTIs | 2,319,272 |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Current Version | BATMAN-TCM 2.0 (2023) |
| Update Frequency | Major version releases |
| Notification | Publications and website |

---

## Dataset Versions

### Current Release: BATMAN-TCM 2.0

| Property | Value |
|----------|-------|
| Version | 2.0 |
| Release Date | 2023-01-01 |
| Total Size | ~500 MB |
| Prediction AUC | 0.9663 |

### Version Contents

| Component | Records | Description |
|-----------|---------|-------------|
| formulas.tsv | 54,832 | TCM formula data |
| herbs.tsv | 8,404 | Herb metadata |
| ingredients.tsv | 39,171 | Compound information |
| known_tti.tsv | 17,068 | Known target interactions |
| predicted_tti.tsv | 2,319,272 | Predicted interactions |

### Previous Versions

| Version | Release | Status |
|---------|---------|--------|
| BATMAN-TCM 1.0 | 2016-01-01 | Archived |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| Base URL | `http://bionet.ncpsb.org.cn/batman-tcm/api` |
| Authentication | None required |
| Rate Limit | 3+ seconds between requests |
| Response Format | JSON, HTML |

### API Parameters

| Parameter | Values | Description |
|-----------|--------|-------------|
| query_type | herb, compound, formula, genes | Query mode |
| query | String | Search term |
| output | json, html | Response format |

### Example Queries

| Query Type | URL |
|------------|-----|
| Herb targets | `?query_type=herb&query=Ginseng&output=json` |
| Compound targets | `?query_type=compound&query=ginsenoside_Rg1&output=json` |
| Reverse query | `?query_type=genes&genes=TP53,EGFR&output=json` |

---

## Notes

- Server may have connectivity issues; implement retry logic
- Use reasonable delays between API calls (3+ seconds)
- Licensed CC BY-NC 4.0 (non-commercial)
- Prediction ROC AUC: 0.9663
- Most comprehensive TCM database with API access
