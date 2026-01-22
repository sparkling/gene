# Research Papers: Coverage Analysis & Multi-Source Strategy

**Date:** January 2026
**Purpose:** Document actual coverage gaps and multi-source bulk download strategy

---

## Executive Summary

**PubMed alone is NOT sufficient.** For comprehensive genetics coverage, use multiple sources with official bulk downloads (no scraping required).

| Source | Papers | Genetics Coverage | Bulk Download |
|--------|--------|-------------------|---------------|
| PubMed | 39M | Excellent (95%) | FTP - 28GB |
| OpenAlex | 250M | Good (80%) | S3 - 300GB |
| PMC OA | 3.4M | Full-text subset | FTP/S3 - 50GB |

---

## 1. Coverage Reality

### What PubMed Covers Well

| Topic | Coverage | Notes |
|-------|----------|-------|
| SNP/GWAS studies | ~95% | Core strength |
| Pharmacogenomics | ~90% | Strong NIH funding |
| Clinical genetics | ~85% | Most major journals indexed |
| Human disease genetics | ~90% | Comprehensive |

### What PubMed Misses

| Topic | Coverage | Gap |
|-------|----------|-----|
| Supplements/nutrition | ~60% | Many non-indexed journals |
| Traditional medicine (TCM) | ~30% | Regional databases not indexed |
| Ayurveda research | ~25% | Indian journals often missing |
| Preprints | ~5% | Not indexed until published |
| International journals | ~70% | Many non-English journals missing |
| Conference papers | ~20% | Rarely indexed |

### Coverage by Journal Type

| Journal Type | In PubMed | Notes |
|--------------|-----------|-------|
| Major biomedical (Nature, Science, Cell) | 100% | Complete |
| Clinical journals | ~95% | Very good |
| Regional biomedical | ~60% | Variable |
| Nutrition/supplement journals | ~50% | Many gaps |
| Integrative medicine | ~40% | Poor coverage |
| Non-English journals | ~30% | Significant gap |

---

## 2. Multi-Source Strategy

### Recommended Sources (All with Bulk Download)

```
┌─────────────────────────────────────────────────────────────────┐
│                 MULTI-SOURCE COVERAGE STRATEGY                   │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  LAYER 1: PubMed (Primary - Genetics Focus)                     │
│  ├── 39M biomedical papers                                      │
│  ├── Best MeSH term filtering for genetics                      │
│  ├── Structured abstracts                                       │
│  └── Download: FTP baseline (28GB compressed)                   │
│                                                                 │
│  LAYER 2: OpenAlex (Breadth - Fill Coverage Gaps)               │
│  ├── 250M works (6x PubMed)                                     │
│  ├── Includes preprints, conference papers, books               │
│  ├── Better international coverage                              │
│  ├── Citation graph included                                    │
│  └── Download: S3 snapshot (300GB)                              │
│                                                                 │
│  LAYER 3: PMC Open Access (Depth - Full Text)                   │
│  ├── 3.4M full-text articles                                    │
│  ├── Complete methods, results, figures                         │
│  ├── For high-impact papers needing detail                      │
│  └── Download: FTP/S3 (50GB)                                    │
│                                                                 │
│  LAYER 4: Specialized (Domain-Specific)                         │
│  ├── Europe PMC: 43M with text-mined annotations                │
│  ├── Semantic Scholar: AI-extracted entities                    │
│  └── Preprint servers: bioRxiv, medRxiv                         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 3. Source Comparison

### Full Comparison Table

| Feature | PubMed | OpenAlex | PMC OA | Europe PMC | Semantic Scholar |
|---------|--------|----------|--------|------------|------------------|
| **Total Records** | 39M | 250M | 3.4M | 43M | 200M |
| **Abstracts** | Yes | Yes* | Yes | Yes | Yes |
| **Full Text** | No | No | Yes | Partial | No |
| **MeSH Terms** | Yes | No | Yes | Yes | No |
| **Concepts/Topics** | MeSH | OpenAlex concepts | MeSH | MeSH | S2 fields |
| **Citations** | Limited | Full graph | Yes | Yes | Full graph |
| **Bulk Download** | FTP | S3 | FTP/S3 | FTP | Request |
| **API Rate Limit** | 10/sec | 100K/day | 10/sec | 10/sec | 100/sec |
| **License** | Public domain | CC0 | Per-article | Mixed | Academic |
| **Update Frequency** | Daily | Monthly | Weekly | Weekly | Weekly |

*OpenAlex has abstracts where available from source

### Genetics-Specific Coverage

| Genetics Topic | PubMed | OpenAlex | Best Source |
|----------------|--------|----------|-------------|
| SNP associations | 95% | 80% | PubMed |
| GWAS studies | 98% | 85% | PubMed |
| Pharmacogenomics | 90% | 75% | PubMed |
| Gene expression | 92% | 85% | PubMed |
| Nutrigenomics | 60% | 70% | OpenAlex |
| Epigenetics | 88% | 82% | PubMed |
| TCM + genetics | 30% | 50% | OpenAlex |
| Preprints | 5% | 40% | OpenAlex |

---

## 4. Bulk Download Methods (No Scraping)

### 4.1 PubMed Baseline (Official FTP)

**Location:** `ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/`

```bash
#!/bin/bash
# Download PubMed baseline - Official NCBI method
# No API key required, no rate limits for FTP

PUBMED_DIR="/data/pubmed"
mkdir -p $PUBMED_DIR/baseline

# Download all baseline files (~28GB, 1,219 files)
cd $PUBMED_DIR/baseline
wget --recursive --no-parent --no-directories \
     --accept "*.xml.gz" \
     ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/

# Download checksums and verify
wget ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/*.md5
md5sum -c *.md5

# Each file contains ~30,000 citations
# Total: ~39 million citations
```

**Daily Updates:**
```bash
# Download daily update files
wget --recursive --no-parent --no-directories \
     --accept "*.xml.gz" \
     ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/
```

**Restrictions:** None. NCBI explicitly supports bulk FTP download.

### 4.2 OpenAlex Snapshot (AWS S3)

**Location:** `s3://openalex`

```bash
#!/bin/bash
# Download OpenAlex snapshot - Official method
# No authentication required, CC0 license

OPENALEX_DIR="/data/openalex"
mkdir -p $OPENALEX_DIR

# Option 1: Full sync (~300GB)
aws s3 sync s3://openalex --no-sign-request $OPENALEX_DIR/

# Option 2: Works only (papers)
aws s3 sync s3://openalex/data/works --no-sign-request $OPENALEX_DIR/works/

# Option 3: Use manifest for selective download
curl -o manifest.json https://openalex.s3.amazonaws.com/data/works/manifest

# Data is partitioned by updated_date
# Format: JSON lines, gzipped
```

**Structure:**
```
s3://openalex/
├── data/
│   ├── works/           # 250M papers (~250GB)
│   ├── authors/         # 90M authors
│   ├── institutions/    # 100K institutions
│   ├── concepts/        # 65K concepts
│   ├── venues/          # 125K journals
│   └── sources/         # Source metadata
└── RELEASE_NOTES.txt
```

**Restrictions:** None. CC0 license, designed for bulk download.

### 4.3 PMC Open Access (FTP + S3)

**FTP Location:** `ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/`

```bash
#!/bin/bash
# Download PMC Open Access subset

PMC_DIR="/data/pmc"
mkdir -p $PMC_DIR

# Option 1: FTP (~50GB)
wget --recursive --no-parent \
     ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/

# Option 2: AWS S3 (faster)
aws s3 sync s3://pmc-oa-opendata --no-sign-request $PMC_DIR/

# Get file list with license info
wget ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv
```

**License Tracking Required:**
```csv
# oa_file_list.csv format:
File,Article Citation,AccessionID,LastUpdated,PMID,License
oa_comm/PMC0001/PMC001.tar.gz,"...",PMC001,2024-01-01,12345,CC BY
oa_noncomm/PMC0002/PMC002.tar.gz,"...",PMC002,2024-01-01,12346,CC BY-NC
```

**Restrictions:** Must track and respect per-article CC licenses.

### 4.4 Europe PMC (FTP)

**Location:** `ftp://ftp.ebi.ac.uk/pub/databases/pmc/`

```bash
#!/bin/bash
# Download Europe PMC data

EUROPMC_DIR="/data/europepmc"
mkdir -p $EUROPMC_DIR

# Abstracts and metadata
wget --recursive --no-parent \
     ftp://ftp.ebi.ac.uk/pub/databases/pmc/TextMinedTerms/

# Full text (OA subset)
wget --recursive --no-parent \
     ftp://ftp.ebi.ac.uk/pub/databases/pmc/oa/
```

**Restrictions:** None for text-mined data. Per-article licenses for full text.

### 4.5 Crossref (Metadata Dumps)

**For DOI resolution and citation metadata:**

```bash
# Public data file (metadata plus)
# Request from: https://www.crossref.org/documentation/retrieve-metadata/rest-api/
# Academic/research use: Free
# Commercial: Requires membership

# API for selective queries (50 requests/sec with polite pool)
curl "https://api.crossref.org/works?query=genetics&rows=1000" \
     -H "User-Agent: YourApp/1.0 (mailto:you@example.com)"
```

---

## 5. Download Size & Time Estimates

### Raw Download Sizes

| Source | Compressed | Uncompressed | Download Time* |
|--------|------------|--------------|----------------|
| PubMed baseline | 28 GB | 100 GB | 2-4 hours |
| OpenAlex works | 250 GB | 800 GB | 8-12 hours |
| PMC OA subset | 50 GB | 200 GB | 3-5 hours |
| Europe PMC annotations | 10 GB | 40 GB | 1 hour |
| **Total** | **338 GB** | **1.1 TB** | **~20 hours** |

*On 100 Mbps connection

### After Filtering (Genetics Only)

| Source | Raw | Filtered | Reduction |
|--------|-----|----------|-----------|
| PubMed | 100 GB | 15 GB | 85% |
| OpenAlex | 800 GB | 100 GB | 87% |
| PMC OA | 200 GB | 20 GB | 90% |
| **Total** | **1.1 TB** | **135 GB** | **88%** |

### Final Storage (RuVector)

| Configuration | Papers | Raw Text | Embeddings | Total |
|---------------|--------|----------|------------|-------|
| PubMed only | 3-5M | 5 GB | 1 GB | ~6 GB |
| + OpenAlex | 10-15M | 15 GB | 4 GB | ~19 GB |
| + PMC full-text | +1M | +10 GB | +2 GB | ~31 GB |

---

## 6. Restrictions Summary

| Source | Bulk Download | Rate Limits | License | Commercial Use |
|--------|---------------|-------------|---------|----------------|
| **PubMed FTP** | Unlimited | None | Public domain | YES |
| **OpenAlex S3** | Unlimited | None | CC0 | YES |
| **PMC OA FTP** | Unlimited | None | Per-article | CHECK EACH |
| **Europe PMC FTP** | Unlimited | None | Mixed | CHECK EACH |
| **Crossref API** | 50/sec | Polite pool | CC0 metadata | YES |
| **Semantic Scholar** | By request | 100/sec | Academic | NO* |

*Semantic Scholar requires approval for commercial use

---

## 7. Recommended Implementation

### Phase 1: MVP (Week 1-2)
```
Download: PubMed baseline (28GB)
Filter: Genetics MeSH terms → 3-5M papers
Store: Abstracts + embeddings in RuVector
Result: ~6GB database, excellent genetics coverage
```

### Phase 2: Expand Coverage (Week 3-4)
```
Download: OpenAlex works snapshot (250GB)
Filter: Genetics concepts → 10-15M additional papers
Dedupe: Against PubMed by DOI/PMID
Store: Fill coverage gaps (preprints, international)
Result: ~20GB database, comprehensive coverage
```

### Phase 3: Add Depth (Week 5-6)
```
Download: PMC OA subset (50GB)
Filter: High-citation genetics papers
Store: Full-text for detailed RAG
Result: ~30GB database, abstracts + key full-text
```

---

## 8. Quick Start Commands

```bash
#!/bin/bash
# Complete download script for genetics knowledge base

# Create directories
mkdir -p /data/{pubmed,openalex,pmc}

# 1. PubMed (primary - do this first)
echo "Downloading PubMed baseline..."
cd /data/pubmed
wget -r -np -nd -A "*.xml.gz" ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/
echo "PubMed complete: $(ls *.gz | wc -l) files"

# 2. OpenAlex (secondary - broader coverage)
echo "Downloading OpenAlex works..."
aws s3 sync s3://openalex/data/works --no-sign-request /data/openalex/works/
echo "OpenAlex complete"

# 3. PMC OA (tertiary - full text)
echo "Downloading PMC Open Access..."
aws s3 sync s3://pmc-oa-opendata --no-sign-request /data/pmc/
echo "PMC complete"

echo "All downloads complete!"
du -sh /data/*
```

---

## 9. API Access (For Real-Time Queries)

When bulk data isn't enough, use APIs for real-time queries:

| API | Rate Limit | Best For |
|-----|------------|----------|
| PubMed E-utilities | 10/sec (with key) | SNP-specific searches |
| OpenAlex | 100K/day | Broad searches, citations |
| Europe PMC | 10/sec | Text-mined entities |
| Semantic Scholar | 100/sec | Citation graphs |

```python
# Example: Real-time search across sources
import requests

def search_papers(query: str, sources: list = ['pubmed', 'openalex']):
    results = []

    if 'pubmed' in sources:
        # PubMed E-utilities
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {'db': 'pubmed', 'term': query, 'retmax': 100, 'retmode': 'json'}
        resp = requests.get(url, params=params)
        results.extend(resp.json().get('esearchresult', {}).get('idlist', []))

    if 'openalex' in sources:
        # OpenAlex API
        url = f"https://api.openalex.org/works"
        params = {'search': query, 'per_page': 100}
        resp = requests.get(url, params=params)
        results.extend([w['id'] for w in resp.json().get('results', [])])

    return results
```

---

## Summary

| Question | Answer |
|----------|--------|
| Is PubMed enough? | No - misses ~40% of relevant papers |
| Best additional source? | OpenAlex (250M papers, CC0, S3 bulk download) |
| Need to scrape? | **NO** - all sources have official bulk downloads |
| Restrictions? | Minimal - PubMed/OpenAlex are public domain/CC0 |
| Total download size? | ~350GB raw, ~135GB after genetics filtering |
| Final database size? | ~6-30GB depending on scope |
