# Data Processing Pipeline for Bulk Downloads

**Last Updated:** January 2026
**Purpose:** Comprehensive guide for downloading, processing, and loading data from all public sources into the database.

---

## Table of Contents

1. [Storage Requirements](#1-storage-requirements)
2. [Processing Architecture](#2-processing-architecture)
3. [Streaming Parsers](#3-streaming-parsers)
4. [Memory-Efficient Processing Patterns](#4-memory-efficient-processing-patterns)
5. [Parallel Processing](#5-parallel-processing)
6. [Database Loading](#6-database-loading)
7. [ID Mapping Pipeline](#7-id-mapping-pipeline)
8. [Incremental Updates](#8-incremental-updates)
9. [Recommended Tools](#9-recommended-tools)
10. [Example Pipeline Scripts](#10-example-pipeline-scripts)

---

## 1. Storage Requirements

### 1.1 Total Download Sizes by Source

| Category | Source | Compressed | Uncompressed | Format |
|----------|--------|------------|--------------|--------|
| **SNP/Variants** | dbSNP (full) | 15 GB | 150 GB | VCF |
| | dbSNP (practical) | 500 MB | 5 GB | JSON |
| | ClinVar | 500 MB | 5 GB | XML/VCF |
| | SNPedia | 100 MB | 1 GB | Wiki export |
| | GWAS Catalog | 100 MB | 1 GB | TSV |
| | gnomAD (summary) | 50 GB | 500 GB | VCF |
| **Genes/Proteins** | UniProt (human) | 500 MB | 5 GB | XML |
| | NCBI Gene | 500 MB | 5 GB | XML |
| | GeneCards | 500 MB | 5 GB | JSON |
| | Ensembl | 5 GB | 50 GB | GFF/VCF |
| **Pathways** | Reactome | 200 MB | 2 GB | BioPAX |
| | WikiPathways | 100 MB | 1 GB | GPML |
| **Pharmacogenomics** | PharmGKB | 50 MB | 500 MB | TSV |
| | ChEMBL | 500 MB | 3 GB | SQL/SDF |
| | Open Targets | 1 GB | 10 GB | Parquet |
| **Supplements** | NIH ODS | 50 MB | 500 MB | JSON |
| | DSLD | 500 MB | 5 GB | JSON |
| **Traditional Medicine** | BATMAN-TCM 2.0 | 300 MB | 3 GB | REST API |
| | IMPPAT | 100 MB | 1 GB | SDF |
| | TCMBank | 300 MB | 3 GB | Download |
| **Research** | PubMed baseline | 100 GB | 500 GB | XML |
| | PubMed (abstracts) | 20 GB | 100 GB | Text |

### 1.2 Total Storage by Scope

| Scope | Download Size | Working Space | Processed Data | Total Required |
|-------|---------------|---------------|----------------|----------------|
| **MVP** | 20 GB | 40 GB | 3 GB | **~65 GB** |
| **Standard** | 50 GB | 100 GB | 12 GB | **~165 GB** |
| **Comprehensive** | 200 GB | 400 GB | 30 GB | **~650 GB** |
| **Full (with gnomAD)** | 700 GB | 1.4 TB | 120 GB | **~2.2 TB** |

### 1.3 Storage Recommendations

| Storage Type | Use Case | Recommended |
|--------------|----------|-------------|
| **SSD (NVMe)** | Active processing, database | 500 GB - 1 TB |
| **HDD** | Raw downloads archive, backups | 2-4 TB |
| **Cloud Object Storage** | Long-term archive, large datasets | S3/GCS |

**Directory Layout:**

```
/data/
├── downloads/          # Raw compressed downloads (HDD ok)
│   ├── ncbi/
│   ├── ebi/
│   └── ...
├── raw/                # Decompressed raw files (SSD preferred)
│   └── {source}/{date}/
├── processed/          # Cleaned, normalized data (SSD)
│   ├── snps/
│   ├── genes/
│   └── ...
├── staging/            # Temporary processing space (SSD required)
└── database/           # Database files (SSD required)
```

### 1.4 Database Storage Requirements

| Configuration | RuVector | PostgreSQL | Total |
|---------------|----------|------------|-------|
| MVP (1M SNPs) | 800 MB | 500 MB | ~1.5 GB |
| Standard (5M SNPs) | 3.5 GB | 2 GB | ~6 GB |
| Comprehensive (10M SNPs) | 7 GB | 5 GB | ~12 GB |

---

## 2. Processing Architecture

### 2.1 Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              DATA SOURCES                                    │
│   NCBI FTP | EBI FTP | AWS S3 | Hugging Face | Web APIs | Web Scraping     │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           DOWNLOAD LAYER                                     │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐     │
│  │  wget/   │  │  AWS CLI │  │ Hugging  │  │  API     │  │  Web     │     │
│  │  curl    │  │  (S3)    │  │  Face    │  │ Clients  │  │ Scraper  │     │
│  └──────────┘  └──────────┘  └──────────┘  └──────────┘  └──────────┘     │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         RAW DATA STORAGE                                     │
│                    /data/downloads/{source}/{date}/                          │
│                         (Compressed archives)                                │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      DECOMPRESSION LAYER                                     │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐                    │
│  │   pigz   │  │   lbzip2 │  │   xz     │  │  unzip   │                    │
│  │ (gzip)   │  │  (bzip2) │  │  (lzma)  │  │  (zip)   │                    │
│  └──────────┘  └──────────┘  └──────────┘  └──────────┘                    │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      STREAMING PARSER LAYER                                  │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐                    │
│  │  ijson   │  │  lxml    │  │  pandas  │  │  VCF     │                    │
│  │  (JSON)  │  │  (XML)   │  │ (CSV/TSV)│  │  Parser  │                    │
│  └──────────┘  └──────────┘  └──────────┘  └──────────┘                    │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        FILTER LAYER                                          │
│  - Remove irrelevant records (e.g., non-human, deprecated)                  │
│  - Apply quality thresholds                                                  │
│  - Deduplicate entries                                                       │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       TRANSFORM LAYER                                        │
│  - Normalize identifiers (rs numbers, gene symbols)                         │
│  - Standardize formats (dates, coordinates)                                  │
│  - Merge multi-source records                                                │
│  - Generate text for embeddings                                              │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         LOAD LAYER                                           │
│  ┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐          │
│  │    PostgreSQL    │  │     RuVector     │  │   Relationship   │          │
│  │   (COPY bulk)    │  │  (batch insert)  │  │     Builder      │          │
│  └──────────────────┘  └──────────────────┘  └──────────────────┘          │
└─────────────────────────────────────┬───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         INDEX LAYER                                          │
│  - Build HNSW vector indexes                                                 │
│  - Create B-tree indexes on properties                                       │
│  - Build graph relationship indexes                                          │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 2.2 Processing Phases

| Phase | Description | Parallelizable | I/O Bound |
|-------|-------------|----------------|-----------|
| Download | Fetch from sources | Yes (per source) | Network |
| Decompress | Uncompress archives | Yes (pigz) | Disk |
| Parse | Stream through files | Partially | Disk |
| Filter | Apply selection criteria | Yes | CPU |
| Transform | Normalize and enrich | Yes | CPU |
| Load | Insert into database | Partially | Disk/Network |
| Index | Build search indexes | No | CPU/Memory |

---

## 3. Streaming Parsers

### 3.1 JSON Streaming (ijson)

For large JSON files (UniProt, API responses):

```python
import ijson
import gzip

def stream_json_array(filepath: str, item_prefix: str = "item"):
    """Stream JSON array items without loading entire file into memory."""

    # Handle gzipped files
    if filepath.endswith('.gz'):
        f = gzip.open(filepath, 'rb')
    else:
        f = open(filepath, 'rb')

    try:
        # Stream items from JSON array
        for item in ijson.items(f, item_prefix):
            yield item
    finally:
        f.close()

# Usage: Process 50GB UniProt JSON with <100MB RAM
for protein in stream_json_array('uniprot_human.json.gz', 'results.item'):
    process_protein(protein)
```

### 3.2 XML Streaming (lxml.iterparse)

For large XML files (PubMed, ClinVar, dbSNP):

```python
from lxml import etree
import gzip

def stream_xml_elements(filepath: str, tag: str):
    """Stream XML elements without loading entire file."""

    if filepath.endswith('.gz'):
        f = gzip.open(filepath, 'rb')
    else:
        f = open(filepath, 'rb')

    context = etree.iterparse(f, events=('end',), tag=tag)

    for event, elem in context:
        yield elem

        # Critical: Clear element to free memory
        elem.clear()

        # Also clear preceding siblings
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    f.close()

# Usage: Process 100GB PubMed XML with <500MB RAM
for article in stream_xml_elements('pubmed24n0001.xml.gz', 'PubmedArticle'):
    pmid = article.find('.//PMID').text
    title = article.find('.//ArticleTitle').text
    process_article(pmid, title)
```

### 3.3 TSV/CSV Streaming (pandas chunked)

For tabular data (GWAS Catalog, PharmGKB):

```python
import pandas as pd

def stream_csv_chunks(filepath: str, chunksize: int = 10000, **kwargs):
    """Stream CSV/TSV in chunks."""

    # Auto-detect compression
    compression = 'infer'

    for chunk in pd.read_csv(
        filepath,
        chunksize=chunksize,
        compression=compression,
        **kwargs
    ):
        yield chunk

# Usage: Process 10GB GWAS catalog with <1GB RAM
for chunk in stream_csv_chunks('gwas_catalog.tsv.gz', sep='\t'):
    # Process 10,000 rows at a time
    filtered = chunk[chunk['P-VALUE'] < 5e-8]
    process_associations(filtered)
```

### 3.4 VCF Streaming (cyvcf2)

For variant call format files (dbSNP, gnomAD, ClinVar):

```python
from cyvcf2 import VCF

def stream_vcf_variants(filepath: str, region: str = None):
    """Stream VCF variants efficiently."""

    vcf = VCF(filepath)

    if region:
        variants = vcf(region)  # e.g., "chr1:1000000-2000000"
    else:
        variants = vcf

    for variant in variants:
        yield {
            'chrom': variant.CHROM,
            'pos': variant.POS,
            'id': variant.ID,
            'ref': variant.REF,
            'alt': variant.ALT,
            'qual': variant.QUAL,
            'filter': variant.FILTER,
            'info': dict(variant.INFO),
        }

    vcf.close()

# Usage: Process gnomAD chromosome by chromosome
for variant in stream_vcf_variants('gnomad.vcf.bgz', region='chr1'):
    if variant['info'].get('AF', 0) > 0.01:  # MAF > 1%
        process_variant(variant)
```

### 3.5 Compressed File Streaming

```python
import gzip
import bz2
import lzma

def open_compressed(filepath: str, mode: str = 'rt'):
    """Open compressed file with appropriate decompressor."""

    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode)
    elif filepath.endswith('.bz2'):
        return bz2.open(filepath, mode)
    elif filepath.endswith('.xz') or filepath.endswith('.lzma'):
        return lzma.open(filepath, mode)
    else:
        return open(filepath, mode.replace('t', ''))

# Usage
with open_compressed('data.json.gz', 'rt') as f:
    for line in f:
        process_line(line)
```

---

## 4. Memory-Efficient Processing Patterns

### 4.1 Generator-Based Pipeline

```python
from typing import Iterator, Generator
from dataclasses import dataclass

@dataclass
class ProcessedSNP:
    rs_number: str
    chromosome: str
    position: int
    gene: str
    significance: str

def parse_snps(filepath: str) -> Generator[dict, None, None]:
    """Parse raw SNP data."""
    for item in stream_json_array(filepath):
        yield item

def filter_snps(snps: Iterator[dict]) -> Generator[dict, None, None]:
    """Filter to clinically relevant SNPs."""
    for snp in snps:
        if snp.get('clinical_significance') not in (None, 'benign'):
            yield snp

def transform_snps(snps: Iterator[dict]) -> Generator[ProcessedSNP, None, None]:
    """Transform to normalized format."""
    for snp in snps:
        yield ProcessedSNP(
            rs_number=f"rs{snp['id']}",
            chromosome=normalize_chrom(snp['seq_id']),
            position=snp['position'],
            gene=snp.get('gene_symbol', ''),
            significance=snp['clinical_significance']
        )

def load_snps(snps: Iterator[ProcessedSNP], batch_size: int = 1000):
    """Load to database in batches."""
    batch = []
    for snp in snps:
        batch.append(snp)
        if len(batch) >= batch_size:
            bulk_insert(batch)
            batch = []
    if batch:
        bulk_insert(batch)

# Complete pipeline - processes 100GB file with <1GB RAM
def process_snp_file(filepath: str):
    raw = parse_snps(filepath)
    filtered = filter_snps(raw)
    transformed = transform_snps(filtered)
    load_snps(transformed)
```

### 4.2 Chunked DataFrame Processing

```python
import pandas as pd
from typing import Callable

def process_large_file(
    filepath: str,
    processor: Callable[[pd.DataFrame], pd.DataFrame],
    output_path: str,
    chunksize: int = 50000,
    **read_kwargs
):
    """Process large file in chunks, writing results incrementally."""

    first_chunk = True
    total_rows = 0

    for chunk in pd.read_csv(filepath, chunksize=chunksize, **read_kwargs):
        # Process chunk
        processed = processor(chunk)

        # Write to output
        processed.to_csv(
            output_path,
            mode='w' if first_chunk else 'a',
            header=first_chunk,
            index=False
        )

        first_chunk = False
        total_rows += len(processed)
        print(f"Processed {total_rows:,} rows")

    return total_rows

# Usage
def filter_gwas(df: pd.DataFrame) -> pd.DataFrame:
    """Filter GWAS associations."""
    return df[
        (df['P-VALUE'] < 5e-8) &
        (df['OR or BETA'].notna())
    ]

process_large_file(
    'gwas_catalog_full.tsv.gz',
    filter_gwas,
    'gwas_filtered.tsv',
    sep='\t'
)
```

### 4.3 Memory-Mapped Processing

```python
import numpy as np
import mmap

def process_with_mmap(filepath: str, dtype: np.dtype, shape: tuple):
    """Process large binary arrays using memory mapping."""

    # Memory-map the file (doesn't load into RAM)
    arr = np.memmap(filepath, dtype=dtype, mode='r', shape=shape)

    # Process in chunks
    chunk_size = 100000
    results = []

    for i in range(0, len(arr), chunk_size):
        chunk = arr[i:i + chunk_size]
        result = process_chunk(chunk)
        results.append(result)

    return np.concatenate(results)

# For embeddings processing
def batch_similarity_search(
    query_embeddings: np.ndarray,
    index_file: str,
    index_shape: tuple,
    batch_size: int = 1000
):
    """Search large embedding index without loading fully into memory."""

    index = np.memmap(index_file, dtype='float32', mode='r', shape=index_shape)

    results = []
    for i in range(0, len(query_embeddings), batch_size):
        batch = query_embeddings[i:i + batch_size]
        # Compute similarities for batch
        sims = np.dot(batch, index.T)
        top_k = np.argsort(sims, axis=1)[:, -10:]
        results.append(top_k)

    return np.vstack(results)
```

### 4.4 Memory Monitoring

```python
import psutil
import os
from functools import wraps

def memory_monitor(func):
    """Decorator to monitor memory usage of a function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())

        mem_before = process.memory_info().rss / 1024 / 1024
        result = func(*args, **kwargs)
        mem_after = process.memory_info().rss / 1024 / 1024

        print(f"{func.__name__}: {mem_before:.1f}MB -> {mem_after:.1f}MB "
              f"(delta: {mem_after - mem_before:+.1f}MB)")

        return result
    return wrapper

def check_memory_threshold(threshold_mb: int = 8000):
    """Raise error if memory usage exceeds threshold."""
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024 / 1024

    if mem_mb > threshold_mb:
        raise MemoryError(
            f"Memory usage ({mem_mb:.0f}MB) exceeds threshold ({threshold_mb}MB)"
        )

# Usage
@memory_monitor
def process_file(filepath):
    for chunk in stream_file(filepath):
        check_memory_threshold(8000)  # 8GB limit
        process_chunk(chunk)
```

---

## 5. Parallel Processing

### 5.1 Multiprocessing for CPU-Bound Tasks

```python
from multiprocessing import Pool, cpu_count
from functools import partial
import os

def process_file_parallel(
    filepaths: list[str],
    processor: callable,
    num_workers: int = None
):
    """Process multiple files in parallel."""

    if num_workers is None:
        num_workers = min(cpu_count(), len(filepaths))

    with Pool(num_workers) as pool:
        results = pool.map(processor, filepaths)

    return results

# Process PubMed XML files in parallel
def process_pubmed_file(filepath: str) -> list[dict]:
    """Process single PubMed XML file."""
    articles = []
    for article in stream_xml_elements(filepath, 'PubmedArticle'):
        articles.append(extract_article(article))
    return articles

# Process all 1,274 baseline files
pubmed_files = glob.glob('/data/pubmed/baseline/*.xml.gz')
all_articles = process_file_parallel(pubmed_files, process_pubmed_file)
```

### 5.2 Asyncio for I/O-Bound Tasks

```python
import asyncio
import aiohttp
import aiofiles
from typing import List

async def download_file(session: aiohttp.ClientSession, url: str, dest: str):
    """Download single file asynchronously."""
    async with session.get(url) as response:
        async with aiofiles.open(dest, 'wb') as f:
            async for chunk in response.content.iter_chunked(8192):
                await f.write(chunk)

async def download_many(urls: List[tuple[str, str]], max_concurrent: int = 10):
    """Download multiple files with concurrency limit."""

    semaphore = asyncio.Semaphore(max_concurrent)

    async def bounded_download(session, url, dest):
        async with semaphore:
            await download_file(session, url, dest)
            print(f"Downloaded: {dest}")

    async with aiohttp.ClientSession() as session:
        tasks = [
            bounded_download(session, url, dest)
            for url, dest in urls
        ]
        await asyncio.gather(*tasks)

# Usage
urls = [
    ('ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz',
     '/data/downloads/dbsnp/00-All.vcf.gz'),
    ('ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz',
     '/data/downloads/clinvar/clinvar.vcf.gz'),
    # ... more files
]

asyncio.run(download_many(urls))
```

### 5.3 Dask for Out-of-Core DataFrames

```python
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster

def setup_dask_cluster(n_workers: int = 4, memory_limit: str = '4GB'):
    """Set up local Dask cluster."""
    cluster = LocalCluster(
        n_workers=n_workers,
        threads_per_worker=2,
        memory_limit=memory_limit
    )
    client = Client(cluster)
    print(f"Dashboard: {client.dashboard_link}")
    return client

def process_with_dask(input_pattern: str, output_path: str):
    """Process large dataset with Dask."""

    client = setup_dask_cluster()

    # Read all matching files
    ddf = dd.read_parquet(input_pattern)

    # Process (lazy evaluation)
    filtered = ddf[ddf['p_value'] < 5e-8]
    aggregated = filtered.groupby('gene_symbol').agg({
        'effect_size': 'mean',
        'study_id': 'count'
    }).reset_index()

    # Execute and save
    aggregated.to_parquet(output_path)

    client.close()

# Usage for gnomAD or Open Targets
process_with_dask(
    '/data/open_targets/*.parquet',
    '/data/processed/gene_associations.parquet'
)
```

### 5.4 GNU Parallel for Shell Tasks

```bash
#!/bin/bash
# Process multiple files with GNU parallel

# Decompress multiple files in parallel
find /data/downloads -name "*.gz" | \
    parallel -j 8 'pigz -d -k {}'

# Process XML files in parallel with Python
find /data/pubmed/baseline -name "*.xml" | \
    parallel -j 4 --progress \
    'python process_pubmed.py {} > /data/processed/pubmed/{/.}.json'

# Download multiple files in parallel
cat download_urls.txt | \
    parallel -j 10 --colsep '\t' \
    'wget -q -O {2} {1}'
```

---

## 6. Database Loading

### 6.1 PostgreSQL COPY for Bulk Inserts

```python
import psycopg2
from io import StringIO
import csv

def bulk_load_postgres(
    connection_string: str,
    table: str,
    data: list[dict],
    columns: list[str]
):
    """Bulk load data using PostgreSQL COPY command."""

    conn = psycopg2.connect(connection_string)
    cursor = conn.cursor()

    # Create CSV buffer
    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction='ignore')

    for row in data:
        writer.writerow(row)

    buffer.seek(0)

    # Use COPY for fast bulk insert
    cursor.copy_expert(
        f"COPY {table} ({','.join(columns)}) FROM STDIN WITH CSV",
        buffer
    )

    conn.commit()
    cursor.close()
    conn.close()

# Usage - 10x faster than INSERT
snps = [...]  # List of SNP dictionaries
bulk_load_postgres(
    'postgresql://user:pass@localhost/gene',
    'snps',
    snps,
    ['rs_number', 'chromosome', 'position', 'ref_allele', 'alt_alleles']
)
```

### 6.2 Streaming COPY from File

```python
def stream_copy_from_file(
    connection_string: str,
    table: str,
    filepath: str,
    columns: list[str],
    delimiter: str = '\t'
):
    """Stream large file directly to PostgreSQL."""

    conn = psycopg2.connect(connection_string)
    cursor = conn.cursor()

    with open(filepath, 'r') as f:
        cursor.copy_expert(
            f"""COPY {table} ({','.join(columns)})
                FROM STDIN
                WITH (FORMAT csv, DELIMITER '{delimiter}', HEADER true)""",
            f
        )

    conn.commit()
    print(f"Loaded {cursor.rowcount:,} rows")
    cursor.close()
    conn.close()

# Load pre-processed TSV directly
stream_copy_from_file(
    'postgresql://localhost/gene',
    'gwas_associations',
    '/data/processed/gwas_filtered.tsv',
    ['snp_id', 'gene', 'trait', 'p_value', 'effect_size']
)
```

### 6.3 SQLite for Development

```python
import sqlite3
from contextlib import contextmanager

@contextmanager
def sqlite_connection(db_path: str):
    """Context manager for SQLite with optimizations."""
    conn = sqlite3.connect(db_path)

    # Performance optimizations
    conn.execute('PRAGMA journal_mode=WAL')
    conn.execute('PRAGMA synchronous=NORMAL')
    conn.execute('PRAGMA cache_size=-64000')  # 64MB cache
    conn.execute('PRAGMA temp_store=MEMORY')

    try:
        yield conn
    finally:
        conn.close()

def bulk_load_sqlite(db_path: str, table: str, data: list[dict]):
    """Bulk load to SQLite with transaction batching."""

    with sqlite_connection(db_path) as conn:
        cursor = conn.cursor()

        # Disable auto-commit for bulk loading
        cursor.execute('BEGIN TRANSACTION')

        columns = list(data[0].keys())
        placeholders = ','.join(['?' for _ in columns])
        sql = f"INSERT OR REPLACE INTO {table} ({','.join(columns)}) VALUES ({placeholders})"

        # Insert in batches
        batch_size = 10000
        for i in range(0, len(data), batch_size):
            batch = data[i:i + batch_size]
            cursor.executemany(sql, [tuple(row.values()) for row in batch])

            if i % 100000 == 0:
                conn.commit()
                cursor.execute('BEGIN TRANSACTION')

        conn.commit()
```

### 6.4 RuVector Batch Insert

```typescript
import { RuVector } from 'ruvector';

async function batchInsertToRuVector(
  db: RuVector,
  collection: string,
  data: any[],
  batchSize: number = 1000
) {
  const total = data.length;
  let inserted = 0;

  for (let i = 0; i < total; i += batchSize) {
    const batch = data.slice(i, i + batchSize);

    await db.batchInsert(collection, batch, {
      generateEmbeddings: true,
      textField: 'text_for_embedding'
    });

    inserted += batch.length;
    console.log(`Inserted ${inserted.toLocaleString()}/${total.toLocaleString()}`);
  }
}

// Usage
const snps = await loadProcessedSNPs();
await batchInsertToRuVector(db, 'snps', snps, 1000);
```

### 6.5 MongoDB for Document Storage

```python
from pymongo import MongoClient, InsertOne
from pymongo.errors import BulkWriteError

def bulk_load_mongodb(
    connection_string: str,
    database: str,
    collection: str,
    documents: list[dict],
    batch_size: int = 1000
):
    """Bulk load documents to MongoDB."""

    client = MongoClient(connection_string)
    db = client[database]
    coll = db[collection]

    operations = [InsertOne(doc) for doc in documents]

    for i in range(0, len(operations), batch_size):
        batch = operations[i:i + batch_size]
        try:
            result = coll.bulk_write(batch, ordered=False)
            print(f"Inserted {result.inserted_count} documents")
        except BulkWriteError as e:
            print(f"Bulk write error: {e.details['nInserted']} inserted")

    client.close()

# For semi-structured traditional medicine data
tcm_compounds = load_tcmsp_data()
bulk_load_mongodb(
    'mongodb://localhost:27017',
    'gene',
    'tcm_compounds',
    tcm_compounds
)
```

---

## 7. ID Mapping Pipeline

### 7.1 Hub-and-Spoke ID Strategy

```
                    ┌─────────────────┐
                    │    UniProt      │
                    │  (Protein Hub)  │
                    └────────┬────────┘
                             │
        ┌────────────────────┼────────────────────┐
        │                    │                    │
        ▼                    ▼                    ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────┐
│  NCBI Gene ID │   │  Ensembl ID   │   │   RefSeq ID   │
└───────────────┘   └───────────────┘   └───────────────┘


                    ┌─────────────────┐
                    │   PubChem CID   │
                    │ (Compound Hub)  │
                    └────────┬────────┘
                             │
        ┌────────────────────┼────────────────────┐
        │                    │                    │
        ▼                    ▼                    ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────┐
│   ChEMBL ID   │   │   DrugBank ID │   │    CAS RN     │
└───────────────┘   └───────────────┘   └───────────────┘
```

### 7.2 Building Mapping Tables

```python
from collections import defaultdict
import pandas as pd

class IDMapper:
    """Centralized ID mapping using hub identifiers."""

    def __init__(self):
        self.protein_map = defaultdict(dict)  # UniProt as hub
        self.compound_map = defaultdict(dict)  # PubChem CID as hub
        self.gene_map = defaultdict(dict)      # NCBI Gene ID as hub

    def load_uniprot_mappings(self, filepath: str):
        """Load UniProt ID mapping file."""

        # UniProt provides mapping files
        # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

        for chunk in pd.read_csv(filepath, sep='\t', chunksize=100000,
                                  names=['uniprot', 'db', 'external_id']):
            for _, row in chunk.iterrows():
                self.protein_map[row['uniprot']][row['db']] = row['external_id']

    def load_pubchem_mappings(self, filepath: str):
        """Load PubChem compound mappings."""

        for chunk in pd.read_csv(filepath, sep='\t', chunksize=100000):
            for _, row in chunk.iterrows():
                cid = row['PubChem_CID']
                if pd.notna(row.get('ChEMBL_ID')):
                    self.compound_map[cid]['chembl'] = row['ChEMBL_ID']
                if pd.notna(row.get('CAS')):
                    self.compound_map[cid]['cas'] = row['CAS']

    def get_gene_id(self, identifier: str, source: str) -> str:
        """Map any gene identifier to NCBI Gene ID."""

        for gene_id, mappings in self.gene_map.items():
            if mappings.get(source) == identifier:
                return gene_id
        return None

    def get_compound_by_any_id(self, identifier: str) -> dict:
        """Find compound by any known identifier."""

        # Check if it's already a PubChem CID
        if identifier.isdigit():
            if identifier in self.compound_map:
                return {'pubchem_cid': identifier, **self.compound_map[identifier]}

        # Search through mappings
        for cid, mappings in self.compound_map.items():
            if identifier in mappings.values():
                return {'pubchem_cid': cid, **mappings}

        return None

# Build mapping tables at pipeline start
mapper = IDMapper()
mapper.load_uniprot_mappings('/data/raw/uniprot/idmapping.dat.gz')
mapper.load_pubchem_mappings('/data/raw/pubchem/compound_mappings.tsv')
```

### 7.3 Cross-Reference Resolution

```python
def resolve_gene_references(snp_data: dict, mapper: IDMapper) -> dict:
    """Resolve all gene identifiers to canonical form."""

    gene_symbol = snp_data.get('gene_symbol')
    gene_id = snp_data.get('gene_id')
    ensembl_id = snp_data.get('ensembl_gene_id')

    # Try to find canonical NCBI Gene ID
    ncbi_gene_id = None

    if gene_id and gene_id.isdigit():
        ncbi_gene_id = gene_id
    elif ensembl_id:
        ncbi_gene_id = mapper.get_gene_id(ensembl_id, 'ensembl')
    elif gene_symbol:
        ncbi_gene_id = mapper.get_gene_id(gene_symbol, 'symbol')

    snp_data['ncbi_gene_id'] = ncbi_gene_id
    return snp_data
```

---

## 8. Incremental Updates

### 8.1 Download State Tracking

```python
import json
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass, asdict

@dataclass
class DownloadState:
    source: str
    last_download: str
    last_version: str
    file_count: int
    total_size: int
    checksum: str

class StateTracker:
    """Track download and processing state."""

    def __init__(self, state_file: str = '/data/state/download_state.json'):
        self.state_file = Path(state_file)
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
        self.state = self._load_state()

    def _load_state(self) -> dict:
        if self.state_file.exists():
            return json.loads(self.state_file.read_text())
        return {}

    def _save_state(self):
        self.state_file.write_text(json.dumps(self.state, indent=2))

    def get_last_download(self, source: str) -> DownloadState:
        if source in self.state:
            return DownloadState(**self.state[source])
        return None

    def update_state(self, state: DownloadState):
        self.state[state.source] = asdict(state)
        self._save_state()

    def needs_update(self, source: str, remote_version: str) -> bool:
        current = self.get_last_download(source)
        if current is None:
            return True
        return current.last_version != remote_version

# Usage
tracker = StateTracker()

if tracker.needs_update('clinvar', get_clinvar_version()):
    download_clinvar()
    tracker.update_state(DownloadState(
        source='clinvar',
        last_download=datetime.now().isoformat(),
        last_version=get_clinvar_version(),
        file_count=1,
        total_size=get_file_size(),
        checksum=compute_checksum()
    ))
```

### 8.2 Diff Processing for Updates

```python
import hashlib
from typing import Set, Tuple

def compute_record_hash(record: dict, key_fields: list[str]) -> str:
    """Compute hash of record for change detection."""
    key_data = '|'.join(str(record.get(f, '')) for f in key_fields)
    return hashlib.md5(key_data.encode()).hexdigest()

def diff_datasets(
    old_data: list[dict],
    new_data: list[dict],
    id_field: str,
    hash_fields: list[str]
) -> Tuple[list[dict], list[dict], list[str]]:
    """Find added, modified, and deleted records."""

    # Build hash maps
    old_hashes = {
        r[id_field]: compute_record_hash(r, hash_fields)
        for r in old_data
    }
    new_hashes = {
        r[id_field]: compute_record_hash(r, hash_fields)
        for r in new_data
    }

    old_ids = set(old_hashes.keys())
    new_ids = set(new_hashes.keys())

    # Find changes
    added_ids = new_ids - old_ids
    deleted_ids = old_ids - new_ids

    modified_ids = {
        id for id in (old_ids & new_ids)
        if old_hashes[id] != new_hashes[id]
    }

    # Get records
    new_by_id = {r[id_field]: r for r in new_data}

    added = [new_by_id[id] for id in added_ids]
    modified = [new_by_id[id] for id in modified_ids]
    deleted = list(deleted_ids)

    return added, modified, deleted

# Usage for ClinVar weekly updates
def process_clinvar_update():
    old_data = load_current_clinvar()
    new_data = parse_new_clinvar_release()

    added, modified, deleted = diff_datasets(
        old_data, new_data,
        id_field='variation_id',
        hash_fields=['clinical_significance', 'review_status', 'last_evaluated']
    )

    print(f"Added: {len(added)}, Modified: {len(modified)}, Deleted: {len(deleted)}")

    # Apply changes
    insert_records(added)
    update_records(modified)
    mark_deleted(deleted)
```

### 8.3 Version Control for Data

```python
from datetime import datetime
import shutil

class DataVersionManager:
    """Manage versioned snapshots of processed data."""

    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.versions_dir = self.base_dir / 'versions'
        self.current_link = self.base_dir / 'current'
        self.versions_dir.mkdir(parents=True, exist_ok=True)

    def create_version(self, source_dir: str, version_tag: str = None) -> str:
        """Create new version snapshot."""

        if version_tag is None:
            version_tag = datetime.now().strftime('%Y%m%d_%H%M%S')

        version_dir = self.versions_dir / version_tag

        # Copy data to versioned directory
        shutil.copytree(source_dir, version_dir)

        # Update current symlink
        if self.current_link.exists():
            self.current_link.unlink()
        self.current_link.symlink_to(version_dir)

        return version_tag

    def rollback(self, version_tag: str):
        """Rollback to previous version."""

        version_dir = self.versions_dir / version_tag
        if not version_dir.exists():
            raise ValueError(f"Version {version_tag} not found")

        if self.current_link.exists():
            self.current_link.unlink()
        self.current_link.symlink_to(version_dir)

    def list_versions(self) -> list[dict]:
        """List all available versions."""

        versions = []
        for d in sorted(self.versions_dir.iterdir()):
            if d.is_dir():
                versions.append({
                    'tag': d.name,
                    'path': str(d),
                    'size': sum(f.stat().st_size for f in d.rglob('*') if f.is_file()),
                    'created': datetime.fromtimestamp(d.stat().st_ctime)
                })
        return versions

    def cleanup_old_versions(self, keep_count: int = 5):
        """Remove old versions, keeping most recent."""

        versions = self.list_versions()
        if len(versions) <= keep_count:
            return

        for v in versions[:-keep_count]:
            shutil.rmtree(v['path'])
            print(f"Removed old version: {v['tag']}")

# Usage
version_mgr = DataVersionManager('/data/processed')

# After processing new data
version_mgr.create_version('/data/staging/processed', 'clinvar_2026_01_15')

# If issues found, rollback
version_mgr.rollback('clinvar_2026_01_08')
```

---

## 9. Recommended Tools

### 9.1 Download Tools

| Tool | Use Case | Installation |
|------|----------|--------------|
| **wget** | FTP/HTTP downloads, resume support | Built-in |
| **curl** | API downloads, flexible options | Built-in |
| **aria2** | Parallel downloads, BitTorrent | `apt install aria2` |
| **aws cli** | S3 downloads (gnomAD, 1000G) | `pip install awscli` |
| **huggingface-hub** | HF datasets | `pip install huggingface-hub` |

```bash
# wget with resume and rate limiting
wget -c --limit-rate=10M -P /data/downloads/ \
    ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz

# aria2 for parallel download
aria2c -x 16 -s 16 -d /data/downloads/ \
    ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

# AWS S3 (no credentials for public data)
aws s3 cp --no-sign-request \
    s3://gnomad-public-us-east-1/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr1.vcf.bgz \
    /data/downloads/gnomad/
```

### 9.2 Decompression Tools

| Tool | Format | Parallel | Installation |
|------|--------|----------|--------------|
| **pigz** | gzip | Yes | `apt install pigz` |
| **lbzip2** | bzip2 | Yes | `apt install lbzip2` |
| **pixz** | xz | Yes | `apt install pixz` |
| **zstd** | zstd | Yes | `apt install zstd` |

```bash
# Parallel gzip decompression (uses all cores)
pigz -d -k large_file.gz

# Decompress multiple files in parallel
find /data/downloads -name "*.gz" | parallel -j 4 pigz -d -k {}

# Keep compressed, decompress to stdout for streaming
pigz -d -c file.gz | python process.py
```

### 9.3 Progress Monitoring

| Tool | Use Case | Installation |
|------|----------|--------------|
| **pv** | Pipe progress | `apt install pv` |
| **tqdm** | Python progress bars | `pip install tqdm` |
| **htop** | System monitoring | `apt install htop` |

```bash
# Monitor decompression progress
pv large_file.gz | pigz -d > large_file

# Monitor pipeline progress
pv -cN source < input.tsv | python process.py | pv -cN output > output.tsv

# Download with progress
curl -# -o file.gz http://example.com/file.gz
```

### 9.4 GNU Parallel for Batch Processing

```bash
# Install
apt install parallel

# Process files in parallel
find /data/raw -name "*.xml" | \
    parallel -j 8 --progress --eta \
    'python process.py {} > /data/processed/{/.}.json'

# Parallel downloads with progress
cat urls.txt | parallel -j 10 --bar 'wget -q -P /data/downloads/ {}'

# Process with memory limit per job
parallel -j 4 --memfree 2G 'python heavy_process.py {}'
```

### 9.5 Python Libraries Summary

```txt
# requirements-pipeline.txt

# Streaming parsers
ijson>=3.2.0          # JSON streaming
lxml>=5.0.0           # XML streaming
cyvcf2>=0.30.0        # VCF streaming

# Data processing
pandas>=2.0.0         # DataFrames
numpy>=1.24.0         # Arrays
dask>=2024.1.0        # Out-of-core processing
polars>=0.20.0        # Fast DataFrames (alternative)

# Parallel processing
aiohttp>=3.9.0        # Async HTTP
aiofiles>=23.0.0      # Async file I/O

# Database
psycopg2-binary>=2.9.0    # PostgreSQL
pymongo>=4.6.0            # MongoDB

# Utilities
tqdm>=4.66.0          # Progress bars
psutil>=5.9.0         # Memory monitoring
```

---

## 10. Example Pipeline Scripts

### 10.1 Bash: Full Download Workflow

```bash
#!/bin/bash
# download-all-sources.sh
# Full download workflow for all data sources

set -euo pipefail

# Configuration
DATA_DIR="/data"
DOWNLOAD_DIR="$DATA_DIR/downloads"
LOG_DIR="$DATA_DIR/logs"
DATE=$(date +%Y%m%d)

mkdir -p "$DOWNLOAD_DIR"/{ncbi,ebi,community,nih} "$LOG_DIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/download_$DATE.log"
}

download_with_retry() {
    local url=$1
    local dest=$2
    local max_retries=3
    local retry=0

    while [ $retry -lt $max_retries ]; do
        if wget -c -q --show-progress -O "$dest" "$url"; then
            log "Downloaded: $dest"
            return 0
        fi
        retry=$((retry + 1))
        log "Retry $retry/$max_retries for $url"
        sleep 10
    done

    log "FAILED: $url"
    return 1
}

# ============================================
# TIER 1: Core Sources
# ============================================

log "=== Starting Tier 1 Downloads ==="

# ClinVar (weekly release)
log "Downloading ClinVar..."
download_with_retry \
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
    "$DOWNLOAD_DIR/ncbi/clinvar/clinvar_$DATE.vcf.gz"

# dbSNP (common variants only for MVP)
log "Downloading dbSNP common variants..."
download_with_retry \
    "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz" \
    "$DOWNLOAD_DIR/ncbi/dbsnp/common_all_$DATE.vcf.gz"

# SNPedia (via bulk export)
log "Downloading SNPedia..."
python3 << 'EOF'
from scripts.etl.extract.snpedia_extractor import SNPediaExtractor
import json

extractor = SNPediaExtractor()
snps = extractor.get_all_snps()

output_path = f"/data/downloads/community/snpedia/snpedia_{DATE}.json"
with open(output_path, 'w') as f:
    json.dump(snps, f)
print(f"Exported {len(snps)} SNPs")
EOF

# PharmGKB (requires registration - manual step noted)
log "PharmGKB requires manual download from https://www.pharmgkb.org/downloads"

# Reactome
log "Downloading Reactome..."
download_with_retry \
    "https://reactome.org/download/current/ReactomePathways.txt" \
    "$DOWNLOAD_DIR/ebi/reactome/pathways_$DATE.txt"
download_with_retry \
    "https://reactome.org/download/current/ReactomePathwaysRelation.txt" \
    "$DOWNLOAD_DIR/ebi/reactome/pathway_relations_$DATE.txt"

# GWAS Catalog
log "Downloading GWAS Catalog..."
download_with_retry \
    "https://www.ebi.ac.uk/gwas/api/search/downloads/full" \
    "$DOWNLOAD_DIR/ebi/gwas/gwas_catalog_$DATE.tsv"

# ============================================
# TIER 2: Enrichment Sources
# ============================================

log "=== Starting Tier 2 Downloads ==="

# UniProt (human proteome)
log "Downloading UniProt human proteome..."
download_with_retry \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=json&query=%28organism_id%3A9606%29" \
    "$DOWNLOAD_DIR/ebi/uniprot/human_proteome_$DATE.json.gz"

# WikiPathways
log "Downloading WikiPathways..."
download_with_retry \
    "https://data.wikipathways.org/current/gpml/wikipathways-current-gpml-Homo_sapiens.zip" \
    "$DOWNLOAD_DIR/community/wikipathways/wikipathways_$DATE.zip"

# ChEMBL (SQLite database)
log "Downloading ChEMBL..."
CHEMBL_VERSION="33"  # Update as needed
download_with_retry \
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_${CHEMBL_VERSION}_sqlite.tar.gz" \
    "$DOWNLOAD_DIR/ebi/chembl/chembl_${CHEMBL_VERSION}_$DATE.tar.gz"

# NIH ODS Supplement Fact Sheets
log "Downloading NIH ODS data..."
python3 << 'EOF'
import requests
import json

# Fetch all supplement fact sheets via API
base_url = "https://ods.od.nih.gov/api/v1"
response = requests.get(f"{base_url}/factsheets")
data = response.json()

output_path = f"/data/downloads/nih/ods/factsheets_{DATE}.json"
with open(output_path, 'w') as f:
    json.dump(data, f)
print(f"Downloaded {len(data)} fact sheets")
EOF

# ============================================
# TIER 3: Large Datasets (Optional)
# ============================================

if [ "${DOWNLOAD_TIER3:-false}" = "true" ]; then
    log "=== Starting Tier 3 Downloads ==="

    # gnomAD (summary statistics only - full dataset is 30TB)
    log "Downloading gnomAD summary (this will take a while)..."
    aws s3 cp --no-sign-request \
        s3://gnomad-public-us-east-1/release/4.1/vcf/joint/ \
        "$DOWNLOAD_DIR/gnomad/" \
        --recursive \
        --exclude "*" \
        --include "gnomad.joint.v4.1.sites.chr*.vcf.bgz"

    # PubMed baseline (100GB)
    log "Downloading PubMed baseline..."
    mkdir -p "$DOWNLOAD_DIR/ncbi/pubmed/baseline"
    wget -c -r -np -nd -P "$DOWNLOAD_DIR/ncbi/pubmed/baseline" \
        "ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"
fi

# ============================================
# Post-Download Processing
# ============================================

log "=== Decompressing files ==="

# Decompress gzip files in parallel
find "$DOWNLOAD_DIR" -name "*.gz" -type f | \
    parallel -j 4 --progress 'pigz -d -k {}'

# Decompress tar.gz archives
find "$DOWNLOAD_DIR" -name "*.tar.gz" -type f | \
    while read f; do
        tar -xzf "$f" -C "$(dirname "$f")"
    done

# Unzip archives
find "$DOWNLOAD_DIR" -name "*.zip" -type f | \
    while read f; do
        unzip -o "$f" -d "$(dirname "$f")"
    done

log "=== Download complete ==="
log "Total size: $(du -sh "$DOWNLOAD_DIR" | cut -f1)"

# Record state
echo "{
    \"date\": \"$DATE\",
    \"completed\": \"$(date -Iseconds)\",
    \"size\": \"$(du -sb "$DOWNLOAD_DIR" | cut -f1)\"
}" > "$DATA_DIR/state/last_download.json"
```

### 10.2 Python: Complete Processing Pipeline

```python
#!/usr/bin/env python3
"""
process-pipeline.py
Complete data processing pipeline for all sources.
"""

import asyncio
import logging
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import Generator, Iterator
import json
import gzip

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'/data/logs/process_{datetime.now():%Y%m%d}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================
# Configuration
# ============================================

@dataclass
class PipelineConfig:
    data_dir: Path = Path('/data')
    download_dir: Path = Path('/data/downloads')
    raw_dir: Path = Path('/data/raw')
    processed_dir: Path = Path('/data/processed')
    staging_dir: Path = Path('/data/staging')
    db_url: str = 'postgresql://localhost/gene'
    batch_size: int = 10000
    memory_limit_mb: int = 8000

config = PipelineConfig()

# ============================================
# Streaming Parsers
# ============================================

def stream_vcf(filepath: Path) -> Generator[dict, None, None]:
    """Stream VCF file variants."""
    from cyvcf2 import VCF

    vcf = VCF(str(filepath))
    for variant in vcf:
        yield {
            'chrom': variant.CHROM,
            'pos': variant.POS,
            'id': variant.ID,
            'ref': variant.REF,
            'alt': variant.ALT,
            'info': dict(variant.INFO)
        }
    vcf.close()

def stream_xml(filepath: Path, tag: str) -> Generator, None, None]:
    """Stream XML elements."""
    from lxml import etree

    opener = gzip.open if str(filepath).endswith('.gz') else open

    with opener(filepath, 'rb') as f:
        context = etree.iterparse(f, events=('end',), tag=tag)
        for event, elem in context:
            yield elem
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]

def stream_json_array(filepath: Path, prefix: str = 'item') -> Generator[dict, None, None]:
    """Stream JSON array items."""
    import ijson

    opener = gzip.open if str(filepath).endswith('.gz') else open

    with opener(filepath, 'rb') as f:
        for item in ijson.items(f, prefix):
            yield item

# ============================================
# Processors
# ============================================

class SNPProcessor:
    """Process SNP data from multiple sources."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.processed_count = 0

    def process_clinvar(self) -> Generator[dict, None, None]:
        """Process ClinVar VCF."""
        logger.info("Processing ClinVar...")

        clinvar_file = list(self.config.download_dir.glob('ncbi/clinvar/*.vcf'))[0]

        for variant in stream_vcf(clinvar_file):
            if variant['id'] and variant['id'].startswith('rs'):
                yield {
                    'rs_number': variant['id'],
                    'chromosome': variant['chrom'].replace('chr', ''),
                    'position': variant['pos'],
                    'ref_allele': variant['ref'],
                    'alt_alleles': ','.join(variant['alt']),
                    'clinical_significance': variant['info'].get('CLNSIG', ''),
                    'source': 'clinvar'
                }
                self.processed_count += 1

                if self.processed_count % 100000 == 0:
                    logger.info(f"Processed {self.processed_count:,} ClinVar variants")

    def process_dbsnp(self) -> Generator[dict, None, None]:
        """Process dbSNP VCF."""
        logger.info("Processing dbSNP...")

        dbsnp_file = list(self.config.download_dir.glob('ncbi/dbsnp/*.vcf'))[0]

        for variant in stream_vcf(dbsnp_file):
            if variant['id']:
                yield {
                    'rs_number': variant['id'],
                    'chromosome': variant['chrom'].replace('chr', ''),
                    'position': variant['pos'],
                    'ref_allele': variant['ref'],
                    'alt_alleles': ','.join(variant['alt']),
                    'global_maf': variant['info'].get('AF', None),
                    'source': 'dbsnp'
                }

    def merge_sources(self, *generators) -> Generator[dict, None, None]:
        """Merge SNPs from multiple sources."""
        from collections import defaultdict

        snp_data = defaultdict(dict)

        for gen in generators:
            for snp in gen:
                rs = snp['rs_number']
                # Merge fields, preferring non-null values
                for key, value in snp.items():
                    if value and key != 'source':
                        snp_data[rs][key] = value

                # Track sources
                if 'sources' not in snp_data[rs]:
                    snp_data[rs]['sources'] = []
                snp_data[rs]['sources'].append(snp['source'])

        for rs, data in snp_data.items():
            yield {'rs_number': rs, **data}


class GeneProcessor:
    """Process gene data."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def process_uniprot(self) -> Generator[dict, None, None]:
        """Process UniProt human proteome."""
        logger.info("Processing UniProt...")

        uniprot_file = list(self.config.download_dir.glob('ebi/uniprot/*.json.gz'))[0]

        for protein in stream_json_array(uniprot_file, 'results.item'):
            gene_names = protein.get('genes', [])
            primary_gene = gene_names[0]['geneName']['value'] if gene_names else None

            yield {
                'uniprot_id': protein['primaryAccession'],
                'gene_symbol': primary_gene,
                'protein_name': protein.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
                'function': self._extract_function(protein),
                'source': 'uniprot'
            }

    def _extract_function(self, protein: dict) -> str:
        """Extract function description from UniProt entry."""
        comments = protein.get('comments', [])
        for comment in comments:
            if comment.get('commentType') == 'FUNCTION':
                texts = comment.get('texts', [])
                if texts:
                    return texts[0].get('value', '')
        return ''


class PathwayProcessor:
    """Process pathway data."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def process_reactome(self) -> Generator[dict, None, None]:
        """Process Reactome pathways."""
        logger.info("Processing Reactome...")

        pathways_file = list(self.config.download_dir.glob('ebi/reactome/pathways*.txt'))[0]

        with open(pathways_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[2] == 'Homo sapiens':
                    yield {
                        'pathway_id': parts[0],
                        'name': parts[1],
                        'species': parts[2],
                        'source': 'reactome'
                    }


# ============================================
# Database Loader
# ============================================

class DatabaseLoader:
    """Load processed data to database."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.conn = None

    def connect(self):
        import psycopg2
        self.conn = psycopg2.connect(self.config.db_url)

    def close(self):
        if self.conn:
            self.conn.close()

    def bulk_load_snps(self, snps: Iterator[dict]):
        """Bulk load SNPs using COPY."""
        from io import StringIO
        import csv

        cursor = self.conn.cursor()
        columns = ['rs_number', 'chromosome', 'position', 'ref_allele',
                   'alt_alleles', 'clinical_significance', 'global_maf']

        buffer = StringIO()
        writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction='ignore')

        batch = []
        total = 0

        for snp in snps:
            batch.append(snp)

            if len(batch) >= self.config.batch_size:
                buffer.seek(0)
                buffer.truncate()

                for row in batch:
                    writer.writerow(row)

                buffer.seek(0)
                cursor.copy_expert(
                    f"COPY snps ({','.join(columns)}) FROM STDIN WITH CSV",
                    buffer
                )
                self.conn.commit()

                total += len(batch)
                logger.info(f"Loaded {total:,} SNPs")
                batch = []

        # Load remaining
        if batch:
            buffer.seek(0)
            buffer.truncate()
            for row in batch:
                writer.writerow(row)
            buffer.seek(0)
            cursor.copy_expert(
                f"COPY snps ({','.join(columns)}) FROM STDIN WITH CSV",
                buffer
            )
            self.conn.commit()
            total += len(batch)

        cursor.close()
        logger.info(f"Total SNPs loaded: {total:,}")


# ============================================
# Main Pipeline
# ============================================

class DataPipeline:
    """Main pipeline orchestrator."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.snp_processor = SNPProcessor(config)
        self.gene_processor = GeneProcessor(config)
        self.pathway_processor = PathwayProcessor(config)
        self.loader = DatabaseLoader(config)

    def run(self):
        """Execute full pipeline."""
        start_time = datetime.now()
        logger.info("=== Starting Data Processing Pipeline ===")

        try:
            # Phase 1: Process SNPs
            logger.info("Phase 1: Processing SNPs...")
            clinvar_snps = self.snp_processor.process_clinvar()
            dbsnp_snps = self.snp_processor.process_dbsnp()
            merged_snps = self.snp_processor.merge_sources(clinvar_snps, dbsnp_snps)

            # Save processed SNPs
            self._save_processed(merged_snps, 'snps.jsonl')

            # Phase 2: Process Genes
            logger.info("Phase 2: Processing Genes...")
            genes = self.gene_processor.process_uniprot()
            self._save_processed(genes, 'genes.jsonl')

            # Phase 3: Process Pathways
            logger.info("Phase 3: Processing Pathways...")
            pathways = self.pathway_processor.process_reactome()
            self._save_processed(pathways, 'pathways.jsonl')

            # Phase 4: Load to Database
            logger.info("Phase 4: Loading to Database...")
            self.loader.connect()

            with open(self.config.processed_dir / 'snps.jsonl') as f:
                snps = (json.loads(line) for line in f)
                self.loader.bulk_load_snps(snps)

            self.loader.close()

            # Phase 5: Generate Embeddings
            logger.info("Phase 5: Generating Embeddings...")
            self._generate_embeddings()

            elapsed = datetime.now() - start_time
            logger.info(f"=== Pipeline Complete in {elapsed} ===")

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise

    def _save_processed(self, data: Iterator[dict], filename: str):
        """Save processed data to JSONL file."""
        output_path = self.config.processed_dir / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)

        count = 0
        with open(output_path, 'w') as f:
            for item in data:
                f.write(json.dumps(item) + '\n')
                count += 1

                if count % 100000 == 0:
                    logger.info(f"Saved {count:,} records to {filename}")

        logger.info(f"Saved {count:,} total records to {filename}")

    def _generate_embeddings(self):
        """Generate embeddings for all entities."""
        from sentence_transformers import SentenceTransformer
        import numpy as np

        model = SentenceTransformer('sentence-transformers/all-MiniLM-L6-v2')

        # Generate SNP embeddings
        logger.info("Generating SNP embeddings...")

        texts = []
        ids = []

        with open(self.config.processed_dir / 'snps.jsonl') as f:
            for line in f:
                snp = json.loads(line)
                text = f"SNP {snp['rs_number']} "
                if snp.get('clinical_significance'):
                    text += f"clinical significance: {snp['clinical_significance']} "
                texts.append(text)
                ids.append(snp['rs_number'])

                # Batch embedding generation
                if len(texts) >= 10000:
                    embeddings = model.encode(texts, show_progress_bar=True)
                    self._save_embeddings(ids, embeddings, 'snp_embeddings.npy')
                    texts = []
                    ids = []

        # Final batch
        if texts:
            embeddings = model.encode(texts, show_progress_bar=True)
            self._save_embeddings(ids, embeddings, 'snp_embeddings.npy')

    def _save_embeddings(self, ids: list, embeddings, filename: str):
        """Append embeddings to file."""
        import numpy as np

        output_path = self.config.processed_dir / 'embeddings' / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save incrementally
        if output_path.exists():
            existing = np.load(output_path)
            embeddings = np.vstack([existing, embeddings])

        np.save(output_path, embeddings)

        # Save ID mapping
        id_path = output_path.with_suffix('.ids.json')
        if id_path.exists():
            with open(id_path) as f:
                existing_ids = json.load(f)
            ids = existing_ids + ids

        with open(id_path, 'w') as f:
            json.dump(ids, f)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Data Processing Pipeline')
    parser.add_argument('--data-dir', default='/data', help='Base data directory')
    parser.add_argument('--db-url', default='postgresql://localhost/gene', help='Database URL')
    parser.add_argument('--batch-size', type=int, default=10000, help='Batch size for loading')
    args = parser.parse_args()

    config = PipelineConfig(
        data_dir=Path(args.data_dir),
        download_dir=Path(args.data_dir) / 'downloads',
        raw_dir=Path(args.data_dir) / 'raw',
        processed_dir=Path(args.data_dir) / 'processed',
        staging_dir=Path(args.data_dir) / 'staging',
        db_url=args.db_url,
        batch_size=args.batch_size
    )

    pipeline = DataPipeline(config)
    pipeline.run()
```

---

## Quick Reference

### Memory Guidelines

| File Size | Chunk Size | Expected RAM | Processing Time |
|-----------|------------|--------------|-----------------|
| 1 GB | 10,000 rows | ~500 MB | ~5 min |
| 10 GB | 50,000 rows | ~1 GB | ~30 min |
| 100 GB | 100,000 rows | ~2 GB | ~5 hours |
| 1 TB | 100,000 rows | ~4 GB | ~50 hours |

### Command Cheatsheet

```bash
# Download with resume
wget -c URL

# Parallel decompression
pigz -d -k *.gz

# Monitor progress
pv file.gz | pigz -d > file

# Parallel file processing
parallel -j 8 'python process.py {}' ::: *.json

# PostgreSQL bulk load
psql -c "\COPY table FROM 'file.csv' CSV HEADER"

# Check memory usage
ps aux --sort=-%mem | head

# Estimate file line count without reading whole file
wc -l < file.txt
```

### Processing Checklist

- [ ] Verify disk space: `df -h`
- [ ] Check download integrity: `md5sum -c checksums.txt`
- [ ] Monitor memory during processing: `htop`
- [ ] Validate row counts after loading
- [ ] Create database indexes after bulk load
- [ ] Run data quality checks
- [ ] Update state tracking
- [ ] Clean up temporary files

---

## References

- [PostgreSQL COPY Documentation](https://www.postgresql.org/docs/current/sql-copy.html)
- [Dask Best Practices](https://docs.dask.org/en/stable/best-practices.html)
- [GNU Parallel Tutorial](https://www.gnu.org/software/parallel/parallel_tutorial.html)
- [ijson Documentation](https://github.com/ICRAR/ijson)
- [cyvcf2 Documentation](https://brentp.github.io/cyvcf2/)
