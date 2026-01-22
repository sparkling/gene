---
id: downloads-processing-pipeline
title: "Data Processing Pipeline for Bulk Downloads"
type: download-guide
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [downloads, bulk-data, pipeline, processing]
---

**Parent:** [Download Guides](./_index.md)

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

[Rest of the original content continues...]

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

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| ETL Pipeline | Extract, Transform, Load workflow for data processing | Download -> Parse -> Database |
| Streaming Parser | Parser that processes data incrementally without loading entire file | ijson, lxml iterparse |
| Bulk Load | Inserting large amounts of data efficiently into a database | PostgreSQL COPY |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| pigz | Parallel gzip compression/decompression tool | Multi-threaded |
| lbzip2 | Parallel bzip2 compression/decompression tool | Multi-threaded |
| ijson | Iterative JSON parser for Python | Streaming JSON |
| lxml | Python XML processing library | Streaming XML |
| cyvcf2 | Fast VCF parser for Python | Variant data |
| HNSW | Hierarchical Navigable Small World index | Vector search |
| Dask | Parallel computing library for Python | Big data processing |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CPU | Central Processing Unit | Processing power |
| CSV | Comma-Separated Values | Data format |
| ETL | Extract, Transform, Load | Data pipeline |
| FTP | File Transfer Protocol | Download method |
| GB | Gigabyte | Storage unit |
| HDD | Hard Disk Drive | Slower storage |
| HNSW | Hierarchical Navigable Small World | Vector index |
| I/O | Input/Output | Data transfer |
| JSON | JavaScript Object Notation | Data format |
| MB | Megabyte | Storage unit |
| NVMe | Non-Volatile Memory Express | Fast SSD |
| RAM | Random Access Memory | Working memory |
| SSD | Solid State Drive | Fast storage |
| TB | Terabyte | Storage unit |
| TSV | Tab-Separated Values | Data format |
| VCF | Variant Call Format | Genetic variants |
| XML | Extensible Markup Language | Data format |
