---
id: integration-size-estimates
title: "Database Size Estimates & Download Planning"
type: integration
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [integration, cross-references, apis]
---

**Parent:** [../_index.md](../_index.md)

# Database Size Estimates & Download Planning

**Created:** January 2026
**Purpose:** Estimate data volumes for capacity planning and download strategy
**Status:** Active reference

---

## Executive Summary

This document provides size estimates for all major traditional medicine, genetic, and pharmaceutical databases to help plan:
- Storage requirements
- Download bandwidth
- Processing time
- Database schema design
- ETL pipeline architecture

**Key Finding:** Most database bulks downloads range from 100 MB to 50 GB. Full integration of all sources requires ~200-500 GB storage.

---

## Part 1: Traditional Medicine Databases

### 1.1 TCM Databases

#### BATMAN-TCM 2.0

**URL:** http://bionet.ncpsb.org.cn/batman-tcm/

| File | Records | Estimated Size | Format |
|------|---------|----------------|--------|
| formula_ingredient.txt | 54,832 formulas | 50 MB | TSV |
| ingredient_target_known.txt | 17,068 interactions | 10 MB | TSV |
| ingredient_target_predicted.txt | 2,319,272 interactions | 500 MB | TSV |
| **Total** | **2.4M+ records** | **~560 MB** | - |

**Download Time:**
- Fast connection (100 Mbps): ~1 minute
- Medium connection (10 Mbps): ~8 minutes

**Storage Recommendation:** 2 GB (raw + indexed)

---

#### TCMBank

**URL:** https://tcmbank.cn/

| Entity | Count | Estimated Size | Format |
|--------|-------|----------------|--------|
| Herbs | 9,192 | 5 MB | JSON/CSV |
| Ingredients | 61,966 | 50 MB | JSON/CSV |
| Targets | 15,179 | 10 MB | JSON/CSV |
| Diseases | 2,710 | 2 MB | JSON/CSV |
| **Total** | **~89,000 entities** | **~67 MB** | - |

**Download Time:** <1 minute (most connections)

**Storage Recommendation:** 200 MB (raw + indexed)

---

#### ETCM (Encyclopedia of TCM)

**URL:** http://www.tcmip.cn/ETCM/

| Entity | Count | Estimated Size | Format |
|--------|-------|----------------|--------|
| Herbs | 403 | 2 MB | Web scraping |
| Ingredients | 7,274 | 30 MB | Web scraping |
| Targets | 2,266 | 5 MB | Web scraping |
| Pathways | 3,260 | 10 MB | Web scraping |
| **Total** | **~13,000 entities** | **~47 MB** | - |

**Note:** No bulk download; requires web scraping (rate-limited)

**Download Time:** 2-4 hours (with rate limiting)

**Storage Recommendation:** 150 MB

---

### 1.2 Ayurveda Databases

#### IMPPAT 2.0

**URL:** https://cb.imsc.res.in/imppat/

| Entity | Count | Estimated Size | Format |
|--------|-------|----------------|--------|
| Medicinal plants | 4,010 | 5 MB | TSV export |
| Phytochemicals | 17,967 | 100 MB | TSV + SDF |
| Predicted targets | 27,365 pairs | 50 MB | TSV |
| Structure files | 17,967 files | 500 MB | SDF/MOL |
| **Total** | **~67,000 records** | **~655 MB** | - |

**Download Time:** 5-10 minutes (fast connection)

**Storage Recommendation:** 2 GB (includes structure files)

---

#### NPACT

**URL:** https://webs.iiitd.edu.in/raghava/npact/

| Entity | Count | Estimated Size | Format |
|--------|-------|----------------|--------|
| Natural products | 1,574 | 10 MB | Web export |
| Cancer targets | ~500 | 2 MB | Web export |
| **Total** | **~2,000 records** | **~12 MB** | - |

**Download Time:** <1 minute

**Storage Recommendation:** 50 MB

---

### 1.3 Kampo Databases

#### KampoDB

**URL:** https://wakanmoview.inm.u-toyama.ac.jp/kampo/

| Entity | Count | Estimated Size | Format |
|--------|-------|----------------|--------|
| Kampo formulas | 298 | 2 MB | Web queries |
| Crude drugs | 180 | 1 MB | Web queries |
| Natural compounds | 3,002 | 20 MB | Web queries |
| Predicted targets | 62,906 proteins | 100 MB | Docking data |
| **Total** | **~66,000 records** | **~123 MB** | - |

**Download Time:** 1-2 hours (web scraping)

**Storage Recommendation:** 500 MB

---

### 1.4 Western Herbal Databases

#### Dr. Duke's Phytochemical Database

**URL:** https://phytochem.nal.usda.gov

| File | Records | Estimated Size | Format |
|------|---------|----------------|--------|
| plants.csv | ~1,000 plants | 1 MB | CSV |
| chemicals.csv | ~80,000 compounds | 50 MB | CSV |
| activities.csv | ~300 activities | 2 MB | CSV |
| plant_chemicals.csv | ~200,000 links | 100 MB | CSV |
| chemical_activities.csv | ~100,000 links | 50 MB | CSV |
| **Total** | **~381,000 records** | **~203 MB** | CSV |

**Download Time:** 2-3 minutes

**Storage Recommendation:** 1 GB (raw + indexed)

---

#### DSLD (NIH Dietary Supplements)

**URL:** https://dsld.od.nih.gov

| Entity | Count | Estimated Size | API Access |
|--------|-------|----------------|------------|
| Product labels | 200,000+ | N/A (API only) | REST API |
| Ingredients | ~50,000 | N/A | REST API |
| Brands | ~10,000 | N/A | REST API |

**Note:** API-only (no bulk download)

**API Rate Limit:** ~1000 requests/hour

**Storage Recommendation:** 5 GB (if caching API responses)

---

## Part 2: Genetic & Pathway Databases

### 2.1 Gene/Protein Databases

#### UniProt

**URL:** https://www.uniprot.org/

| File | Records | Compressed Size | Uncompressed Size | Format |
|------|---------|-----------------|-------------------|--------|
| UniProtKB/Swiss-Prot | ~570,000 reviewed | 250 MB | 1.5 GB | XML/FASTA |
| UniProtKB/TrEMBL | ~245 million | 50 GB | 200 GB | XML/FASTA |
| Human proteins only | ~20,000 | 50 MB | 200 MB | XML/FASTA |

**Recommendation:** Download human proteins only (50 MB compressed)

**Download Time:** 1 minute (human only)

**Storage Recommendation:** 500 MB (human proteins + indexes)

---

#### NCBI Gene

**URL:** https://ftp.ncbi.nlm.nih.gov/gene/

| File | Records | Compressed Size | Format |
|------|---------|-----------------|--------|
| gene_info (all species) | ~60 million genes | 2 GB | TSV |
| gene_info (human only) | ~60,000 genes | 20 MB | TSV |
| gene2pubmed | ~16 million links | 500 MB | TSV |
| gene2go | ~600,000 links | 50 MB | TSV |

**Recommendation:** Filter to human genes

**Download Time:** 5 minutes (human subset)

**Storage Recommendation:** 200 MB (human genes)

---

### 2.2 Pathway Databases

#### Reactome

**URL:** https://reactome.org/download-data

| File | Records | Compressed Size | Format |
|------|---------|-----------------|--------|
| All pathways (BioPAX) | ~2,500 pathways | 100 MB | XML |
| Protein-pathway links | ~11,000 proteins | 20 MB | TSV |
| Reactions | ~13,000 reactions | 30 MB | TSV |
| **Total** | **~26,500 entities** | **~150 MB** | - |

**Download Time:** 2-3 minutes

**Storage Recommendation:** 1 GB (raw + graph database)

---

#### KEGG

**URL:** https://www.kegg.jp/kegg/download/

| Entity | Count | Estimated Size | Access |
|--------|-------|----------------|--------|
| Pathways (human) | ~500 pathways | 50 MB | API only |
| Genes (human) | ~20,000 genes | 20 MB | API only |
| Compounds | ~18,000 compounds | 30 MB | API only |

**Note:** KEGG FTP requires paid license. Use API for academic use.

**API Rate Limit:** No official limit, but respect 1 request/second

**Storage Recommendation:** 200 MB (API cache)

---

#### WikiPathways

**URL:** https://www.wikipathways.org/

| File | Records | Compressed Size | Format |
|------|---------|-----------------|--------|
| All pathways (GPML) | ~3,000 pathways | 200 MB | XML |
| Human pathways only | ~800 pathways | 50 MB | XML |
| GMT gene sets | ~800 pathways | 5 MB | GMT |

**Download Time:** 3-5 minutes

**Storage Recommendation:** 500 MB

---

### 2.3 Chemical Databases

#### PubChem

**URL:** https://ftp.ncbi.nlm.nih.gov/pubchem/

| File | Records | Compressed Size | Uncompressed Size | Format |
|------|---------|-----------------|-------------------|--------|
| Compound (all) | ~110 million | 100 GB | 500 GB | SDF |
| Substance (all) | ~280 million | 200 GB | 1 TB | SDF |
| BioAssay | ~1.5 million | 50 GB | 200 GB | XML |

**Recommendation:** Use PubChem API for on-demand queries (don't download full dumps)

**API Rate Limit:** 5 requests/second (no API key), 10 requests/second (with key)

**Storage Recommendation:** 10 GB (API cache for relevant compounds only)

---

#### ChEMBL

**URL:** https://ftp.ebi.ac.uk/pub/databases/chembl/

| File | Records | Compressed Size | Format |
|------|---------|-----------------|--------|
| chembl_32_postgresql.tar.gz | 2.3 million compounds | 4 GB | PostgreSQL dump |
| chembl_32_mysql.tar.gz | 2.3 million compounds | 5 GB | MySQL dump |
| chembl_32.sdf.gz | 2.3 million compounds | 3 GB | SDF |

**Download Time:** 30-60 minutes (depending on connection)

**Storage Recommendation:** 20 GB (database + indexes)

---

## Part 3: Cross-Reference & Integration Databases

### 3.1 Identifier Mapping

#### UniProt ID Mapping

**URL:** https://www.uniprot.org/id-mapping

| Entity | Count | Estimated Size | Access |
|--------|-------|----------------|--------|
| UniProt ↔ Ensembl | ~20,000 (human) | 5 MB | API/Web |
| UniProt ↔ HGNC | ~20,000 (human) | 2 MB | API/Web |
| UniProt ↔ PDB | ~15,000 | 3 MB | API/Web |

**Note:** API-only (no bulk download)

**Storage Recommendation:** 50 MB (API cache)

---

#### Gene Ontology

**URL:** http://geneontology.org/docs/download-ontology/

| File | Records | Compressed Size | Format |
|------|---------|-----------------|--------|
| go-basic.obo | ~45,000 terms | 10 MB | OBO |
| goa_human.gaf | ~600,000 annotations | 50 MB | GAF |

**Download Time:** 1-2 minutes

**Storage Recommendation:** 200 MB

---

### 3.2 Wikidata

**URL:** https://dumps.wikimedia.org/wikidatawiki/entities/

| File | Records | Compressed Size | Uncompressed Size | Format |
|------|---------|-----------------|-------------------|--------|
| latest-all.json.gz | ~100 million items | 100 GB | 1 TB | JSON |
| latest-all.json.bz2 | ~100 million items | 60 GB | 1 TB | JSON |

**Recommendation:** Use SPARQL API instead of downloading full dump

**SPARQL Rate Limit:** 60-second timeout per query

**Storage Recommendation:** 1 GB (SPARQL result cache)

---

## Part 4: Aggregate Size Estimates

### 4.1 Minimal Integration (MVP)

**Databases:**
- BATMAN-TCM 2.0 (predicted targets)
- IMPPAT 2.0 (Ayurveda)
- UniProt (human proteins)
- Reactome (pathways)

| Category | Raw Data | Indexed/Processed | Total |
|----------|----------|-------------------|-------|
| Traditional Medicine | 1.2 GB | 2.5 GB | 3.7 GB |
| Proteins/Genes | 200 MB | 500 MB | 700 MB |
| Pathways | 150 MB | 1 GB | 1.15 GB |
| **Total** | **1.55 GB** | **4 GB** | **~5.5 GB** |

**Download Time:** 15-30 minutes

**Processing Time:** 2-4 hours (parsing + indexing)

---

### 4.2 Comprehensive Integration

**Databases:**
- All TCM databases (BATMAN-TCM, TCMBank, ETCM)
- All Ayurveda databases (IMPPAT, NPACT)
- KampoDB (Kampo)
- Dr. Duke's (Western herbal)
- ChEMBL (compound-target)
- UniProt + NCBI Gene
- Reactome + WikiPathways + KEGG
- PubChem (API cache)

| Category | Raw Data | Indexed/Processed | Total |
|----------|----------|-------------------|-------|
| Traditional Medicine | 3 GB | 8 GB | 11 GB |
| Proteins/Genes | 500 MB | 2 GB | 2.5 GB |
| Pathways | 1 GB | 3 GB | 4 GB |
| Compounds (ChEMBL) | 5 GB | 15 GB | 20 GB |
| PubChem Cache | 5 GB | 5 GB | 10 GB |
| **Total** | **14.5 GB** | **33 GB** | **~47.5 GB** |

**Download Time:** 2-4 hours

**Processing Time:** 24-48 hours (full ETL pipeline)

---

### 4.3 Full Integration (Research-Grade)

**Adds:**
- Full ChEMBL database
- PubChem subset (natural products)
- BindingDB
- Full WikiPathways
- Clinical trials data
- DrugBank

| Category | Raw Data | Indexed/Processed | Total |
|----------|----------|-------------------|-------|
| Traditional Medicine | 5 GB | 15 GB | 20 GB |
| Proteins/Genes | 2 GB | 5 GB | 7 GB |
| Pathways | 2 GB | 8 GB | 10 GB |
| Compounds | 50 GB | 100 GB | 150 GB |
| Clinical Data | 10 GB | 20 GB | 30 GB |
| **Total** | **69 GB** | **148 GB** | **~217 GB** |

**Download Time:** 1-2 days

**Processing Time:** 3-7 days (distributed processing recommended)

---

## Part 5: Download Strategy

### 5.1 Parallel Downloads

**Recommended approach:**

```bash
# Create download directory
mkdir -p data/{raw,processed}

# Parallel download using xargs or GNU parallel
cat << EOF | parallel -j 4
wget -P data/raw/ http://bionet.ncpsb.org.cn/batman-tcm/download/formula_ingredient.txt
wget -P data/raw/ https://cb.imsc.res.in/imppat/export/compounds.tsv
wget -P data/raw/ https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_human.dat.gz
wget -P data/raw/ https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
EOF
```

**Estimated parallel download time:**
- 4 parallel connections: 25% of sequential time
- 8 parallel connections: 15% of sequential time

---

### 5.2 Incremental Updates

**Strategy:**
1. Initial bulk download (one-time)
2. Periodic incremental updates (weekly/monthly)

**Update sizes (typical):**
- BATMAN-TCM: ~50 MB/month (new predictions)
- IMPPAT: ~10 MB/quarter (new plants)
- UniProt: ~100 MB/month (new annotations)
- Reactome: ~50 MB/quarter (new pathways)

**Total incremental:** ~200 MB/month

---

### 5.3 Caching Strategy

**For API-based sources (PubChem, KEGG, DSLD):**

```python
import requests_cache

# Create SQLite cache with 30-day expiration
session = requests_cache.CachedSession(
    'api_cache',
    expire_after=2592000  # 30 days in seconds
)

# Cached request
response = session.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/969516/JSON')
```

**Estimated cache size:**
- PubChem: 5-10 GB (for ~100,000 compounds)
- KEGG: 200 MB (for human pathways/genes)
- DSLD: 3 GB (for ~200,000 product labels)

---

## Part 6: Database Schema Size Estimates

### 6.1 PostgreSQL Schema Sizes

**Estimated sizes for normalized relational schema:**

| Table | Records | Table Size | Index Size | Total |
|-------|---------|------------|------------|-------|
| compounds | 100,000 | 50 MB | 30 MB | 80 MB |
| plants | 10,000 | 5 MB | 3 MB | 8 MB |
| proteins | 20,000 | 10 MB | 8 MB | 18 MB |
| pathways | 5,000 | 20 MB | 10 MB | 30 MB |
| compound_targets | 2.5M | 500 MB | 300 MB | 800 MB |
| compound_plants | 500,000 | 100 MB | 60 MB | 160 MB |
| protein_pathways | 50,000 | 20 MB | 15 MB | 35 MB |
| **Total** | **~3.2M records** | **705 MB** | **426 MB** | **~1.13 GB** |

**With full-text search indexes:** Add 50% overhead = ~1.7 GB

**With caching/materialized views:** Add 30% overhead = ~2.2 GB

**Total PostgreSQL database:** ~5 GB (includes WAL, temp files)

---

### 6.2 Graph Database Size Estimates (Neo4j)

**Estimated sizes for graph schema:**

| Node Type | Count | Store Size | Index Size | Total |
|-----------|-------|------------|------------|-------|
| Compound | 100,000 | 100 MB | 50 MB | 150 MB |
| Plant | 10,000 | 20 MB | 10 MB | 30 MB |
| Protein | 20,000 | 30 MB | 20 MB | 50 MB |
| Pathway | 5,000 | 50 MB | 30 MB | 80 MB |

| Relationship Type | Count | Store Size |
|-------------------|-------|------------|
| TARGETS | 2.5M | 2 GB |
| FOUND_IN | 500,000 | 400 MB |
| PARTICIPATES_IN | 50,000 | 50 MB |

**Total Neo4j database:** ~8 GB (includes transaction logs)

---

## Part 7: Processing Time Estimates

### 7.1 ETL Pipeline Times

**Hardware assumption:** 8-core CPU, 32 GB RAM, SSD storage

| Stage | Records | Time (Parallel) | Time (Sequential) |
|-------|---------|-----------------|-------------------|
| Download | All | 2 hours | 8 hours |
| Parse CSV/TSV | 3M rows | 30 minutes | 2 hours |
| Parse SDF | 100,000 compounds | 1 hour | 4 hours |
| Normalize IDs | 2.5M pairs | 2 hours | 10 hours |
| Load database | 3M records | 1 hour | 5 hours |
| Build indexes | - | 1 hour | 3 hours |
| **Total** | - | **~7.5 hours** | **~32 hours** |

---

### 7.2 Query Performance Estimates

**After indexing (PostgreSQL):**

| Query Type | Records Returned | Query Time |
|------------|------------------|------------|
| Compound by name | 1 | <10 ms |
| Compound by ID | 1 | <5 ms |
| Targets for compound | 10-100 | 20-50 ms |
| Pathways for gene | 5-20 | 30-100 ms |
| Compounds in plant | 10-50 | 20-80 ms |
| Full-text search | 100 | 100-500 ms |
| Complex join (compound → target → pathway) | 100 | 500-2000 ms |

**After caching (Redis):**
- Cached queries: <5 ms
- Cache hit rate: 70-90% (typical)

---

## Part 8: Recommendations

### 8.1 For MVP (Small Project)

**Budget:**
- Storage: 10 GB
- Download time: 1 hour
- Processing time: 4 hours

**Databases:**
1. BATMAN-TCM (predicted targets)
2. IMPPAT (Ayurveda)
3. UniProt (human proteins)
4. Reactome (pathways)

---

### 8.2 For Production (Medium Project)

**Budget:**
- Storage: 100 GB
- Download time: 4 hours
- Processing time: 24 hours

**Databases:**
- All TCM databases
- All Ayurveda databases
- KampoDB + Dr. Duke's
- ChEMBL (compound-target)
- Reactome + WikiPathways
- PubChem API cache

---

### 8.3 For Research (Large Project)

**Budget:**
- Storage: 500 GB
- Download time: 2 days
- Processing time: 1 week

**Databases:**
- All of the above
- Full ChEMBL database
- BindingDB
- DrugBank
- Clinical trials data
- PubChem subset

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Compressed Size | File size after applying compression (gzip, bzip2) | 500 MB gzipped |
| Uncompressed Size | Original file size before compression | 5 GB raw |
| Working Space | Temporary disk space needed during data processing | 2x download size |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| MVP | Minimum Viable Product - smallest usable version | ~5.5 GB storage |
| Bulk Download | Retrieving entire datasets rather than individual queries | FTP downloads |
| Processing Pipeline | Series of steps to transform raw data into usable format | ETL workflow |
| Index Size | Storage required for database indexes | Search optimization |
| Embedding | Vector representation of text for semantic search | 384-768 dimensions |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ETL | Extract, Transform, Load | Data pipeline |
| FTP | File Transfer Protocol | Download method |
| GB | Gigabyte | Storage unit (10^9 bytes) |
| MB | Megabyte | Storage unit (10^6 bytes) |
| MVP | Minimum Viable Product | Development phase |
| SSD | Solid State Drive | Fast storage |
| TB | Terabyte | Storage unit (10^12 bytes) |

---

*Document compiled January 2026. For download scripts and ETL pipelines, see [integration-guide.md](./integration-guide.md)*
