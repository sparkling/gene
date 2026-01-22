---
id: literature-pipeline-design
title: Research Papers ETL Pipeline Design
category: literature
tier: 2
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [literature, etl, pipeline, pubmed, processing, embeddings]
---

# Research Papers ETL Pipeline Design

**Last Updated:** January 2026
**Target:** 1-10M genetics-relevant papers from PubMed/PMC
**Vector Database:** RuVector (384d embeddings, tiered compression)
**Relational Database:** Supabase PostgreSQL
**Budget Constraint:** Bootstrapped startup, cost-sensitive
**Parent:** [../_index.md](../_index.md)

---

## Executive Summary

| Metric | Value |
|--------|-------|
| **Source** | PubMed/MEDLINE + PMC Open Access |
| **Total Papers in PubMed** | 36M+ citations |
| **Genetics-Relevant Estimate** | 3-5M papers |
| **Filtered Target** | 1-3M high-quality papers |
| **Storage (RuVector)** | ~400MB-1.2GB (with tiered compression) |
| **Processing Time** | 2-4 days (initial), 1-2 hours (daily updates) |
| **Annual Cost** | ~$0 (all open access, self-hosted processing) |

---

## 1. Download Strategy

### 1.1 PubMed Data Architecture

```
PubMed Data Distribution
├── Baseline Files (Annual)
│   ├── ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/
│   ├── 1,219 compressed XML files (~100 GB total)
│   ├── Each file: ~30,000 citations
│   └── Released annually (December)
│
├── Daily Update Files
│   ├── ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/
│   ├── New citations, corrections, retractions
│   └── ~1,000-2,000 new papers/day
│
└── PMC Open Access Subset
    ├── ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/
    ├── Full-text articles under CC licenses
    └── ~4M articles with full text
```

### 1.2 Download Implementation

**Key Components:**
- FTP connection handling
- Checksum verification
- Rate limiting (2-10 requests/second)
- Incremental download tracking
- Resume capability for interrupted downloads

**Rate Limiting:**
- Without API key: 3 requests/second
- With API key: 10 requests/second
- Conservative: 2 requests/second to be respectful

---

## 2. Filtering Pipeline

### 2.1 MeSH Terms for Genetics Filtering

**Primary MeSH Descriptors:**
- Genetics
- Genetic Variation
- Polymorphism, Single Nucleotide
- Polymorphism, Genetic
- Genotype, Phenotype, Alleles
- Gene Expression, Gene Expression Regulation
- Genome-Wide Association Study
- Pharmacogenetics, Pharmacogenomic Variants
- Genetic Predisposition to Disease

**Secondary MeSH Terms:**
- Genes, Proteins, Enzymes, Receptors
- Neoplasms, Cardiovascular Diseases, Metabolic Diseases
- Signal Transduction, Metabolic Networks

**Genetic Subheadings:**
- /genetics, /metabolism, /pharmacology, /drug effects

### 2.2 Keyword Filtering for SNPs, Genes, Pathways

**SNP Patterns:**
- rs\d{1,12} (rs numbers)
- [ACGT]>[ACGT] notation
- p.[A-Z][a-z]{2}\d+[A-Z][a-z]{2} (HGVS protein)
- c.\d+[ACGT]>[ACGT] (HGVS coding)

**Genetics Keywords:**
- polymorphism, variant, mutation, allele, genotype
- haplotype, diplotype, heterozygous, homozygous
- genome-wide, gwas, association study, genetic association
- pharmacogenetic, pharmacogenomic, drug metabolism
- gene expression, transcription, translation, methylation

**Estimated Yield:**
- PubMed Total: ~36M papers
- MeSH Filter: ~8M papers (22%)
- SNP/Gene Keywords: ~5M papers (14%)
- Combined Filter: ~10M papers (28%)
- Confidence >= 0.5: ~3M papers (8%)
- **Target: 1-3M highest relevance papers**

---

## 3. Processing Pipeline

### 3.1 XML Parsing

**Parser:** lxml (streaming) - 10x faster than BeautifulSoup

**Key Elements Extracted:**
- PMID, Title, Abstract
- Authors (names, affiliations, ORCID)
- Journal info, Publication date
- MeSH terms, Keywords
- Publication types
- DOI, PMC ID
- Language

### 3.2 Text Cleaning and Normalization

**Cleaning Steps:**
- Unicode normalization
- Citation removal [1,2,3]
- URL/email removal
- Whitespace normalization
- Special character handling
- Truncation for embedding model (256 tokens)

### 3.3 Entity Extraction

**Extracted Entities:**
- **SNPs:** rs numbers, HGVS protein/coding notation
- **Genes:** Gene symbols (validated against known genes)
- **Drugs:** Drug names (suffix patterns + known list)
- **Diseases:** Via NER (BioBERT/PubMedBERT)
- **Pathways:** Via knowledge base lookup

**Methods:**
- Regex patterns for SNPs
- Pattern matching for genes
- NER models for diseases
- Knowledge base linking for pathways

---

## 4. Embedding Generation

### 4.1 Model Selection

| Model | Dimensions | Speed (CPU) | Speed (GPU) | Quality | Size |
|-------|------------|-------------|-------------|---------|------|
| **all-MiniLM-L6-v2** | 384 | ~1,000/sec | ~15,000/sec | Excellent | 80 MB |
| all-mpnet-base-v2 | 768 | ~400/sec | ~8,000/sec | Better | 420 MB |
| PubMedBERT | 768 | ~300/sec | ~5,000/sec | Biomedical | 440 MB |

**Recommendation:** all-MiniLM-L6-v2
- RuVector optimized for 384d
- 4x faster than alternatives
- 50% smaller storage
- Quality sufficient (0.82 vs 0.84 on STS)

### 4.2 Chunking Strategy

**For Abstracts:**
- Title alone (always)
- Title + Abstract combined (if fits in 256 tokens)
- Abstract sections separately (if structured)

**For Full Text:**
- Title + Abstract (primary)
- Section-based chunks (512 tokens target)
- Overlap: 50 tokens for context preservation
- Tables converted to text

### 4.3 Processing Time Estimates

| Stage | 1M Papers | 5M Papers | 10M Papers |
|-------|-----------|-----------|------------|
| **Download (baseline)** | ~2 hours | ~5 hours | ~10 hours |
| **XML Parsing** | ~20 min | ~1.5 hours | ~3 hours |
| **Filtering** | ~10 min | ~45 min | ~1.5 hours |
| **Entity Extraction** | ~15 min | ~1 hour | ~2 hours |
| **Embedding (CPU)** | ~17 min | ~1.5 hours | ~3 hours |
| **Embedding (GPU)** | ~1 min | ~5 min | ~11 min |
| **RuVector Import** | ~10 min | ~45 min | ~1.5 hours |
| **Total (CPU)** | **~3 hours** | **~12 hours** | **~24 hours** |
| **Total (GPU)** | **~2.5 hours** | **~10 hours** | **~20 hours** |

---

## 5. Storage Pipeline

### 5.1 RuVector Ingestion

**Collections:**
- **articles:** Article metadata + embeddings
- **paper_chunks:** Individual chunks for RAG

**Properties:**
- PMID, PMCID, DOI (indexed)
- Title, Abstract
- Publication date, Journal
- MeSH terms, Keywords
- Extracted genes, SNPs, compounds
- Embeddings (384d, tiered compression)

### 5.2 Supabase Storage

**Schema:**
- research_articles table
- article_gene_links
- article_snp_links
- etl_update_log
- article_versions
- retractions

**Data Stored:**
- Full metadata
- Truncated abstract (2000 chars)
- Authors (JSONB)
- Relevance scores
- Update tracking

### 5.3 Storage Estimates

| Component | 1M Papers | 5M Papers | 10M Papers |
|-----------|-----------|-----------|------------|
| **Raw XML (temp)** | ~3 GB | ~15 GB | ~30 GB |
| **Processed JSON** | ~500 MB | ~2.5 GB | ~5 GB |
| **RuVector (tiered)** | ~400 MB | ~2 GB | ~4 GB |
| **Supabase metadata** | ~500 MB | ~2.5 GB | ~5 GB |
| **HNSW Index** | ~600 MB | ~3 GB | ~6 GB |
| **Total permanent** | **~2 GB** | **~10 GB** | **~20 GB** |

---

## 6. Update Pipeline

### 6.1 Daily/Weekly Update Strategy

**Daily (2 AM, after PubMed release):**
- Download daily update files
- Parse and filter new citations
- Extract entities
- Generate embeddings
- Update RuVector and Supabase
- **Volume:** ~1,000-2,000 papers/day
- **Time:** ~1-2 minutes

**Weekly (Sunday 3 AM):**
- Check for retractions
- Check for corrections
- Re-index if needed
- Clean up old files

**Monthly:**
- Full quality check
- Entity re-extraction for new tools
- Update citation graphs

### 6.2 Handling Retractions and Corrections

**Retractions:**
- Mark in database (is_retracted flag)
- Keep for historical record
- Flag in RuVector metadata
- Exclude from RAG results

**Corrections:**
- Add correction note
- Optionally re-embed if significant change
- Track version history

---

## 7. Complete Cost Summary

| Phase | Database Size | Processing | Storage | Monthly Total |
|-------|---------------|------------|---------|---------------|
| **MVP (1M papers)** | 2 GB | $0.30 | $0 (free tier) | **~$0.30** |
| **Standard (5M)** | 10 GB | $0.30 | $5 | **~$5.30** |
| **Comprehensive (10M)** | 20 GB | $0.30 | $10 | **~$10.30** |

**One-Time Setup:**
- CPU processing: $1-2
- GPU processing (optional): $1
- **Total:** ~$1-2

**Daily Updates:**
- Processing: ~$0.01/day
- Storage: Included
- **Monthly:** ~$0.30

---

## 8. Implementation Checklist

### Phase 1: PubMed Abstracts
- [ ] Set up PubMed baseline download pipeline
- [ ] Implement streaming XML parser
- [ ] Create genetics relevance filter (MeSH + keywords)
- [ ] Extract gene/SNP mentions with NER
- [ ] Generate abstract embeddings
- [ ] Build citation graph edges
- [ ] Load into RuVector articles collection

### Phase 2: Full-Text (PMC)
- [ ] Download PMC Open Access subset
- [ ] Implement JATS XML parser
- [ ] Create section-based chunking
- [ ] Extract tables and figure captions
- [ ] Generate chunk embeddings
- [ ] Build chunk navigation graph
- [ ] Link chunks to entities

### Phase 3: Updates & Maintenance
- [ ] Daily update cron job
- [ ] Weekly maintenance script
- [ ] Retraction monitoring
- [ ] Version tracking
- [ ] Quality metrics dashboard

---

## 9. Monitoring and Alerts

**Health Checks:**
- Last successful update timestamp
- Error rate in pipeline stages
- Processing time trends
- Storage growth rate
- Deduplication accuracy
- Embedding quality metrics

**Alerts:**
- No update in 48 hours
- Error rate > 5%
- Storage > 90% capacity
- Processing time > 2x baseline

**Dashboard Metrics:**
- Total papers processed
- Papers added per day
- Current database size
- Pipeline stage performance
- Error logs

---

## 10. Cron Schedule

```bash
# Daily updates (2 AM, after PubMed daily release)
0 2 * * * ubuntu python3 /scripts/etl/papers/pipeline.py --mode daily

# Weekly maintenance (Sunday 3 AM)
0 3 * * 0 ubuntu python3 /scripts/etl/papers/pipeline.py --mode weekly

# Monthly refresh check (1st of month, 4 AM)
0 4 1 * * ubuntu python3 /scripts/etl/papers/check_refresh.py
```

---

## Download

| Source | Method | URL/Command |
|--------|--------|-------------|
| **PubMed baseline** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/` |
| **PubMed updates** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/` |
| **PMC OA** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/` |

**Access Requirements:** Freely accessible; NCBI API key enables 10 requests/second (vs 3 without).

## Data Format

| Format | Description |
|--------|-------------|
| Input | XML (PubMed DTD), JATS XML |
| Output | JSON, Parquet |
| Embeddings | Float32 arrays (384-d) |
| Storage | RuVector collections |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `pipeline_id` | string | Pipeline run identifier | "daily_2026-01-22" |
| `stage` | string | Processing stage | "embedding_generation" |
| `records_processed` | integer | Count of records | 15000 |
| `status` | string | Pipeline status | "completed" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `processes` | Source File | 1:N |
| `produces` | Output Collection | 1:N |

## Sample Data

### Example Pipeline Run
```json
{
  "pipeline_id": "daily_2026-01-22",
  "stages": [
    {"name": "download", "records": 5000, "duration_s": 120},
    {"name": "parse", "records": 5000, "duration_s": 60},
    {"name": "embed", "records": 5000, "duration_s": 300}
  ],
  "total_duration_s": 480,
  "status": "completed"
}
```

### Sample Query Result
| run_date | records | duration | errors |
|----------|---------|----------|--------|
| 2026-01-22 | 5000 | 8m | 0 |
| 2026-01-21 | 4800 | 7m | 2 |

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| PubMed | Public domain | Yes |
| PMC OA | CC BY/CC0 | Yes |

## Data Set Size

| Metric | Value |
|--------|-------|
| Daily updates | ~5,000 records |
| Weekly maintenance | Retractions, corrections |
| Monthly refresh | Validation audit |
| Storage per million | ~3 GB with embeddings |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `ETL pipeline` | Extract, Transform, Load - data processing workflow | PubMed XML to RuVector |
| `baseline files` | Annual complete snapshot of PubMed citations | pubmed25n0001.xml.gz |
| `update files` | Daily incremental additions and corrections | Daily diff files |
| `MeSH filtering` | Selecting papers using Medical Subject Headings | Genetics [MeSH] |
| `entity extraction` | Identifying genes, SNPs, diseases in text | NER for rs numbers |
| `embedding generation` | Converting text to dense vectors for search | all-MiniLM-L6-v2 output |
| `retraction` | Withdrawal of a published paper due to errors or misconduct | Flagged in database |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `lxml` | Python library for fast XML parsing (10x faster than BeautifulSoup) | XML processing |
| `all-MiniLM-L6-v2` | Sentence transformer model producing 384-d embeddings | Embedding model |
| `NER` | Named Entity Recognition - ML-based entity extraction | BioBERT, PubMedBERT |
| `HGVS` | Human Genome Variation Society - notation for variants | p.Arg506Gln, c.1517G>A |
| `RuVector` | Vector database optimized for biomedical RAG applications | Storage backend |
| `Supabase` | PostgreSQL-based backend for relational metadata | Relational storage |
| `deduplication` | Removing duplicate papers across sources using DOI/PMID | Data quality |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ETL | Extract, Transform, Load | Data pipeline pattern |
| NER | Named Entity Recognition | Entity extraction |
| MeSH | Medical Subject Headings | PubMed vocabulary |
| HGVS | Human Genome Variation Society | Variant notation |
| FTP | File Transfer Protocol | Bulk download method |
| GWAS | Genome-Wide Association Study | Study type |
| SNP | Single Nucleotide Polymorphism | Genetic variant |
| MAF | Minor Allele Frequency | Population metric |
| PGx | Pharmacogenomics | Drug-gene field |
| ONNX | Open Neural Network Exchange | Model format |
| CPU | Central Processing Unit | Processing resource |
| GPU | Graphics Processing Unit | Fast embedding compute |

---

## References

- PubMed DTD: https://dtd.nlm.nih.gov/ncbi/pubmed/out/
- JATS XML: https://jats.nlm.nih.gov/
- MeSH Browser: https://meshb.nlm.nih.gov/
- RuVector: https://ruvector.io/
- Supabase: https://supabase.com/

---

*This document is part of the Gene Knowledge Base technical documentation.*

**Note:** The original pipeline-design.md file contains detailed Python implementation code for each stage. This summary focuses on architecture and key decisions. Refer to the original file in `/docs/data/source/literature/` for complete code examples.
