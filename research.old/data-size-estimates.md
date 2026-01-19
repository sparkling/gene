# Data Size Estimates

**Last Updated:** January 2026
**Vector Database:** RuVector (distributed, self-learning, with tiered compression)
**Purpose:** Comprehensive size estimates for source data downloads and target database

---

## Executive Summary

| Scope | Source Data | Target Database | RAM Required |
|-------|-------------|-----------------|--------------|
| **MVP** | ~20 GB | **~2-3 GB** | 2-4 GB |
| **Standard** | ~50 GB | **~8-12 GB** | 8-16 GB |
| **Comprehensive** | ~200 GB | **~20-30 GB** | 16-32 GB |
| **Full (with gnomAD)** | ~700 GB | **~80-120 GB** | 64-128 GB |

*RuVector's tiered compression reduces storage by 75-90% compared to traditional vector databases.*

---

## RuVector Compression Technology

### Tiered Compression System

RuVector automatically manages data across compression tiers based on access frequency:

| Tier | Access Freq | Format | Compression | 384d Vector | 768d Vector |
|------|-------------|--------|-------------|-------------|-------------|
| **Hot** | >80% | f32 | 1x | 1,536 bytes | 3,072 bytes |
| **Warm** | 40-80% | f16 | 2x | 768 bytes | 1,536 bytes |
| **Cool** | 10-40% | PQ8 | 8x | 192 bytes | 384 bytes |
| **Cold** | 1-10% | PQ4 | 16x | 96 bytes | 192 bytes |
| **Archive** | <1% | Binary | 32x | 48 bytes | 96 bytes |

### Effective Storage Multiplier

For typical access patterns (Pareto distribution):

| Tier | % of Data | Compression | Contribution |
|------|-----------|-------------|--------------|
| Hot | 10% | 1x | 0.100 |
| Warm | 20% | 2x | 0.100 |
| Cool | 30% | 8x | 0.038 |
| Cold | 25% | 16x | 0.016 |
| Archive | 15% | 32x | 0.005 |
| **Weighted Average** | 100% | - | **~0.26x** |

**Result:** ~74% average storage reduction compared to raw f32 vectors.

### Comparison with Other Vector Databases

| Database | 1M Vectors (384d) | Compression | Notes |
|----------|-------------------|-------------|-------|
| **RuVector** | ~200 MB | Automatic tiered | Self-learning |
| pgvector | ~1.5 GB | Manual halfvec only | PostgreSQL extension |
| Pinecone | ~2 GB | Limited | Managed service |
| Qdrant | ~1.5 GB | Scalar quantization | Self-hosted option |
| ChromaDB | ~3 GB | None | Simple setup |

---

## Part 1: Source Data Sizes

### Tier 1 - Core Foundation Sources

| Database | Records | Compressed | Uncompressed | Format | Notes |
|----------|---------|------------|--------------|--------|-------|
| **dbSNP** (full) | 1.2B SNPs | ~15 GB | ~150 GB | VCF | JSON significantly larger |
| **dbSNP** (practical subset) | 1-10M SNPs | ~500 MB | ~5 GB | JSON | Key variants only |
| **ClinVar** | 3M variants | ~500 MB | ~5 GB | XML/VCF | Weekly updates |
| **SNPedia** | 109K SNPs | ~100 MB | ~1 GB | Wiki export | Human-curated |
| **PharmGKB** | 11K interactions | ~50 MB | ~500 MB | TSV/JSON | Clean structured |
| **Reactome** | 2.5K pathways | ~200 MB | ~2 GB | BioPAX | CC0 license |
| **GWAS Catalog** | 500K associations | ~100 MB | ~1 GB | TSV | Trait-SNP links |
| **Tier 1 Total** | - | **~17 GB** | **~165 GB** | - | - |

### Tier 2 - Enrichment Sources

| Database | Records | Compressed | Uncompressed | Format | Notes |
|----------|---------|------------|--------------|--------|-------|
| **GeneCards** | ~20K genes | ~500 MB | ~5 GB | JSON | Gene summaries |
| **UniProt** (human) | ~20K proteins | ~500 MB | ~5 GB | XML/JSON | Human proteome |
| **WikiPathways** | 3K pathways | ~100 MB | ~1 GB | GPML | Community curated |
| **Open Targets** | Gene-disease | ~1 GB | ~10 GB | Parquet | Comprehensive |
| **ChEMBL** | 2M compounds | ~500 MB | ~3 GB | SQL/SDF | Drug data |
| **TCMSP** | 29K ingredients | ~200 MB | ~2 GB | Web scrape | TCM compounds |
| **NIH ODS** | Supplements | ~50 MB | ~500 MB | JSON | Official facts |
| **DSLD** | 140K products | ~500 MB | ~5 GB | JSON | Product labels |
| **Tier 2 Total** | - | **~4 GB** | **~32 GB** | - | - |

### Tier 3 - Expansion Sources

| Database | Records | Compressed | Uncompressed | Format | Notes |
|----------|---------|------------|--------------|--------|-------|
| **gnomAD** (summary) | Frequencies | ~50 GB | ~500 GB | VCF | Per-variant freq |
| **gnomAD** (full) | 807K individuals | ~500 GB | ~30 TB | VCF | AWS hosted |
| **PubMed** (baseline) | 36M citations | ~100 GB | ~500 GB | XML | 1,274 files |
| **PubMed** (abstracts only) | 36M abstracts | ~20 GB | ~100 GB | Text | For embeddings |
| **IMPPAT** | 18K phytochemicals | ~100 MB | ~1 GB | SDF/MOL | Ayurveda |
| **TCMBank** | 62K ingredients | ~300 MB | ~3 GB | Download | Largest TCM |
| **OMIM** | 16K genes | ~100 MB | ~1 GB | TXT | Gene-disease |
| **HPO** | 18K phenotypes | ~50 MB | ~500 MB | OBO/JSON | Phenotype ontology |
| **DisGeNET** | 1M+ associations | ~500 MB | ~5 GB | TSV | Gene-disease |
| **Tier 3 Total (no gnomAD)** | - | **~170 GB** | **~610 GB** | - | - |

### Total Source Data by Scope

| Scope | Tiers | Compressed | Uncompressed |
|-------|-------|------------|--------------|
| **MVP** | Tier 1 only | ~17 GB | ~165 GB |
| **Standard** | Tier 1 + 2 | ~21 GB | ~200 GB |
| **Comprehensive** | Tier 1 + 2 + 3 (partial) | ~100 GB | ~500 GB |
| **Full** | All tiers + gnomAD summary | ~220 GB | ~1.1 TB |

---

## Part 2: Target Database Sizes (RuVector)

### Vector Storage with RuVector Compression

**Using 384-dimensional embeddings (all-MiniLM-L6-v2):**

| Entity | Count | Raw f32 | RuVector (tiered) | With HNSW Index |
|--------|-------|---------|-------------------|-----------------|
| SNPs | 100K | 150 MB | **~40 MB** | ~60-80 MB |
| SNPs | 500K | 750 MB | **~200 MB** | ~300-400 MB |
| SNPs | 1M | 1.5 GB | **~400 MB** | ~600-800 MB |
| SNPs | 5M | 7.5 GB | **~2 GB** | ~3-4 GB |
| SNPs | 10M | 15 GB | **~4 GB** | ~6-8 GB |
| Genes | 20K | 30 MB | **~8 MB** | ~15 MB |
| Pathways | 5K | 8 MB | **~2 MB** | ~4 MB |
| Compounds | 100K | 150 MB | **~40 MB** | ~80 MB |
| Phenotypes | 50K | 75 MB | **~20 MB** | ~40 MB |
| Research | 1M | 1.5 GB | **~400 MB** | ~600-800 MB |
| Research | 5M | 7.5 GB | **~2 GB** | ~3-4 GB |
| Research | 10M | 15 GB | **~4 GB** | ~6-8 GB |

### Relational Data Storage

| Table | Records | Row Size (avg) | Total Size | Notes |
|-------|---------|----------------|------------|-------|
| **snps** | 1M | ~500 bytes | ~500 MB | Core SNP data |
| **snps** | 10M | ~500 bytes | ~5 GB | Comprehensive |
| **snp_interpretations** | 500K | ~2 KB | ~1 GB | Per-genotype info |
| **genes** | 20K | ~5 KB | ~100 MB | Gene summaries |
| **pathways** | 5K | ~10 KB | ~50 MB | Pathway definitions |
| **pathway_genes** (junction) | 500K | ~50 bytes | ~25 MB | Gene-pathway links |
| **compounds** | 100K | ~2 KB | ~200 MB | Supplements/drugs |
| **compound_gene_links** | 500K | ~100 bytes | ~50 MB | Drug targets |
| **phenotypes** | 50K | ~1 KB | ~50 MB | Conditions/symptoms |
| **snp_phenotype_links** | 5M | ~100 bytes | ~500 MB | Associations |
| **research_articles** | 1M | ~5 KB | ~5 GB | Citations |
| **research_articles** | 10M | ~5 KB | ~50 GB | Full PubMed |
| **entity_research_links** | 50M | ~50 bytes | ~2.5 GB | Citation links |
| **users** | 10K | ~1 KB | ~10 MB | User accounts |
| **user_genotypes** | 100M | ~50 bytes | ~5 GB | 10K users × 10K SNPs |

### Database Size by Configuration (RuVector)

#### MVP Configuration (~2-3 GB)

| Component | Size |
|-----------|------|
| SNPs (1M) + RuVector vectors | **800 MB** |
| Genes + vectors | 25 MB |
| Pathways + vectors | 10 MB |
| Compounds + vectors | 120 MB |
| Research (1M) + vectors | **800 MB** |
| Relational data & junctions | 500 MB |
| User data (1K users) | 100 MB |
| **Total** | **~2.3 GB** |
| With indexes & overhead | **~3 GB** |

#### Standard Configuration (~8-12 GB)

| Component | Size |
|-----------|------|
| SNPs (5M) + RuVector vectors | **3.5 GB** |
| Genes + vectors | 25 MB |
| Pathways + vectors | 10 MB |
| Compounds + vectors | 120 MB |
| Phenotypes + vectors | 60 MB |
| Research (5M) + vectors | **3.5 GB** |
| Relational data & junctions | 2 GB |
| User data (5K users) | 500 MB |
| **Total** | **~10 GB** |
| With indexes & overhead | **~12 GB** |

#### Comprehensive Configuration (~20-30 GB)

| Component | Size |
|-----------|------|
| SNPs (10M) + RuVector vectors | **7 GB** |
| Genes + vectors | 25 MB |
| Pathways + vectors | 10 MB |
| Compounds + vectors | 120 MB |
| Phenotypes + vectors | 60 MB |
| Research (10M) + vectors | **7 GB** |
| Relational data & junctions | 5 GB |
| User data (10K users) | 1 GB |
| **Total** | **~20 GB** |
| With indexes & overhead | **~25-30 GB** |

---

## Part 3: RAM Requirements (RuVector)

### Database RAM Guidelines

| Database Size | Minimum RAM | Recommended RAM | Notes |
|---------------|-------------|-----------------|-------|
| 3 GB (MVP) | 2 GB | 4 GB | Supabase Free/Pro |
| 12 GB (Standard) | 4 GB | 8 GB | Supabase Pro |
| 30 GB (Comprehensive) | 8 GB | 16 GB | Supabase Pro+ |
| 100 GB (Full) | 32 GB | 64 GB | Dedicated server |

### RuVector HNSW Index RAM

RuVector's HNSW implementation is memory-efficient:

| Vectors | Traditional HNSW | RuVector HNSW | Savings |
|---------|------------------|---------------|---------|
| 1M | ~4.5 GB | ~800 MB | 82% |
| 10M | ~45 GB | ~8 GB | 82% |

*RuVector applies compression to index structures, not just base vectors.*

### Working Memory

For query execution:
- **Hot cache:** 10% of vectors at f32 = minimal overhead
- **GNN layer:** ~50 MB for self-learning components
- **Query buffer:** ~100 MB for concurrent queries

---

## Part 4: Storage Costs (Revised)

### Supabase Pricing with RuVector

| Configuration | Database Size | Monthly Cost | Notes |
|---------------|---------------|--------------|-------|
| MVP | 3 GB | **$0** | Free tier (500 MB) or Pro |
| Standard | 12 GB | **$25** | Pro tier (8 GB included) |
| Comprehensive | 30 GB | **$25-50** | Pro + storage overage |
| Full | 100 GB | **$125-200** | Team tier |

### Comparison: pgvector vs RuVector Costs

| Configuration | pgvector Storage | pgvector Cost | RuVector Storage | RuVector Cost |
|---------------|------------------|---------------|------------------|---------------|
| MVP | 20 GB | $25/mo | 3 GB | **$0-25/mo** |
| Standard | 75 GB | $125/mo | 12 GB | **$25/mo** |
| Comprehensive | 200 GB | $300/mo | 30 GB | **$50/mo** |

**Annual Savings:** $900-3,000+ depending on configuration.

---

## Part 5: Data Transfer Estimates

### Initial Download Times

| Data Size | 100 Mbps | 1 Gbps | Notes |
|-----------|----------|--------|-------|
| 20 GB (MVP) | ~30 min | ~3 min | Tier 1 |
| 50 GB | ~1 hour | ~7 min | Tier 1+2 |
| 200 GB | ~4 hours | ~30 min | Comprehensive |

### RuVector Compression Time

| Operation | 1M Vectors | 10M Vectors |
|-----------|------------|-------------|
| Initial compression (PQ8) | ~2 min | ~20 min |
| Tier promotion/demotion | Real-time | Real-time |
| Full recompression | ~5 min | ~50 min |

---

## Part 6: Embedding Generation Time

### RuVector with Local ONNX Runtime

RuVector includes a pure JavaScript ONNX runtime - no Python or external APIs needed.

| Hardware | Model | Speed |
|----------|-------|-------|
| CPU (modern) | all-MiniLM-L6-v2 | ~1,000/sec |
| GPU (RTX 3080) | all-MiniLM-L6-v2 | ~15,000/sec |
| Browser (WASM) | all-MiniLM-L6-v2 | ~200/sec |

### Estimated Generation Times

| Records | CPU | GPU |
|---------|-----|-----|
| 100K | 2 min | 7 sec |
| 1M | 17 min | 1 min |
| 10M | 3 hours | 11 min |

---

## Part 7: Growth Projections

### User Data Growth

| Users | Genotypes/User | Total Genotypes | Storage |
|-------|----------------|-----------------|---------|
| 100 | 10K SNPs | 1M rows | 50 MB |
| 1K | 10K SNPs | 10M rows | 500 MB |
| 10K | 10K SNPs | 100M rows | 5 GB |

### Database Growth Rate (RuVector)

| Month | Users | Knowledge Base | User Data | Total |
|-------|-------|----------------|-----------|-------|
| 1 | 100 | 2 GB | 50 MB | **2 GB** |
| 3 | 500 | 3 GB | 250 MB | **3.3 GB** |
| 6 | 2K | 8 GB | 1 GB | **9 GB** |
| 12 | 10K | 20 GB | 5 GB | **25 GB** |

---

## Part 8: Recommendations

### MVP (Months 1-3)

```
Source Data Download: ~20 GB
Target Database: ~3 GB (RuVector)
RAM: 2-4 GB
Storage: 10 GB (headroom)
Cost: $0-25/month (Supabase Free/Pro)
```

**Focus:**
- 1M curated SNPs (SNPedia + ClinVar + key dbSNP)
- 20K genes
- 5K pathways (Reactome + WikiPathways)
- 100K compounds
- 1M research articles

### Standard (Months 4-8)

```
Source Data Download: ~50 GB
Target Database: ~12 GB (RuVector)
RAM: 8-16 GB
Storage: 30 GB (headroom)
Cost: $25-50/month (Supabase Pro)
```

**Add:**
- 5M SNPs
- gnomAD frequencies (summary)
- Full pathway coverage
- Traditional medicine compounds
- 5M research articles

### Comprehensive (Months 9-12)

```
Source Data Download: ~200 GB
Target Database: ~30 GB (RuVector)
RAM: 16-32 GB
Storage: 80 GB (headroom)
Cost: $50-100/month (Supabase Pro/Team)
```

**Add:**
- 10M+ SNPs
- Full gnomAD summary
- Complete PubMed coverage
- All specialized databases

---

## Quick Reference Card

| Metric | MVP | Standard | Comprehensive |
|--------|-----|----------|---------------|
| **SNPs** | 1M | 5M | 10M |
| **Research Articles** | 1M | 5M | 10M |
| **Source Download** | 20 GB | 50 GB | 200 GB |
| **Database Size (RuVector)** | **3 GB** | **12 GB** | **30 GB** |
| **Database Size (pgvector)** | 20 GB | 75 GB | 200 GB |
| **RAM Needed** | 2-4 GB | 8-16 GB | 16-32 GB |
| **Monthly Cost** | **$0-25** | **$25-50** | **$50-100** |
| **Compression Savings** | 85% | 84% | 85% |

---

## RuVector CLI Commands

```bash
# Install
npm install ruvector

# Quick start
npx ruvector

# Compress embeddings
npx ruvector gnn compress -f embeddings.json -l pq8 -o compressed.json

# Available compression levels
# none (f32), half (f16), pq8, pq4, binary

# Initialize for Claude Code
npx @ruvector/cli hooks init
```
