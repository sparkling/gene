---
id: research-literature-priority
title: "Research Papers Strategy: Final Recommendations"
type: research-priority
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [research, priorities, swarm-synthesis, literature, pubmed, pmc, abstracts]
---

# Research Papers Strategy: Final Recommendations

**Document Type:** Synthesis & Recommendations
**Date:** January 2026
**Swarm Analysis:** 5 parallel research agents
**Total Research Output:** 183KB across 5 detailed documents
**Parent:** [_index.md](./_index.md)

---

## Executive Summary

After comprehensive analysis by a claude-flow research swarm, here are the definitive recommendations for research paper integration in the Gene platform's knowledge base.

### Key Decision: ABSTRACTS-ONLY with PubMed/PMC

| Aspect | Recommendation | Rationale |
|--------|----------------|-----------|
| **Primary Source** | PubMed/PMC | Public domain, 39M+ papers, legal clarity |
| **Content Scope** | Abstracts only (initially) | 74% storage reduction, simpler processing |
| **Alternative Sources** | DO NOT USE | Extreme legal risk for commercial platforms |
| **Genetics Subset** | 3-5M relevant papers | Filter by MeSH terms |
| **Storage (RuVector)** | ~400MB - 1.2GB | With tiered compression |
| **Annual Cost** | ~$0 | Self-hosted, open access data |

---

## 1. Data Source Recommendation

### USE: Public Resources (Multi-Source Strategy)

**Note:** PubMed alone misses ~40% of relevant papers. Use multiple sources for comprehensive coverage.

| Source | Records | Bulk Download | License | Risk | Coverage Gap |
|--------|---------|---------------|---------|------|--------------|
| **PubMed** | 39M+ | FTP (28GB) | Public Domain | NONE | Primary |
| **OpenAlex** | 250M+ | S3 (300GB) | CC0 | NONE | Fills gaps |
| **PMC Open Access** | 3.4M | FTP/S3 (50GB) | CC licenses | LOW | Full-text |
| **Europe PMC** | 43M+ | FTP | Open | LOW | Annotations |

See `papers-coverage-analysis.md` for detailed multi-source strategy and download scripts.

### DO NOT USE: Shadow Libraries

| Source | Status | Legal Risk | Recommendation |
|--------|--------|------------|----------------|
| **Sci-Hub** | Active | EXTREME | Avoid completely |
| **Anna's Archive** | Active | EXTREME | Avoid completely |
| **Library Genesis** | Active | EXTREME | Avoid completely |
| **Z-Library** | Active | EXTREME | Avoid completely |

**Legal Precedent Warning:** Multiple AI companies (Anthropic, Meta, OpenAI) face $1B+ lawsuits specifically for downloading from shadow libraries. A bootstrapped startup cannot absorb this risk.

---

## 2. Abstracts vs Full-Text: Definitive Recommendation

### Recommendation: ABSTRACTS-ONLY (Phase 1)

**Rationale:**

| Factor | Abstracts | Full-Text | Winner |
|--------|-----------|-----------|--------|
| **Availability** | 36M (100%) | 5-8M (15-20%) | Abstracts |
| **Legal Clarity** | Public domain | Mixed licenses | Abstracts |
| **Storage per paper** | 1.5 KB | 30-100 KB | Abstracts |
| **Processing complexity** | Simple XML | PDF extraction | Abstracts |
| **RAG effectiveness** | Good (key findings) | Better (details) | Full-text |
| **Time to implement** | Days | Weeks | Abstracts |

**Storage Comparison (1M papers with embeddings):**
- Abstracts + 1 embedding: ~3 GB
- Abstracts + 5 chunk embeddings: ~9 GB
- Full-text + 20 chunk embeddings: ~60 GB

### Phase 2: Hybrid Approach (Future)

After MVP launch, add PMC Open Access full-text for:
- High-impact papers (citation count > 100)
- Papers directly linked to key SNPs
- Review articles and meta-analyses

---

## 3. Recommended Architecture

### Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                    PubMed FTP (Baseline + Updates)              │
│                    ftp://ftp.ncbi.nlm.nih.gov/pubmed/           │
└─────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                         FILTER PIPELINE                          │
│  MeSH Terms: Genetics, Genomics, Pharmacogenomics, SNP          │
│  Result: 39M → 3-5M relevant papers                             │
└─────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                       PROCESSING PIPELINE                        │
│  1. Parse XML (lxml)                                            │
│  2. Extract: PMID, title, abstract, MeSH, authors, date         │
│  3. Entity extraction: gene names, SNP rs numbers, drugs        │
│  4. Generate embeddings (all-MiniLM-L6-v2, 384d)                │
└─────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                          RUVECTOR STORAGE                        │
│  Collection: articles                                           │
│  Vectors: 384d with tiered compression                          │
│  Graph: CITED_IN relationships to SNPs, genes, compounds        │
│  Size: ~400MB - 1.2GB (1-3M papers)                             │
└─────────────────────────────────────────────────────────────────┘
```

### MeSH Terms for Filtering

**High Priority (must include):**
- Genetics [MeSH]
- Genomics [MeSH]
- Polymorphism, Single Nucleotide [MeSH]
- Pharmacogenetics [MeSH]
- Gene Expression [MeSH]

**Medium Priority (include if space allows):**
- Genetic Predisposition to Disease [MeSH]
- Genotype [MeSH]
- Alleles [MeSH]
- Mutation [MeSH]

**Domain-Specific:**
- Autoimmune Diseases/genetics [MeSH]
- Methylation [MeSH]
- Vitamins/genetics [MeSH]

---

## 4. Implementation Timeline

### Week 1-2: Infrastructure
- [ ] Set up PubMed FTP download scripts
- [ ] Implement XML parser for PubMed format
- [ ] Create MeSH filtering pipeline
- [ ] Set up RuVector collection schema

### Week 3-4: Initial Load
- [ ] Download PubMed baseline (1,219 files, ~100 GB compressed)
- [ ] Filter to genetics subset (~3-5M papers)
- [ ] Generate embeddings (batch processing)
- [ ] Load into RuVector with graph relationships

### Week 5+: Integration
- [ ] Build RAG retrieval API
- [ ] Create SNP-to-paper linking
- [ ] Implement daily update pipeline
- [ ] Add citation extraction

---

## 5. Cost Analysis

### One-Time Setup

| Item | Cost |
|------|------|
| PubMed baseline download | $0 (FTP) |
| Processing compute (initial) | $10-50 (spot instances) |
| Embedding generation | $0 (local model) |
| **Total Setup** | **~$50** |

### Ongoing Monthly

| Item | Cost |
|------|------|
| Storage (RuVector self-hosted) | $0 (included in VPS) |
| Daily updates processing | $0 (cron job) |
| API calls | $0 (bulk download) |
| **Total Monthly** | **~$0** |

### Comparison: If Using Full-Text + Commercial Sources

| Item | Cost |
|------|------|
| DrugBank license | $10,000+/year |
| NatMed Pro | $5,000+/year |
| Storage (10x more) | $50+/month |
| Legal risk (shadow libraries) | Potentially $1B+ |
| **Total** | **Unaffordable** |

---

## 6. RuVector Schema for Papers

```typescript
// RuVector collection configuration
const articlesCollection = {
  name: 'articles',
  dimension: 384,
  distanceMetric: 'cosine',
  compression: 'auto',  // Tiered: f32 → f16 → PQ8 → PQ4
  properties: {
    pmid: { type: 'number', indexed: true },
    title: { type: 'string' },
    abstract: { type: 'string' },
    publication_date: { type: 'date', indexed: true },
    journal: { type: 'string', indexed: true },
    mesh_terms: { type: 'string[]', indexed: true },
    authors: { type: 'string[]' },
    // Extracted entities
    mentioned_genes: { type: 'string[]', indexed: true },
    mentioned_snps: { type: 'string[]', indexed: true },
    mentioned_compounds: { type: 'string[]', indexed: true },
  }
};

// Graph relationships (Cypher)
// (:SNP)-[:CITED_IN]->(:Article)
// (:Gene)-[:CITED_IN]->(:Article)
// (:Compound)-[:CITED_IN]->(:Article)
// (:Article)-[:CITES]->(:Article)
```

---

## 7. Key Queries for RAG

### Find papers about a SNP

```cypher
MATCH (snp:SNP {rs_number: 'rs1801133'})-[:CITED_IN]->(article:Article)
RETURN article.pmid, article.title, article.abstract
ORDER BY article.publication_date DESC
LIMIT 10
```

### Semantic search for user question

```typescript
// User asks: "What supplements help with MTHFR mutations?"
const results = await ruvector.search('articles', {
  query: "MTHFR mutation supplements folate methylfolate treatment",
  filter: { mesh_terms: { $contains: 'Folic Acid' } },
  topK: 10
});
```

### Find papers connecting genes to conditions

```cypher
MATCH (gene:Gene {symbol: 'MTHFR'})-[:CITED_IN]->(article:Article)
      <-[:CITED_IN]-(phenotype:Phenotype {type: 'disease'})
WHERE article.publication_date > date('2020-01-01')
RETURN DISTINCT phenotype.name, COUNT(article) as paper_count
ORDER BY paper_count DESC
```

---

## 8. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| PubMed API rate limiting | Medium | Low | Use bulk FTP, respect limits |
| Storage growth | Low | Medium | Tiered compression, archival |
| Data quality issues | Medium | Medium | Validation pipeline |
| Missing recent papers | Low | Low | Daily update sync |
| Legal (public sources) | Very Low | Low | All public domain |
| Legal (shadow libraries) | N/A | EXTREME | Not using |

---

## 9. Summary of Research Documents

| Document | Size | Key Content |
|----------|------|-------------|
| `papers-public-sources.md` | 23 KB | PubMed, PMC, Europe PMC APIs, FTP, formats |
| `papers-alternative-sources.md` | 15 KB | Sci-Hub, Anna's Archive - legal warnings |
| `papers-data-structures.md` | 49 KB | XML schemas, parsing code, embedding strategy |
| `papers-pipeline-design.md` | 78 KB | Full ETL pipeline with Python code |
| `papers-abstracts-vs-fulltext.md` | 18 KB | Detailed comparison, storage analysis |
| **Total Research** | **183 KB** | Comprehensive analysis |

---

## 10. Final Recommendation

### For a bootstrapped genetics platform:

1. **Use PubMed abstracts only** - Legal, free, comprehensive
2. **Filter to 3-5M genetics papers** - MeSH term filtering
3. **Store in RuVector** - ~1GB with tiered compression
4. **Build daily update pipeline** - Stay current
5. **DO NOT use shadow libraries** - Legal risk is existential

### Future expansion path:

1. **Phase 1 (MVP):** PubMed abstracts (3-5M papers)
2. **Phase 2:** Add PMC Open Access full-text (high-impact papers)
3. **Phase 3:** Add Europe PMC for broader coverage
4. **Never:** Shadow libraries for commercial use

---

## Download

| Resource | Method | URL |
|----------|--------|-----|
| PubMed | FTP | ftp://ftp.ncbi.nlm.nih.gov/pubmed/ |
| PubMed Central | FTP/API | https://www.ncbi.nlm.nih.gov/pmc/ |
| OpenAlex | S3 | https://www.openalex.org/download |
| Europe PMC | FTP | https://europepmc.org/downloads |

**Access Requirements:** Public domain data, no registration required for bulk downloads.

## Data Format

| Format | Description |
|--------|-------------|
| XML | PubMed article records |
| TSV | Tab-delimited metadata |
| JSON | Structured data from APIs |
| UTF-8 | Text encoding standard |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `pmid` | number | PubMed identifier | 12345678 |
| `title` | string | Article title | "Genetics of diabetes" |
| `abstract` | string | Article summary | "Abstract text..." |
| `mesh_terms` | array | Medical Subject Headings | ["Genetics", "Genomics"] |
| `publication_date` | date | Publication date | "2023-06-15" |

## Sample Data

### Example Article Record
```json
{
  "pmid": 36418627,
  "title": "Pharmacogenomic insights into drug metabolism and response",
  "abstract": "This study examines genetic variants affecting drug response...",
  "journal": "Nature Genetics",
  "publication_date": "2023-01-15",
  "mesh_terms": ["Pharmacogenetics", "Genetic Variation", "Drug Response"],
  "authors": ["Smith J", "Johnson M"],
  "mentioned_genes": ["CYP2D6", "MTHFR"],
  "mentioned_snps": ["rs1801133", "rs5030655"],
  "embedding_vector": [0.125, -0.034, ...]
}
```

## License

| Resource | License | Notes |
|----------|---------|-------|
| PubMed | Public Domain | Free for all uses |
| PMC Open Access | CC licenses | Varies by article |
| OpenAlex | CC0 | Public domain |
| Europe PMC | Open | CC license terms |

## Data Set Size

| Metric | Value |
|--------|-------|
| Total PubMed records | ~39M |
| Genetics subset | ~3-5M |
| With embeddings | ~400MB-1.2GB |
| Weekly new articles | ~20K |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Abstract | Summary of a research paper's key findings | 200-300 word synopsis |
| Full-Text | Complete research paper content | PDF article |
| RAG | Retrieval-Augmented Generation using retrieved documents for AI responses | Semantic search + LLM |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| PubMed | NCBI's biomedical literature database | 39M+ citations |
| PMC | PubMed Central full-text archive | Open access papers |
| Europe PMC | European version of PMC | Extended coverage |
| OpenAlex | Open scholarly metadata database | 250M+ works |
| MeSH Terms | Medical Subject Headings vocabulary | Literature indexing |
| Embedding | Vector representation of text | Semantic search |
| RuVector | Project's vector database | Embedding storage |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | Data access |
| CC | Creative Commons | License type |
| CC0 | Creative Commons Zero (Public Domain) | No restrictions |
| FTP | File Transfer Protocol | Bulk downloads |
| GB | Gigabyte | Storage unit |
| MeSH | Medical Subject Headings | Controlled vocabulary |
| MVP | Minimum Viable Product | Development phase |
| NCBI | National Center for Biotechnology Information | US agency |
| PDF | Portable Document Format | Document format |
| PMC | PubMed Central | Full-text archive |
| RAG | Retrieval-Augmented Generation | AI technique |
| XML | Extensible Markup Language | Data format |

---

*This recommendation synthesizes findings from 5 parallel research agents analyzing data sources, formats, pipelines, and legal considerations for research paper integration.*
