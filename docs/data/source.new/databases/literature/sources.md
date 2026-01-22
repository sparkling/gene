---
id: literature-sources
title: Literature & Research Paper Sources
category: literature
tier: 1
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [literature, papers, pubmed, pmc, openalex, research]
---

# Literature & Research Paper Sources

**Document ID:** 43-61-LITERATURE-SOURCES
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

For a bootstrapped genetics platform, use **PubMed abstracts as primary source** (39M+ papers, public domain, free bulk download) supplemented by **OpenAlex** (250M works, CC0) for coverage gaps and **PMC Open Access** (3.4M articles) for full-text depth. Do NOT use shadow libraries (Sci-Hub, LibGen, Anna's Archive) due to extreme legal risk. Abstracts-only approach provides 80% RAG effectiveness at 5% of storage cost, with a clear path to hybrid full-text enrichment.

---

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **Primary Source** | PubMed/MEDLINE | Public domain, best genetics coverage (95%), MeSH filtering |
| **Coverage Gap Fill** | OpenAlex | 250M works, CC0 license, preprints + international journals |
| **Full-Text Source** | PMC Open Access | 3.4M articles, free bulk download, known licenses |
| **Content Scope** | Abstracts-first | 74% storage reduction, simpler processing, legal clarity |
| **Shadow Libraries** | DO NOT USE | Extreme legal risk ($1B+ lawsuits against AI companies) |
| **Embedding Model** | all-MiniLM-L6-v2 | 384d, RuVector optimized, local ONNX runtime |
| **Target Subset** | 3-5M genetics papers | MeSH term filtering from 39M total |
| **Storage Estimate** | 400MB - 1.2GB | RuVector with tiered compression |

---

## Database Catalog

### 1. PubMed / MEDLINE

| Attribute | Value |
|-----------|-------|
| **Provider** | National Center for Biotechnology Information (NCBI), NLM, NIH |
| **URL** | https://pubmed.ncbi.nlm.nih.gov/ |
| **Total Records** | 39+ million citations |
| **Annual Growth** | ~1 million new citations/year |
| **Date Range** | 1809 to present (comprehensive from 1960s) |
| **Content Type** | Abstracts and metadata only (no full text) |
| **Genetics Relevance** | ~3-5M papers (~10-12% of total) |
| **License** | Public domain (US Government) |
| **Commercial Use** | YES - abstracts are public domain |

**API Access:**

| Endpoint | Function | Rate Limit |
|----------|----------|------------|
| `esearch.fcgi` | Text search, returns PMIDs | 3/sec (anon), 10/sec (API key) |
| `efetch.fcgi` | Retrieve full records | Same |
| `esummary.fcgi` | Document summaries | Same |
| `elink.fcgi` | Related records | Same |

**Bulk Download:**

| Resource | URL | Size | Update |
|----------|-----|------|--------|
| Baseline | ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/ | ~28 GB compressed | Annual (December) |
| Updates | ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/ | ~5 GB rolling | Daily |

**Data Format:** XML (PubMed DTD), JSON for esearch/esummary only

**Key Fields:**
- PMID (unique identifier)
- Title, Abstract
- Authors (names, affiliations, ORCID)
- Journal info (name, ISSN, volume, issue, pages)
- Publication date
- MeSH terms (Medical Subject Headings)
- Keywords, Publication type
- DOI, PMC ID links

**Genetics Subset Estimate:**
- MeSH "Genetics" and "Genomics": ~3-5M articles
- SNP-related articles: ~200,000+
- Pharmacogenomics: ~50,000+
- GWAS studies: ~250,000+

---

### 2. PubMed Central (PMC)

| Attribute | Value |
|-----------|-------|
| **Provider** | NCBI at NLM |
| **URL** | https://pmc.ncbi.nlm.nih.gov/ |
| **Total Records** | 10+ million full-text articles |
| **Open Access Subset** | ~3.4 million reusable articles |
| **Date Range** | Late 1700s to present |
| **Content Type** | Full-text articles with figures, tables |
| **Genetics Relevance** | ~1-2M full-text genetics articles |
| **License** | Per-article CC licenses |
| **Commercial Use** | ~35% (CC0, CC BY, CC BY-SA, CC BY-ND) |

**Open Access License Breakdown:**

| License Category | Articles | Commercial Use |
|------------------|----------|----------------|
| CC0, CC BY, CC BY-SA, CC BY-ND | ~1.5M | YES |
| CC BY-NC, CC BY-NC-SA, CC BY-NC-ND | ~1.5M | NO |
| Other/Custom | ~400K | CHECK EACH |

**Bulk Download:**

| Resource | URL | Size |
|----------|-----|------|
| FTP | ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/ | ~50 GB compressed |
| AWS S3 | s3://pmc-oa-opendata | Same |
| File List | ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv | License tracking |

**Directory Structure:**
```
/pub/pmc/
├── oa_bulk/
│   ├── oa_comm/        # Commercial use allowed
│   ├── oa_noncomm/     # Non-commercial only
│   └── oa_other/       # Custom/no license
├── manuscript/
└── historical_ocr/
```

**Data Format:** JATS XML (Journal Article Tag Suite, NISO Z39.96)

**Key Elements:**
- Complete article body with sections
- Tables, figures with captions
- Supplementary materials
- Full reference lists with linked PMIDs
- Author affiliations, funding info

---

### 3. Europe PMC

| Attribute | Value |
|-----------|-------|
| **Provider** | European Bioinformatics Institute (EMBL-EBI) |
| **URL** | https://europepmc.org/ |
| **Total Abstracts** | 43+ million (includes all PubMed) |
| **Full-Text Articles** | 9+ million |
| **Open Access** | 6.5 million |
| **Preprints** | 880,000+ (from 35 servers) |
| **PhD Theses** | 176,000+ (UK EThOS) |
| **Genetics Relevance** | ~4-5M genetics abstracts |
| **License** | Open API, per-article content licenses |
| **Commercial Use** | Check per-article |

**Unique Value:**
- Text-mined annotations (genes, proteins, chemicals, diseases)
- Database cross-references (UniProt, ENA, PDB, ChEMBL)
- Preprint-to-publication tracking
- Citation network (19.4M+ reference lists)

**API Access:**

| Feature | Details |
|---------|---------|
| Base URL | https://www.ebi.ac.uk/europepmc/webservices/rest/ |
| Authentication | None required |
| Rate Limit | ~10 requests/second |
| Format | JSON, XML, Dublin Core |

**Key Endpoints:**
- `/search` - Query publications
- `/{source}/{id}/citations` - Citation data
- `/{source}/{id}/references` - Reference lists
- `/{source}/{id}/textMinedTerms` - Extracted entities
- `/{id}/fullTextXML` - Full-text XML

**Bulk Download:**

| Dataset | Format | Update |
|---------|--------|--------|
| Open Access Subset | XML, PDF | Weekly |
| Metadata | XML | Weekly |
| PMID-PMCID-DOI Mappings | CSV | Monthly |
| Accession Numbers | CSV | Weekly |

---

### 4. OpenAlex

| Attribute | Value |
|-----------|-------|
| **Provider** | OurResearch (nonprofit) |
| **URL** | https://openalex.org/ |
| **Total Works** | 250+ million |
| **Daily Growth** | ~50,000 works/day |
| **Content Type** | Metadata + abstracts where available |
| **Coverage** | PubMed + Crossref + arXiv + repositories |
| **License** | CC0 (public domain) |
| **Commercial Use** | YES |

**Fills PubMed Gaps:**

| Content Type | PubMed | OpenAlex |
|--------------|--------|----------|
| Preprints | 5% | 40% |
| Conference papers | 20% | 60% |
| International journals | 70% | 85% |
| TCM/Integrative | 30% | 50% |

**Bulk Download:**

| Resource | URL | Size |
|----------|-----|------|
| Full Snapshot | s3://openalex | ~300 GB |
| Works only | s3://openalex/data/works | ~250 GB |

```bash
# Download with AWS CLI (no authentication)
aws s3 sync s3://openalex/data/works --no-sign-request /data/openalex/
```

**API Access:**

| Feature | Details |
|---------|---------|
| Base URL | https://api.openalex.org/ |
| Rate Limit | 100,000 requests/day |
| Authentication | None (polite pool with email) |
| Format | JSON |

**Data Structure:**
- Works (papers)
- Authors (90M)
- Institutions (100K)
- Concepts (65K topics)
- Venues/Sources (125K journals)

---

### 5. Semantic Scholar

| Attribute | Value |
|-----------|-------|
| **Provider** | Allen Institute for AI |
| **URL** | https://www.semanticscholar.org/ |
| **Total Papers** | 200+ million |
| **Content Type** | Metadata, abstracts, citations |
| **Unique Value** | AI-extracted entities, influence scores |
| **License** | Academic use (commercial requires approval) |
| **Commercial Use** | BY REQUEST ONLY |

**API Access:**

| Feature | Details |
|---------|---------|
| Rate Limit | 100 requests/second |
| Bulk Access | Request required |
| Features | Citation context, TLDR summaries |

**Recommendation:** Use for supplementary citation analysis, not primary source.

---

### 6. Preprint Servers

#### bioRxiv / medRxiv

| Attribute | bioRxiv | medRxiv |
|-----------|---------|---------|
| **Provider** | Cold Spring Harbor Laboratory |
| **Focus** | Biology | Medicine/Health |
| **Total** | 200K+ | 50K+ |
| **API** | Yes | Yes |
| **Bulk** | Yes (with agreement) |
| **License** | CC BY (most) |
| **Commercial** | YES |

**API Base:** https://api.biorxiv.org/

**Value:** Access to research 6-12 months before peer-reviewed publication.

---

### 7. Crossref

| Attribute | Value |
|-----------|-------|
| **Provider** | Crossref (nonprofit) |
| **URL** | https://www.crossref.org/ |
| **Total Records** | 156+ million metadata records |
| **Content Type** | DOI metadata only (no abstracts) |
| **License** | CC0 for most metadata |
| **Commercial Use** | YES |

**API Access:**

| Feature | Details |
|---------|---------|
| Rate Limit | 50 requests/second (polite pool) |
| Bulk | Annual data file available |
| Use Case | DOI resolution, citation metadata |

---

### 8. Unpaywall

| Attribute | Value |
|-----------|-------|
| **Provider** | OurResearch |
| **URL** | https://unpaywall.org/ |
| **Coverage** | 20+ million free legal PDF locations |
| **Success Rate** | ~52% of papers have free legal version |
| **License** | Open |
| **Commercial Use** | YES |

**Use Case:** Find legal full-text locations for papers discovered via metadata sources.

---

## Coverage Analysis

### What PubMed Covers Well

| Topic | Coverage | Notes |
|-------|----------|-------|
| SNP/GWAS studies | ~95% | Core strength |
| Pharmacogenomics | ~90% | Strong NIH funding |
| Clinical genetics | ~85% | Most major journals indexed |
| Human disease genetics | ~90% | Comprehensive |

### What PubMed Misses

| Topic | Coverage | Gap Source |
|-------|----------|------------|
| Supplements/nutrition | ~60% | Non-indexed journals |
| Traditional medicine (TCM) | ~30% | Regional databases |
| Ayurveda research | ~25% | Indian journals |
| Preprints | ~5% | Not indexed until published |
| International journals | ~70% | Non-English journals |
| Conference papers | ~20% | Rarely indexed |

### Multi-Source Coverage Strategy

```
LAYER 1: PubMed (Primary - Genetics Focus)
├── 39M biomedical papers
├── Best MeSH term filtering
└── Download: FTP baseline (28GB compressed)

LAYER 2: OpenAlex (Breadth - Fill Coverage Gaps)
├── 250M works (6x PubMed)
├── Includes preprints, conference papers, books
├── Better international coverage
└── Download: S3 snapshot (300GB)

LAYER 3: PMC Open Access (Depth - Full Text)
├── 3.4M full-text articles
├── Complete methods, results, figures
└── Download: FTP/S3 (50GB)

LAYER 4: Specialized (Domain-Specific)
├── Europe PMC: 43M with text-mined annotations
├── Semantic Scholar: AI-extracted entities
└── Preprint servers: bioRxiv, medRxiv
```

---

## Full-Text vs Abstracts Analysis

### Recommendation: ABSTRACTS-FIRST

| Factor | Abstracts | Full-Text | Winner |
|--------|-----------|-----------|--------|
| **Availability** | 36M (100%) | 5-8M (15-20%) | Abstracts |
| **Legal Clarity** | Public domain | Mixed licenses | Abstracts |
| **Storage per paper** | 1.5 KB | 30-100 KB | Abstracts |
| **Processing complexity** | Simple XML | PDF extraction | Abstracts |
| **RAG effectiveness** | Good (key findings) | Better (details) | Full-text |
| **Time to implement** | Days | Weeks | Abstracts |

### Storage Comparison (1M Papers)

| Configuration | Storage |
|---------------|---------|
| Abstract only (text) | 1.5 GB |
| Abstract + 1 embedding | 3 GB |
| Abstract + 5 chunk embeddings | 9 GB |
| Full text + 20 chunk embeddings | 60 GB |
| Full text + embeddings + metadata | 100 GB |

### What Abstracts Can Answer

| Query Type | Quality | Notes |
|------------|---------|-------|
| "What diseases are associated with rs12345?" | Excellent | Key findings in abstract |
| "Which SNPs affect drug X metabolism?" | Good | Major associations mentioned |
| "Is gene Y implicated in condition Z?" | Excellent | Core hypothesis stated |
| Drug-gene interactions | Good | Existence yes, details no |

### What Requires Full-Text

| Query Type | Abstract Quality | Full-Text Quality |
|------------|------------------|-------------------|
| "What was the exact p-value for rs12345?" | Poor | Excellent |
| "What covariates were adjusted for?" | Very Poor | Excellent |
| "What population was studied?" | Moderate | Excellent |
| "What were the effect sizes?" | Poor | Excellent |
| "What methods were used for genotyping?" | Very Poor | Excellent |

---

## Shadow Libraries (DO NOT USE)

### Legal Risk Assessment

| Source | Status | Legal Risk | Recommendation |
|--------|--------|------------|----------------|
| **Sci-Hub** | Active | EXTREME | AVOID |
| **Anna's Archive** | Active | EXTREME | AVOID |
| **Library Genesis** | Active | EXTREME | AVOID |
| **Z-Library** | Active | EXTREME | AVOID |

### Legal Precedents (2024-2026)

- **Bartz v. Anthropic (2025):** Judge ruled downloading from LibGen = willful infringement. Settlement: $1.5 billion.
- **Meta Litigation (2025):** Unsealed emails revealed 81+ TB downloaded from Anna's Archive torrents.
- **OpenAI (2024-2025):** Court allowed plaintiffs to pursue shadow library downloading as separate claim.

### Bottom Line

For a bootstrapped startup, using shadow library content for commercial purposes carries **existential legal risk**. All recommended sources above provide legal alternatives with sufficient coverage.

---

## Integration Recommendations

### Phase 1: MVP (Weeks 1-4)

| Component | Scope | Cost |
|-----------|-------|------|
| PubMed genetics abstracts | 3-5M papers | Free |
| Embeddings (384d) | 3-5M vectors | ~$50 compute |
| Vector storage (RuVector) | 4.5 GB | ~$5/month |
| **Total** | - | **~$50 setup, ~$5/month** |

**Implementation:**
1. Download PubMed baseline via FTP
2. Filter by genetics MeSH terms
3. Parse XML, extract abstracts
4. Generate embeddings (all-MiniLM-L6-v2)
5. Load into RuVector with graph relationships

### Phase 2: Coverage Expansion (Months 2-3)

| Component | Scope | Cost |
|-----------|-------|------|
| Phase 1 content | 3-5M abstracts | Included |
| OpenAlex genetics filter | 5-10M additional | Free |
| Deduplication by DOI/PMID | - | - |
| **Total** | 8-15M papers | **~$100 compute** |

### Phase 3: Full-Text Enrichment (Months 4-6)

| Component | Scope | Cost |
|-----------|-------|------|
| PMC OA genetics (CC BY) | 300K full-texts | Free |
| Additional chunk embeddings | 5M chunks | ~$200 compute |
| **Total** | 15 GB storage | **~$15/month** |

### MeSH Terms for Filtering

**High Priority:**
- Genetics [MeSH]
- Genomics [MeSH]
- Polymorphism, Single Nucleotide [MeSH]
- Pharmacogenetics [MeSH]
- Gene Expression [MeSH]
- Genome-Wide Association Study [MeSH]

**Medium Priority:**
- Genetic Predisposition to Disease [MeSH]
- Genotype [MeSH]
- Alleles [MeSH]
- Mutation [MeSH]
- Epigenesis, Genetic [MeSH]

---

## Technical Specifications

### RuVector Schema

```typescript
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
    mentioned_genes: { type: 'string[]', indexed: true },
    mentioned_snps: { type: 'string[]', indexed: true },
    mentioned_compounds: { type: 'string[]', indexed: true },
  }
};

// Graph Relationships (Cypher)
// (:SNP)-[:CITED_IN]->(:Article)
// (:Gene)-[:CITED_IN]->(:Article)
// (:Article)-[:CITES]->(:Article)
```

### Processing Pipeline

```
PubMed FTP → XML Parser → MeSH Filter → Entity Extraction →
Embedding Generation → RuVector Storage → Graph Linking
```

| Stage | Time (1M papers) | Error Rate |
|-------|------------------|------------|
| Download | 2-4 hours | <0.1% |
| XML parsing | 1-2 hours | <0.1% |
| Text cleaning | 30 min | <1% |
| Embedding (CPU) | 17 min | <0.5% |
| Embedding (GPU) | 1 min | <0.5% |
| **Total (CPU)** | **~4 hours** | **<2%** |

### Update Schedule

```bash
# Daily updates (2 AM, after PubMed daily release)
0 2 * * * /scripts/papers-daily-update.py

# Weekly maintenance (Sunday 3 AM)
0 3 * * 0 /scripts/papers-weekly-check.py
```

---

## Dependencies

### Upstream Dependencies

| Dependency | Purpose | Risk if Unavailable |
|------------|---------|---------------------|
| PubMed FTP | Primary data source | HIGH - no alternative |
| NCBI E-utilities | Real-time queries | MEDIUM - use bulk data |
| OpenAlex S3 | Coverage expansion | LOW - not critical for MVP |

### Downstream Dependents

| Dependent | Usage |
|-----------|-------|
| RAG Chat System | Paper retrieval for user queries |
| SNP Evidence | Citations for SNP-phenotype associations |
| Recommendation Engine | Evidence backing for recommendations |
| Knowledge Graph | CITED_IN relationships |

---

## Cost Summary

### One-Time Setup

| Item | Cost |
|------|------|
| PubMed baseline download | $0 (FTP) |
| Processing compute (initial) | $10-50 |
| Embedding generation | $0 (local model) |
| **Total Setup** | **~$50** |

### Ongoing Monthly

| Item | Cost |
|------|------|
| Storage (RuVector self-hosted) | $0 (included in VPS) |
| Daily updates processing | $0 (cron job) |
| API calls | $0 (bulk download) |
| **Total Monthly** | **~$0** |

### Comparison: Full-Text + Commercial Sources

| Item | Cost |
|------|------|
| DrugBank license | $10,000+/year |
| NatMed Pro | $5,000+/year |
| Storage (10x more) | $50+/month |
| Legal risk (shadow libraries) | $1B+ potential |
| **Total** | **Unaffordable** |

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| PubMed API rate limiting | Medium | Low | Use bulk FTP, respect limits |
| Storage growth | Low | Medium | Tiered compression, archival |
| Data quality issues | Medium | Medium | Validation pipeline |
| Missing recent papers | Low | Low | Daily update sync |
| Legal (public sources) | Very Low | Low | All public domain/CC0 |
| Legal (shadow libraries) | N/A | EXTREME | Not using |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `PMID` | PubMed Identifier - unique integer assigned to each PubMed citation | PMID: 12345678 |
| `DOI` | Digital Object Identifier - persistent identifier for digital content | 10.1038/s41586-021-03205-y |
| `MeSH term` | Medical Subject Heading - controlled vocabulary for indexing biomedical literature | "Polymorphism, Single Nucleotide" |
| `abstract` | Brief summary of a research paper's objectives, methods, results, and conclusions | 250-word paper summary |
| `full-text` | Complete content of a research article including all sections, figures, and tables | PDF or JATS XML article |
| `RAG` | Retrieval-Augmented Generation - combining search retrieval with LLM generation | Vector search + LLM response |
| `embedding` | Dense vector representation of text for semantic similarity search | 384-dimensional float array |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `PubMed` | NCBI database of 39M+ biomedical literature citations and abstracts | Primary literature source |
| `PMC` | PubMed Central - free full-text archive of biomedical journal articles | Full-text access |
| `MEDLINE` | NLM bibliographic database that is the primary component of PubMed | Curated citations |
| `OpenAlex` | Open catalog of 250M+ scholarly works with CC0 license | Coverage expansion |
| `Europe PMC` | European mirror of PMC with text-mined annotations | Enriched metadata |
| `E-utilities` | NCBI's Entrez Programming Utilities for programmatic PubMed access | API access |
| `JATS XML` | Journal Article Tag Suite - standard XML format for full-text articles | PMC format |
| `RuVector` | Vector database optimized for 384-dimensional embeddings | Storage backend |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PMID | PubMed Identifier | Unique citation ID |
| PMC | PubMed Central | Full-text archive |
| DOI | Digital Object Identifier | Persistent paper ID |
| MeSH | Medical Subject Headings | Controlled vocabulary |
| NLM | National Library of Medicine | PubMed maintainer |
| NCBI | National Center for Biotechnology Information | NLM division |
| RAG | Retrieval-Augmented Generation | Search + LLM pattern |
| JATS | Journal Article Tag Suite | XML standard |
| CC0 | Creative Commons Zero | Public domain license |
| CC BY | Creative Commons Attribution | Open license |
| FTP | File Transfer Protocol | Bulk download method |
| ONNX | Open Neural Network Exchange | Embedding model format |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial consolidation from research.old papers-*.md files |
