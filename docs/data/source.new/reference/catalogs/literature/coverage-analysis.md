---
id: literature-coverage-analysis
title: Literature Coverage Analysis for Genetics Knowledge Base
category: literature
tier: 2
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [literature, coverage, gap-analysis, international, pubmed]
---

# Literature Coverage Analysis for Genetics Knowledge Base

> **Last Updated**: January 2026
> **Purpose**: Assess the coverage gaps in PubMed and other major literature databases for genetics/health topics to inform multi-source data strategy.
> **Parent:** [../_index.md](../_index.md)

---

## Executive Summary

PubMed provides **excellent coverage (90-95%)** of SNP/GWAS studies, pharmacogenomics, and core genetics research in English-language journals. Major gaps exist for:

- **Alternative medicine research** (60-70% missing): Supplements, TCM, Ayurveda not indexed
- **Non-English journals** (50-60% missing): European, Asian, South American journals underrepresented
- **Preprints** (95% missing): PubMed doesn't index preprints until peer-reviewed publication
- **Conference papers** (80% missing): Most conference proceedings not indexed
- **Gray literature** (90% missing): Technical reports, theses, government documents

**Recommendation**: Use OpenAlex + Europe PMC for gap-filling (~10-20% additional coverage), especially for:
- International genetics research
- Preprints (early access to emerging findings)
- Supplements/nutrition research
- Traditional medicine

---

## 1. PubMed Coverage Strengths

### 1.1 High Coverage Areas (90-95%)

| Topic | Coverage | Total Papers (est.) | PubMed Count (est.) |
|-------|----------|---------------------|---------------------|
| **SNP/GWAS studies** | 95% | 300,000 | ~285,000 |
| **Pharmacogenomics (core)** | 92% | 60,000 | ~55,000 |
| **Clinical genetics** | 90% | 500,000 | ~450,000 |
| **Human disease genetics** | 93% | 800,000 | ~750,000 |
| **Gene expression** | 88% | 1,200,000 | ~1,050,000 |
| **Molecular diagnostics** | 85% | 200,000 | ~170,000 |

**Why PubMed Excels**:
- NIH-funded research mandatory deposit
- Major journals indexed: Nature Genetics, NEJM, JAMA, etc.
- Strong U.S./Western Europe coverage
- Comprehensive MeSH indexing for genetics terms

### 1.2 Medium Coverage Areas (70-85%)

| Topic | Coverage | Gap Source |
|-------|----------|------------|
| **Cancer genomics** | 82% | Specialist conference proceedings not indexed |
| **Microbiome genetics** | 75% | Emerging field, many preprints |
| **Epigenetics** | 78% | Preprints, non-traditional journals |
| **Behavioral genetics** | 70% | Overlaps with psychology journals (less indexed) |
| **Agricultural genetics** | 65% | Many papers in Agricola, not PubMed |

---

## 2. PubMed Coverage Gaps

### 2.1 Alternative Medicine & Supplements (40-70% Missing)

| Sub-field | PubMed Coverage | Missing Source |
|-----------|-----------------|----------------|
| **Dietary supplements** | 60% | Industry-funded research, trade journals |
| **Herbal medicine** | 55% | Regional journals (Asia, South America) |
| **Traditional Chinese Medicine (TCM)** | 30% | CNKI (China), not indexed |
| **Ayurveda** | 25% | Indian journals (many not indexed) |
| **Homeopathy** | 40% | Specialized journals |
| **Integrative medicine** | 50% | Many non-MEDLINE journals |

**Example**: A search for "MTHFR methylation folate" returns 5,000+ PubMed results. But **folate supplements** research has significant coverage outside PubMed:
- Trade journals (e.g., Natural Products Insider)
- Practitioner journals (e.g., Townsend Letter)
- Regional databases (e.g., IndMED for Ayurveda)

**Impact on Gene KB**:
- SNP interpretations for "methylation support" recommendations may miss supplement research
- Alternative protocol evidence scattered across non-indexed sources

**Gap-Filling Sources**:

| Source | Coverage Added | URL |
|--------|----------------|-----|
| **AMED (Allied & Complementary Medicine)** | TCM, Ayurveda, homeopathy | https://www.bl.uk/collection-guides/amed |
| **IndMED** | Indian biomedical research | https://indmed.nic.in/ |
| **CNKI (China)** | Chinese medical research | https://www.cnki.net/ (paid) |
| **NAPRALERT** | Natural products research | https://www.napralert.org/ |

### 2.2 Non-English Literature (30-50% Missing)

| Region | PubMed Coverage | Major Gaps |
|--------|-----------------|------------|
| **Western Europe** | 85% | German, French journals underrepresented |
| **Eastern Europe** | 50% | Russian, Polish journals |
| **Asia** | 45% | Chinese, Japanese, Korean journals |
| **South America** | 40% | Spanish, Portuguese journals |
| **Africa** | 35% | Regional health journals |

**Example**: A study on genetic variations in Asian populations may be published in:
- Japanese Journal of Human Genetics (partially indexed)
- Chinese Journal of Medical Genetics (not indexed in PubMed)
- Korean Journal of Genetics (not indexed)

**Gap-Filling Sources**:

| Database | Region | Free Access | URL |
|----------|--------|-------------|-----|
| **SciELO** | Latin America | Yes | https://scielo.org/ |
| **LILACS** | Latin America | Yes | https://lilacs.bvsalud.org/ |
| **J-STAGE** | Japan | Yes | https://www.jstage.jst.go.jp/ |
| **KoreaMed** | South Korea | Yes | https://koreamed.org/ |
| **African Index Medicus** | Africa | Yes | https://indexmedicus.afro.who.int/ |
| **IndMED** | India | Yes | https://indmed.nic.in/ |

### 2.3 Preprints (95% Missing)

PubMed does NOT index preprints until they're peer-reviewed and published. This creates a **6-12 month lag** for emerging research.

| Preprint Server | Papers | Genetics Relevant |
|----------------|--------|-------------------|
| **bioRxiv** | 200,000+ | ~50,000 |
| **medRxiv** | 50,000+ | ~10,000 |
| **arXiv (q-bio)** | 30,000+ | ~5,000 |
| **Research Square** | 100,000+ | ~15,000 |

**Why This Matters**:
- GWAS findings often appear on bioRxiv first
- COVID-19 genomics research primarily published as preprints
- Pharmacogenomics updates faster via preprints

**Gap-Filling Sources**:

| Source | Content | URL |
|--------|---------|-----|
| **Europe PMC** | 880,000+ preprints from 35 servers | https://europepmc.org/ |
| **OpenAlex** | Preprint tracking + publication linking | https://openalex.org/ |
| **bioRxiv API** | Direct preprint access | https://api.biorxiv.org/ |

### 2.4 Conference Papers & Proceedings (80% Missing)

| Conference Type | PubMed Indexing |
|----------------|-----------------|
| Major conferences (ASHG, ASCO) | ~30% (abstracts only) |
| Regional conferences | ~5% |
| Workshop proceedings | ~0% |
| Poster presentations | ~0% |

**Example**: American Society of Human Genetics (ASHG) annual meeting abstracts contain novel SNP findings, but only published abstracts (not posters) are indexed.

**Gap-Filling Sources**:

| Source | Coverage |
|--------|----------|
| **OpenAlex** | Conference papers + preprints |
| **Semantic Scholar** | Conference proceedings mining |

### 2.5 Gray Literature (90% Missing)

| Content Type | PubMed Coverage |
|--------------|-----------------|
| PhD theses | ~0% (except published papers from theses) |
| Technical reports | ~5% |
| Government reports | ~10% |
| Industry white papers | ~0% |
| Clinical trial data (unpublished) | ~0% |

**Impact on Gene KB**:
- Negative findings often in theses/reports, not published
- Regulatory data (FDA, EMA) not indexed
- Industry pharmacogenomics data unpublished

**Gap-Filling Sources**:

| Source | Content | URL |
|--------|---------|-----|
| **Europe PMC** | 176,000+ UK theses (EThOS) | https://europepmc.org/ |
| **ProQuest Dissertations** | Global theses | https://www.proquest.com/ (paid) |
| **ClinicalTrials.gov** | Trial results (some unpublished) | https://clinicaltrials.gov/ |
| **OpenGrey** | European gray literature | https://www.opengrey.eu/ |

---

## 3. Multi-Source Coverage Strategy

### 3.1 Layered Approach

```
┌───────────────────────────────────────────────────────────────┐
│ LAYER 1: PubMed/MEDLINE (Primary - 90% core genetics)        │
│ - 39M biomedical papers                                       │
│ - Best MeSH term filtering                                    │
│ - Strong SNP/GWAS/pharmacogenomics coverage                   │
└───────────────────────────────────────────────────────────────┘
                              ↓
┌───────────────────────────────────────────────────────────────┐
│ LAYER 2: OpenAlex (Breadth - Fill International/Preprint Gaps│
│ - 250M works (6x PubMed)                                      │
│ - 40% preprint coverage                                       │
│ - Better non-English journal coverage                         │
│ - Conference papers included                                  │
└───────────────────────────────────────────────────────────────┘
                              ↓
┌───────────────────────────────────────────────────────────────┐
│ LAYER 3: Europe PMC (Enrichment - Text Mining + Preprints)   │
│ - 43M abstracts (includes PubMed)                             │
│ - 880K preprints                                              │
│ - Text-mined gene/protein/disease annotations                 │
└───────────────────────────────────────────────────────────────┘
                              ↓
┌───────────────────────────────────────────────────────────────┐
│ LAYER 4: Specialized (Domain-Specific Gaps)                  │
│ - SciELO (Latin America)                                      │
│ - IndMED (India/Ayurveda)                                     │
│ - CNKI (China/TCM)                                            │
│ - AMED (Alternative medicine)                                 │
└───────────────────────────────────────────────────────────────┘
```

### 3.2 Estimated Coverage Gains

| Source Combination | Genetics Papers | Coverage |
|--------------------|-----------------|----------|
| **PubMed only** | 3-4M | 85% |
| **PubMed + OpenAlex** | 4-5M | 92% |
| **PubMed + OpenAlex + Europe PMC** | 5-6M | 95% |
| **+ Specialized (CNKI, IndMED, etc.)** | 6-7M | 98% |

### 3.3 Cost-Benefit Analysis

| Source | Setup Cost | Annual Cost | Coverage Added |
|--------|------------|-------------|----------------|
| **PubMed** | $0 | $0 | 85% (baseline) |
| **OpenAlex** | $0 | $0 | +7% |
| **Europe PMC** | $0 | $0 | +3% |
| **SciELO** | $0 | $0 | +1% (Latin America) |
| **IndMED** | $0 | $0 | +0.5% (India) |
| **CNKI** | $5,000 | $5,000 | +2% (China) |
| **Semantic Scholar** | $0 | Request | +1% (AI analysis) |

**Recommendation**: For a bootstrapped startup, use PubMed + OpenAlex + Europe PMC (all free) for **95% coverage**. Defer paid sources (CNKI) until revenue stage.

---

## 4. Coverage by Genetics Sub-Domain

### 4.1 Pharmacogenomics

| Source | Coverage | Notes |
|--------|----------|-------|
| PubMed | 90% | Strong PharmGKB-linked papers |
| OpenAlex | +5% | International pharmacogenomics societies |
| Europe PMC | +2% | Preprints for emerging drugs |
| **Gap**: Industry data | ~20% of pharmacogenomics data never published |

### 4.2 SNP-Disease Associations

| Source | Coverage | Notes |
|--------|----------|-------|
| PubMed | 95% | GWAS catalog sources well-indexed |
| OpenAlex | +2% | Replication studies in regional journals |
| Europe PMC | +1% | Preprints for novel associations |
| **Gap**: Negative findings | Many non-significant SNP studies unpublished |

### 4.3 Nutrigenomics

| Source | Coverage | Notes |
|--------|----------|-------|
| PubMed | 70% | Nutrition journals less indexed |
| OpenAlex | +10% | Food science journals included |
| Europe PMC | +5% | Preprints from nutrition conferences |
| AMED | +10% | Supplement research |
| **Gap**: Industry | Supplement company research proprietary |

### 4.4 Epigenetics

| Source | Coverage | Notes |
|--------|----------|-------|
| PubMed | 78% | Emerging field, growing fast |
| bioRxiv/medRxiv | +15% | Many preprints not yet published |
| Europe PMC | +5% | Text-mined DNA methylation data |
| **Gap**: Methodological papers in specialized journals |

---

## 5. Implementation Priorities

### Phase 1: MVP (Weeks 1-4)

**Focus**: PubMed genetics subset only

- **Coverage**: 85% of core genetics literature
- **Cost**: $0
- **Effort**: Low (established pipeline)

### Phase 2: Breadth (Months 2-3)

**Add**: OpenAlex

- **Coverage**: +7% (92% total)
- **Cost**: $0
- **Effort**: Medium (API integration)
- **Value**: Preprints, international journals, conference papers

### Phase 3: Enrichment (Months 4-6)

**Add**: Europe PMC text mining

- **Coverage**: +3% (95% total)
- **Cost**: $0
- **Effort**: Medium (text-mined annotations)
- **Value**: Gene/protein/disease entity extraction

### Phase 4: Specialization (Year 2)

**Add**: Targeted gap-filling

- SciELO (Latin America)
- IndMED (India/Ayurveda)
- AMED (Alternative medicine)

**Coverage**: +3% (98% total)
**Cost**: $0 (all free)
**Effort**: High (disparate APIs)

---

## 6. Deduplication Strategy

### 6.1 Cross-Database Overlap

| Database Pair | Overlap |
|---------------|---------|
| PubMed ↔ PMC | ~30% (PMC is subset) |
| PubMed ↔ Europe PMC | ~95% (Europe PMC mirrors PubMed) |
| PubMed ↔ OpenAlex | ~95% (OpenAlex indexes PubMed) |
| PMC ↔ Europe PMC | ~100% (identical full-text sources) |

### 6.2 Deduplication Keys

**Primary**: DOI (Digital Object Identifier)
- Most reliable unique identifier
- ~80% of papers have DOIs

**Secondary**: PMID (PubMed ID)
- Unique within PubMed/PMC/Europe PMC
- ~90% of relevant papers have PMIDs

**Tertiary**: Title + First Author + Year
- For papers without DOI/PMID
- Requires normalization (case, punctuation, spacing)

### 6.3 Deduplication Pipeline

```
┌─────────────────────────────────────────────────────────────┐
│ Step 1: Normalize Identifiers                              │
│ - Lowercase DOIs                                            │
│ - Remove "doi:" prefix                                      │
│ - Standardize PMID format                                   │
└─────────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 2: Match by DOI                                        │
│ - Hash DOI → lookup table                                   │
│ - Match rate: ~75-80%                                       │
└─────────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 3: Match by PMID (for papers without DOI)             │
│ - Hash PMID → lookup table                                  │
│ - Match rate: ~10-15% additional                            │
└─────────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 4: Fuzzy Match by Title + Author + Year               │
│ - Levenshtein distance < 5 characters                       │
│ - Match rate: ~5% additional                                │
└─────────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ Result: 90-95% deduplication accuracy                       │
└─────────────────────────────────────────────────────────────┘
```

---

## 7. Quality Indicators

### 7.1 Source Reliability

| Source | Data Quality | Update Frequency | Reliability |
|--------|--------------|------------------|-------------|
| **PubMed** | Excellent | Daily | Very High |
| **PMC** | Excellent | Daily | Very High |
| **Europe PMC** | Excellent | Weekly | Very High |
| **OpenAlex** | Good | Daily | High |
| **bioRxiv** | Variable | Real-time | Medium (not peer-reviewed) |
| **Preprints** | Variable | Real-time | Medium |
| **CNKI** | Good | Weekly | Medium (language barrier) |

### 7.2 Metadata Completeness

| Source | Title | Abstract | Full Text | MeSH | Authors | DOI | PMID |
|--------|-------|----------|-----------|------|---------|-----|------|
| **PubMed** | 100% | 95% | 30% | 90% | 99% | 85% | 100% |
| **PMC** | 100% | 100% | 100% | 90% | 100% | 90% | 100% |
| **Europe PMC** | 100% | 98% | 21% | 85% | 99% | 88% | 95% |
| **OpenAlex** | 100% | 70% | 15% | 0% | 95% | 80% | 75% |
| **bioRxiv** | 100% | 100% | 100% | 0% | 100% | 100% | 0% |

---

## 8. Monitoring and Continuous Improvement

### 8.1 Coverage Metrics to Track

| Metric | Target | Measurement |
|--------|--------|-------------|
| **Genetics paper count** | 5-7M | Quarterly audit |
| **SNP mention coverage** | 95% | Sample 1000 rs numbers, check citations |
| **Drug-gene interaction coverage** | 90% | PharmGKB cross-reference |
| **Recent paper lag** | <7 days | Daily update monitoring |
| **Deduplication accuracy** | 95% | Manual review of 500 random samples |

### 8.2 User Feedback Loop

Track user queries that return 0 results:
- May indicate coverage gaps
- Guide addition of specialized sources
- Identify emerging topics not yet indexed

---

## Download

| Source | Method | URL/Command |
|--------|--------|-------------|
| **PubMed** | FTP baseline | `ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/` |
| **OpenAlex** | API/Snapshot | `https://docs.openalex.org/download-all-data` |
| **Europe PMC** | FTP | `https://europepmc.org/downloads` |
| **SciELO** | OAI-PMH | `https://www.scielo.org/` |

**Access Requirements:** Most sources are freely accessible; some regional databases may require registration.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | XML, JSON, CSV |
| Alternative | Parquet, TSV |
| Identifiers | PMID, DOI, OpenAlex ID |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `source_db` | string | Database origin | "pubmed" |
| `coverage_pct` | float | Coverage percentage | 0.95 |
| `gap_type` | string | Type of coverage gap | "regional" |
| `last_audit` | date | Audit date | "2026-01-15" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `supplements` | Source | N:M |
| `has_gap` | Gap Type | 1:N |

## Sample Data

### Example Coverage Record
```json
{
  "topic": "GWAS",
  "pubmed_coverage": 0.95,
  "openalex_coverage": 0.98,
  "gaps_identified": ["non-English publications", "preprints"],
  "recommendation": "Add bioRxiv/medRxiv"
}
```

### Sample Query Result
| topic | pubmed | openalex | gap |
|-------|--------|----------|-----|
| GWAS | 95% | 98% | Preprints |
| TCM | 45% | 70% | Chinese sources |

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| PubMed | Public domain | Yes |
| OpenAlex | CC0 | Yes |
| Europe PMC | CC BY | Yes |

## Data Set Size

| Metric | Value |
|--------|-------|
| PubMed genetics papers | ~5-7M |
| OpenAlex works | 250M+ |
| Coverage gap records | ~500 documented |
| Audit frequency | Quarterly |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `coverage` | The proportion of relevant literature captured by a database | 95% coverage of GWAS papers |
| `coverage gap` | Literature systematically missing from a database | TCM research not in PubMed |
| `preprint` | Research paper shared before formal peer review | bioRxiv preprint |
| `gray literature` | Publications outside commercial/academic publishing channels | PhD theses, technical reports |
| `deduplication` | Process of identifying and removing duplicate records across sources | DOI-based matching |
| `indexing lag` | Time between publication and database indexing | 6-12 months for preprints |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `MEDLINE` | NLM's curated bibliographic database indexed with MeSH terms | Core PubMed content |
| `OpenAlex` | Open catalog of 250M+ works filling PubMed gaps | Coverage expansion |
| `SciELO` | Scientific Electronic Library Online - Latin American literature | Regional coverage |
| `CNKI` | China National Knowledge Infrastructure - Chinese academic database | TCM research |
| `IndMED` | Indian biomedical literature database | Ayurveda research |
| `AMED` | Allied and Complementary Medicine Database | Alternative medicine |
| `bioRxiv` | Preprint server for biology research | Early access |
| `medRxiv` | Preprint server for health sciences | Clinical preprints |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GWAS | Genome-Wide Association Study | SNP-disease discovery |
| TCM | Traditional Chinese Medicine | Major coverage gap area |
| CNKI | China National Knowledge Infrastructure | Chinese literature |
| SciELO | Scientific Electronic Library Online | Latin America |
| AMED | Allied and Complementary Medicine | CAM database |
| CAM | Complementary and Alternative Medicine | Coverage gap category |
| NIH | National Institutes of Health | US research funder |
| NLM | National Library of Medicine | PubMed maintainer |
| EThOS | Electronic Theses Online Service | UK thesis repository |
| DOI | Digital Object Identifier | Deduplication key |
| PMID | PubMed Identifier | Deduplication key |

---

## References

- [PubMed Journal Selection](https://www.nlm.nih.gov/medline/medline_overview.html)
- [OpenAlex Coverage Analysis](https://docs.openalex.org/about-the-data/coverage)
- [Europe PMC Content Sources](https://europepmc.org/About)
- [SciELO Network](https://scielo.org/en/about-scielo/scielo-network/)
- [Preprint Growth Statistics (ASAPbio)](https://asapbio.org/preprint-info/preprint-statistics)

---

*This document is part of the Gene Knowledge Base technical documentation.*
