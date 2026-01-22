---
id: literature-public-sources
title: Public Research Paper Sources for Genetics/Health Knowledge Base
category: literature
tier: 1
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [literature, pubmed, pmc, europe-pmc, api, public-access]
---

# Public Research Paper Sources for Genetics/Health Knowledge Base

> **Last Updated**: January 2026
> **Purpose**: Document public sources for research papers to support RAG systems, SNP interpretations, evidence-based recommendations, and knowledge base construction for genetics/health platform.
> **Parent:** [../_index.md](../_index.md)

## Executive Summary

| Source | Total Records | Full Text | Open Access | Genetics Relevant | API Rate Limit |
|--------|---------------|-----------|-------------|-------------------|----------------|
| **PubMed** | 39+ million | ~30% linked | N/A (abstracts) | ~3-5 million | 3-10 req/sec |
| **PMC** | 10+ million | 100% | ~3.4 million | ~1-2 million | 3-10 req/sec |
| **Europe PMC** | 43+ million abstracts, 9+ million full text | ~21% | ~6.5 million | ~2-3 million | ~10 req/sec |

---

## 1. PubMed

### 1.1 Overview

**What it is**: PubMed is a free search engine accessing primarily the MEDLINE database of references and abstracts on life sciences and biomedical topics. It is maintained by the National Center for Biotechnology Information (NCBI), at the U.S. National Library of Medicine (NLM), located at the National Institutes of Health (NIH).

**Who maintains it**: U.S. National Library of Medicine (NLM) / National Center for Biotechnology Information (NCBI)

**Website**: https://pubmed.ncbi.nlm.nih.gov/

### 1.2 Content Scope

- **Total citations**: Over 39 million citations for biomedical literature
- **Annual growth**: ~1 million new citations per year (981,270 in 2022; 1,063,140 in 2021)
- **Coverage**: MEDLINE, life science journals, and online books
- **Date range**: 1809 to present (comprehensive from 1960s onward)

**Genetics-Relevant Subset Estimate**:
- MeSH terms "Genetics" and "Genomics" categories: ~3-5 million articles
- SNP-related articles: ~200,000+
- Pharmacogenomics: ~50,000+
- Autoimmune genetics: ~100,000+

### 1.3 API Access - E-utilities

**Base URL**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

#### Available Endpoints

| Endpoint | Function | Use Case |
|----------|----------|----------|
| `esearch.fcgi` | Text search, returns UIDs | Finding PMIDs for genetics terms |
| `efetch.fcgi` | Retrieve full records | Getting abstracts, metadata |
| `esummary.fcgi` | Document summaries | Quick metadata retrieval |
| `einfo.fcgi` | Database statistics | Field information |
| `elink.fcgi` | Related records | Finding linked articles |
| `epost.fcgi` | Upload UID list | Batch operations |
| `egquery.fcgi` | Global search | Cross-database queries |
| `espell.fcgi` | Spelling suggestions | Query assistance |
| `ecitmatch.fcgi` | Citation matching | Reference resolution |

#### Rate Limits

| Tier | Rate | Requirements |
|------|------|--------------|
| Anonymous | 3 requests/second | None |
| API Key | 10 requests/second | Free NCBI account + API key |
| Enhanced | >10 requests/second | Written approval from NCBI |

**Best Practices**:
- Register tool name and email in requests
- Schedule large jobs for weekends or 9 PM - 5 AM Eastern
- Use `&api_key=YOUR_KEY` parameter

#### Query Parameters

```
&db=pubmed          # Target database
&term=QUERY         # Search terms (supports Boolean: AND, OR, NOT)
&id=PMID            # Specific record ID
&rettype=abstract   # Return type (abstract, medline, xml)
&retmode=xml        # Return format (xml, text) - NO JSON for efetch
&retmax=10000       # Maximum records (default 20)
&retstart=0         # Offset for pagination
&api_key=KEY        # Authentication
```

#### Example Queries

```bash
# Search for SNP articles
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=SNP+genetics&retmax=100&api_key=YOUR_KEY

# Fetch abstract by PMID
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12345678&rettype=abstract&retmode=xml&api_key=YOUR_KEY

# Get document summary (supports JSON)
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=12345678&retmode=json&api_key=YOUR_KEY
```

### 1.4 Bulk Download Options

**FTP Locations**:
- Baseline: `https://ftp.ncbi.nlm.nih.gov/pubmed/baseline`
- Updates: `https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles`

**Baseline Files**:
- Released annually (December/January)
- 2025 baseline: files `pubmed25n0001` through `pubmed25n1274`
- Complete snapshot of all PubMed data

**Update Files**:
- Released daily
- Include new, revised, and deleted citations
- Load after baseline files in numerical order

**Update Frequency**:
- Baseline: Once per year (complete refresh)
- Updates: Daily incremental updates

### 1.5 Data Formats

**Primary Format**: XML (PubMed DTD)

**DTD Documentation**: https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_250101.dtd

**Key XML Elements**:

```xml
<PubmedArticle>
  <MedlineCitation>
    <PMID>12345678</PMID>
    <Article>
      <Journal>
        <ISSN>1234-5678</ISSN>
        <JournalIssue>
          <Volume>10</Volume>
          <Issue>2</Issue>
          <PubDate><Year>2024</Year></PubDate>
        </JournalIssue>
        <Title>Journal of Genetics</Title>
      </Journal>
      <ArticleTitle>Study on SNP rs12345</ArticleTitle>
      <Abstract>
        <AbstractText>Full abstract text...</AbstractText>
      </Abstract>
      <AuthorList>
        <Author><LastName>Smith</LastName><ForeName>John</ForeName></Author>
      </AuthorList>
    </Article>
    <MeshHeadingList>
      <MeshHeading>
        <DescriptorName>Polymorphism, Single Nucleotide</DescriptorName>
      </MeshHeading>
    </MeshHeadingList>
    <KeywordList>
      <Keyword>SNP</Keyword>
      <Keyword>genetics</Keyword>
    </KeywordList>
  </MedlineCitation>
</PubmedArticle>
```

**Available Metadata Fields**:
- PMID (unique identifier)
- Title, Abstract
- Authors (names, affiliations)
- Journal info (name, ISSN, volume, issue, pages)
- Publication date
- MeSH terms (Medical Subject Headings)
- Keywords
- Publication type
- Grant information
- DOI, PII
- Language
- References (linked PMIDs)

**JSON Support**:
- `esearch` and `esummary`: JSON supported via `retmode=json`
- `efetch`: **NO JSON support** - XML or text only

### 1.6 Licensing

**PubMed Abstracts**:
- NLM does not claim copyright on abstracts
- Journal publishers or authors may hold copyright
- Generally usable for non-commercial research with attribution

**Bulk Data License**:
- Requires NLM Data License Agreement for commercial use
- Terms documented at: `ftp.ncbi.nlm.nih.gov/pubmed/baseline/README.txt`
- Must display NCBI disclaimer and copyright notice

**Usage Terms**:
- Free for non-commercial research
- Attribution to NLM requested
- Must comply with individual publisher copyrights
- No restrictions on NCBI-generated data itself

### 1.7 Full-Text Availability

- PubMed is **abstracts only** (metadata + abstracts)
- ~30% of citations link to free full text via PMC or publishers
- Full text requires accessing PMC or publisher sites

### 1.8 Size Estimates

**Historical Size Progression**:

| Year | Records | Compressed | Uncompressed |
|------|---------|------------|--------------|
| 2016 | 24.4M | 16.9 GB | 122 GB |
| 2020 | ~32M | ~22 GB | ~160 GB |
| 2024 | ~37M | ~26 GB | ~190 GB |
| 2025 | ~39M | ~28 GB | ~200 GB |

**Genetics Subset Estimate**:
- ~10-12% of total (~4 million records)
- Compressed: ~3 GB
- Uncompressed XML: ~20-25 GB

---

## 2. PubMed Central (PMC)

### 2.1 Overview

**What it is**: PubMed Central (PMC) is a free digital archive of full-text biomedical and life sciences journal literature at the U.S. National Institutes of Health's National Library of Medicine (NIH/NLM).

**Who maintains it**: National Center for Biotechnology Information (NCBI) at NLM

**Website**: https://pmc.ncbi.nlm.nih.gov/

**Key Distinction**: PMC provides **full-text articles**, while PubMed provides abstracts/citations.

### 2.2 Content Scope

- **Total articles**: Over 10 million full-text article records
- **Date range**: Late 1700s to present
- **Content types**:
  - Formally published journal articles
  - Author manuscripts (NIH Public Access Policy)
  - Preprints

**Growth**: Doubled from 5.2 million (2018) to 10+ million (2024)

**Genetics-Relevant Subset Estimate**:
- Full-text genetics articles: ~1-2 million
- SNP-related full text: ~50,000+
- Pharmacogenomics full text: ~20,000+

### 2.3 API Access

PMC uses the same E-utilities as PubMed with `db=pmc`:

```bash
# Search PMC
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=SNP+genetics&api_key=YOUR_KEY

# Fetch full text XML
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMC1234567&rettype=xml&api_key=YOUR_KEY
```

**Additional PMC APIs**:

| Service | Purpose | URL |
|---------|---------|-----|
| OA Web Service | Open Access articles | `https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi` |
| BioC API | Text mining format | `https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/` |
| ID Converter | PMID/PMCID/DOI mapping | `https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/` |

**Rate Limits**: Same as PubMed (3-10 requests/second)

### 2.4 Bulk Download Options

**FTP Base**: `https://ftp.ncbi.nlm.nih.gov/pub/pmc/`

**Available Datasets**:

| Dataset | Description | Formats |
|---------|-------------|---------|
| **PMC Open Access Subset** | Reusable articles | XML, Text, PDF |
| **Author Manuscript Dataset** | NIH-funded manuscripts | XML, Text |
| **Historical OCR Dataset** | Older scanned articles | Text only |

**Directory Structure**:
```
/pub/pmc/
├── oa_bulk/
│   ├── oa_comm/        # Commercial use allowed (CC0, CC BY, CC BY-SA, CC BY-ND)
│   ├── oa_noncomm/     # Non-commercial only (CC BY-NC, etc.)
│   └── oa_other/       # Custom/no license
├── manuscript/
│   ├── xml/
│   └── txt/
└── historical_ocr/
    └── txt/
```

**Baseline Schedule**:
- Created at least twice per year (mid-June, mid-December)
- Additional baselines for legal article suppressions

**Cloud Access** (AWS):
- S3 URL: `s3://pmc-oa-opendata`
- HTTPS: Direct download without login
- Registry: https://registry.opendata.aws/ncbi-pmc/

### 2.5 Data Formats

**Primary Format**: JATS XML (Journal Article Tag Suite)

**JATS Standard**: NISO Z39.96-2012

**Article Structure**:
```xml
<article xmlns:xlink="http://www.w3.org/1999/xlink" article-type="research-article">
  <front>
    <journal-meta>...</journal-meta>
    <article-meta>
      <article-id pub-id-type="pmcid">PMC1234567</article-id>
      <article-id pub-id-type="pmid">12345678</article-id>
      <article-id pub-id-type="doi">10.1234/example</article-id>
      <title-group>
        <article-title>SNP Association Study</article-title>
      </title-group>
      <contrib-group>
        <contrib contrib-type="author">
          <name><surname>Smith</surname><given-names>John</given-names></name>
        </contrib>
      </contrib-group>
      <abstract><p>Abstract text...</p></abstract>
      <kwd-group>
        <kwd>SNP</kwd>
        <kwd>genetics</kwd>
      </kwd-group>
    </article-meta>
  </front>
  <body>
    <sec>
      <title>Introduction</title>
      <p>Full text content...</p>
    </sec>
    <sec>
      <title>Methods</title>
      <p>Methodology...</p>
    </sec>
    <sec>
      <title>Results</title>
      <p>Findings including SNP rs12345...</p>
    </sec>
  </body>
  <back>
    <ref-list>
      <ref id="ref1">
        <mixed-citation>Reference citation...</mixed-citation>
      </ref>
    </ref-list>
  </back>
</article>
```

**Full-Text Elements**:
- Complete article body with sections
- Tables, figures (with captions)
- Supplementary materials
- References with full citation data
- Author affiliations
- Funding information
- Conflict of interest statements

### 2.6 Licensing

**PMC Open Access Subset Categories**:

| Category | Licenses | Commercial Use | Count (est.) |
|----------|----------|----------------|--------------|
| Commercial Use Allowed | CC0, CC BY, CC BY-SA, CC BY-ND | Yes | ~1.5 million |
| Non-Commercial Only | CC BY-NC, CC BY-NC-SA, CC BY-NC-ND | No | ~1.5 million |
| Other | Custom, no license, other | Varies | ~400,000 |

**Total Open Access Subset**: ~3.4 million articles

**Key Restrictions**:
- Bulk download only via FTP or Cloud Service
- Systematic retrieval via other methods prohibited
- Must comply with individual article licenses
- PDF download limited to non-commercial licensed articles

### 2.7 Full-Text Availability

- **100% full text** (by definition - PMC is a full-text archive)
- ~34% in Open Access Subset (reusable)
- ~66% read-only (copyright restrictions)

### 2.8 Size Estimates

**PMC Open Access Subset**:
- Total articles: ~3.4 million
- XML bulk packages: ~50-100 GB compressed
- Uncompressed: ~500 GB - 1 TB (including media)

**Full PMC Archive**:
- Total articles: 10+ million
- Estimated storage: 2-5 TB uncompressed

**Genetics Subset Estimate**:
- ~15% of OA subset (~500,000 articles)
- Compressed XML: ~10-15 GB
- Uncompressed: ~100-150 GB

---

## 3. Europe PMC

### 3.1 Overview

**What it is**: Europe PMC is an open science platform providing access to life science publications and preprints. It mirrors and extends PubMed/PMC content with additional European and international sources.

**Who maintains it**: European Bioinformatics Institute (EMBL-EBI), funded by 34 European research funders

**Website**: https://europepmc.org/

**Key Distinction**: Combines PubMed abstracts, PMC full text, preprints, patents, and additional sources in one searchable platform.

### 3.2 Content Scope

- **Total abstracts**: 43+ million (includes all PubMed)
- **Full-text articles**: 9+ million (10.2 million per some sources)
- **Open access articles**: 6.5 million
- **Preprints**: 880,000+ (from 35 preprint servers)
- **PhD Theses**: 176,000+ (from UK EThOS)
- **Additional sources**: Agricola, EPO patents, NICE guidelines

**Genetics-Relevant Subset Estimate**:
- Genetics-related abstracts: ~4-5 million
- Genetics full text: ~1.5 million
- Genomics preprints: ~50,000+

### 3.3 API Access

**Base URLs**:
- Production: `https://www.ebi.ac.uk/europepmc/webservices/rest/`
- Test: `https://www.ebi.ac.uk/europepmc/webservices/test/rest/`

**No authentication required** - fully open API

#### Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/search` | GET/POST | Query publications |
| `/fields` | GET | Available search fields |
| `/{source}/{id}/citations` | GET | Citation data |
| `/{source}/{id}/references` | GET | Reference lists |
| `/{source}/{id}/databaseLinks` | GET | Database cross-references |
| `/{source}/{id}/textMinedTerms` | GET | Text-mined annotations |
| `/{id}/fullTextXML` | GET | Full-text XML |
| `/{id}/supplementaryFiles` | GET | Supplementary materials |

#### Query Parameters

```
?query=SEARCH_TERM      # Search query (max 1500 chars)
?format=json|xml|dc     # Response format
?resultType=idlist|lite|core  # Detail level
?cursorMark=*           # Pagination cursor
?pageSize=25            # Results per page (max 1000)
?sort=relevance|date|cited  # Sort order
```

#### Rate Limits

- **Standard**: ~10 requests/second
- **No API key required**
- Contact epmc-webservices@ebi.ac.uk for higher limits

#### Example Queries

```bash
# Search for SNP articles
https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=SNP%20genetics&format=json&resultType=core

# Get full text XML
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1234567/fullTextXML

# Get citations
https://www.ebi.ac.uk/europepmc/webservices/rest/MED/12345678/citations?format=json

# Get text-mined terms (genes, diseases, chemicals)
https://www.ebi.ac.uk/europepmc/webservices/rest/MED/12345678/textMinedTerms?format=json
```

### 3.4 Bulk Download Options

**FTP Base**: `https://europepmc.org/ftp/`

**Available Downloads**:

| Dataset | Description | Format | Update |
|---------|-------------|--------|--------|
| Open Access Subset | Full text + supplementary | XML, PDF | Weekly |
| Author Manuscripts | Funder-deposited manuscripts | XML, Text | Weekly |
| Preprints Subset | Open access preprints | XML | Weekly |
| Metadata | All article metadata | XML | Weekly |
| PMID-PMCID-DOI Mappings | ID crosswalk | CSV | Monthly |
| Accession Numbers | Text-mined database IDs | CSV | Weekly |

**Open Access Subset Details**:
- ~5.7 million articles with CC licenses
- Organized in files of 10K articles each
- PDF bulk download available (in addition to XML)
- Quarterly archives

### 3.5 Data Formats

**Response Formats**:

| Format | Description | Use Case |
|--------|-------------|----------|
| XML | Default, matches SOAP service | Full compatibility |
| JSON | Structured data | Modern applications |
| DC (Dublin Core) | RDF/XML metadata | Linked data |

**Result Types**:

| Type | Content | Fields |
|------|---------|--------|
| `idlist` | IDs only | id, source |
| `lite` | Basic metadata | title, authors, journal, date |
| `core` | Full metadata | + abstract, MeSH, links, annotations |

**Unique Features**:

- **Text-mined annotations**: Genes, proteins, chemicals, diseases, organisms
- **Database cross-references**: UniProt, ENA, PDB, ChEMBL, etc.
- **Citation network**: 19.4+ million reference lists
- **Preprint-publication links**: Tracking preprint to peer-review

### 3.6 Licensing

**Content Licensing**:
- All content free to read
- Open Access subset: CC licenses (reusable)
- Non-OA content: Subject to publisher copyright

**API Usage**:
- Free, no registration required
- Must accept EBI Privacy Notice
- Attribution appreciated

**Bulk Download**:
- Open Access subset freely downloadable
- Non-OA content restricted to API access

### 3.7 Full-Text Availability

| Metric | Count | Percentage |
|--------|-------|------------|
| Total abstracts | 43+ million | 100% |
| Full text | 9-10 million | ~21% |
| Open access | 6.5 million | ~15% |
| Reusable (text mining) | 2.5+ million | ~6% |

**Access via Unpaywall**: Links to 13+ million additional free full-text articles on publisher sites.

### 3.8 Size Estimates

**Full Database**:
- Abstracts/metadata: ~50 GB compressed
- Full text XML: ~200-300 GB compressed
- PDFs: ~500 GB+

**Open Access Subset**:
- ~5.7 million articles
- XML: ~150 GB compressed
- PDFs: ~300 GB

**Genetics Subset Estimate**:
- ~15% of content
- Full text XML: ~30-45 GB compressed
- Abstracts: ~7-8 GB compressed

---

## Comparison for Genetics/Health Knowledge Base

### Content Overlap and Uniqueness

```
┌─────────────────────────────────────────────────────────────┐
│                     Europe PMC (43M+ abstracts)             │
│  ┌───────────────────────────────────────────────────────┐  │
│  │              PubMed (39M+ abstracts)                  │  │
│  │  ┌─────────────────────────────────────────────────┐  │  │
│  │  │         PMC (10M+ full text)                    │  │  │
│  │  │  ┌───────────────────────────────────────────┐  │  │  │
│  │  │  │    Open Access (3.4M reusable)            │  │  │  │
│  │  │  └───────────────────────────────────────────┘  │  │  │
│  │  └─────────────────────────────────────────────────┘  │  │
│  └───────────────────────────────────────────────────────┘  │
│  + Preprints (880K) + Patents + NICE + Agricola + Theses    │
└─────────────────────────────────────────────────────────────┘
```

### Recommended Strategy by Use Case

| Use Case | Primary Source | Secondary | Notes |
|----------|----------------|-----------|-------|
| **RAG for AI Chat** | Europe PMC API | PMC OA Subset | Use text-mined annotations |
| **SNP Evidence Linking** | PubMed (comprehensive) | PMC (full context) | Search by rs numbers |
| **Recommendation Citations** | PubMed (abstracts) | Europe PMC (full text) | Include DOIs |
| **Knowledge Base Construction** | PMC OA Subset (bulk) | Europe PMC (enriched) | Use JATS XML |

### API Comparison

| Feature | PubMed E-utils | PMC | Europe PMC |
|---------|----------------|-----|------------|
| Authentication | API key (optional) | API key (optional) | None required |
| Rate limit | 3-10 req/sec | 3-10 req/sec | ~10 req/sec |
| JSON support | Partial (search only) | Partial | Full |
| Full text | No | Yes | Yes |
| Text mining | No | Limited | Extensive |
| Cross-references | PMC links | PubMed links | UniProt, ENA, etc. |

### Storage Requirements Summary

**Minimum (Genetics-focused)**:
- PubMed genetics abstracts: ~3 GB compressed
- PMC OA genetics full text: ~15 GB compressed
- Europe PMC annotations: ~5 GB compressed
- **Total**: ~25 GB compressed, ~250 GB uncompressed

**Comprehensive**:
- Full PubMed baseline: ~28 GB compressed
- Full PMC OA subset: ~100 GB compressed
- Europe PMC enrichments: ~50 GB compressed
- **Total**: ~200 GB compressed, ~2 TB uncompressed

---

## Implementation Recommendations

### Phase 1: Initial Data Collection

1. **Europe PMC API** for real-time search and citations
   - No auth required, immediate access
   - Rich annotations for genetics terms
   - JSON responses for easy parsing

2. **PubMed E-utilities** for comprehensive PMID coverage
   - Register for API key (10 req/sec)
   - Use ESummary with JSON for metadata
   - Use EFetch with XML for abstracts

### Phase 2: Bulk Data Pipeline

1. **PMC Open Access Subset** via FTP/AWS
   - Download commercial-use allowed subset first
   - Parse JATS XML for full text
   - Extract supplementary data files

2. **PubMed Baseline** for abstract corpus
   - Annual baseline + daily updates
   - Parse XML for structured metadata
   - Build PMID-to-content mapping

### Phase 3: Enrichment

1. **Europe PMC Text Mining**
   - Gene/protein mentions
   - Disease annotations
   - Chemical compounds
   - Database accession numbers

2. **Cross-Reference Integration**
   - PMID-PMCID-DOI mappings
   - UniProt, ENA, dbSNP links
   - Citation networks

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `E-utilities` | NCBI's suite of server-side programs for searching and retrieving Entrez data | esearch, efetch, einfo |
| `DTD` | Document Type Definition - schema defining XML document structure | pubmed_250101.dtd |
| `JATS XML` | Journal Article Tag Suite - standardized XML format for full-text articles | PMC article format |
| `rate limit` | Maximum number of API requests allowed per time unit | 10 requests/second with API key |
| `API key` | Authentication token enabling higher rate limits and tracking | NCBI API key for 10 req/sec |
| `bulk download` | Retrieving large datasets via FTP rather than individual API calls | PubMed baseline files |
| `text mining` | Automated extraction of structured information from unstructured text | Gene/disease entity extraction |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `esearch` | E-utility endpoint for text searches returning UIDs | Finding PMIDs |
| `efetch` | E-utility endpoint for retrieving full records | Getting abstracts |
| `esummary` | E-utility endpoint for document summaries | Quick metadata |
| `Open Access Subset` | PMC articles with reusable CC licenses (~3.4M articles) | Commercial use |
| `BioC` | Text mining interchange format for sharing annotated text | NLP pipelines |
| `text-mined terms` | Entities automatically extracted from article text | Genes, diseases, chemicals |
| `oa_comm` | PMC directory containing commercially usable articles | CC0, CC BY, CC BY-SA |
| `oa_noncomm` | PMC directory containing non-commercial only articles | CC BY-NC licenses |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NCBI | National Center for Biotechnology Information | API provider |
| NLM | National Library of Medicine | PubMed host |
| UID | Unique Identifier | Generic ID (PMID, PMCID) |
| DTD | Document Type Definition | XML schema |
| JATS | Journal Article Tag Suite | Full-text XML format |
| OA | Open Access | Freely available content |
| CC | Creative Commons | License type |
| API | Application Programming Interface | Programmatic access |
| FTP | File Transfer Protocol | Bulk download method |
| REST | Representational State Transfer | API architecture style |
| XML | Extensible Markup Language | Data format |
| JSON | JavaScript Object Notation | Data format |

---

## References

- [NCBI E-utilities Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25497/)
- [PubMed Download Data](https://pubmed.ncbi.nlm.nih.gov/download/)
- [PMC FTP Service](https://pmc.ncbi.nlm.nih.gov/tools/ftp/)
- [PMC Open Access Subset](https://pmc.ncbi.nlm.nih.gov/tools/openftlist/)
- [Europe PMC REST API](https://europepmc.org/RestfulWebService)
- [Europe PMC Bulk Downloads](https://europepmc.org/downloads)
- [JATS Standard (NISO)](https://www.niso.org/standards-committees/jats)
- [NCBI Usage Policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- [NLM Data License](https://www.nlm.nih.gov/databases/download.html)
- [Europe PMC 2024 Year in Review](https://blog.europepmc.org/2025/01/europe-pmc-2024-a-year-in-review.html)
