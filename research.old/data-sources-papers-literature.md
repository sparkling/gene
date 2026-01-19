# Scientific Paper and Literature Data Sources for Gene Platform

> **Last Updated**: January 2026
> **Purpose**: Comprehensive reference for scientific literature APIs and databases for full-text access, open access detection, and preprint repositories.

---

## Executive Summary

| Source | Total Records | Full-Text | Open Access | API Rate Limit | Auth Required |
|--------|---------------|-----------|-------------|----------------|---------------|
| **PubMed Central (PMC)** | 11+ million | 100% | ~3.4M (OA Subset) | 3-10 req/sec | API key optional |
| **Europe PMC** | 33+ million | 10.2M | 6.5M | ~10 req/sec | None |
| **Semantic Scholar** | 214+ million | Links only | Varies | 1-100 req/sec | API key optional |
| **OpenAlex** | 271+ million | Links only | Tracks OA status | 100K/day | API key free |
| **CORE** | 391+ million | 38M+ hosted | 391M+ | Tiered | API key optional |
| **Unpaywall** | 95M+ DOIs | Links to OA | ~24-27M | 100K/day | Email required |
| **bioRxiv/medRxiv** | 250K+ | 100% | 100% | No limit | None |

---

## 1. PubMed Central (PMC) Open Access Subset

### Overview

**Description**: PMC is a free full-text archive of biomedical and life sciences journal literature maintained by the NIH/NLM. The Open Access Subset contains articles that allow reuse.

**Maintainer**: National Center for Biotechnology Information (NCBI) / National Library of Medicine (NLM)

**Website**: https://pmc.ncbi.nlm.nih.gov/

### Content Statistics

| Metric | Value |
|--------|-------|
| Total Articles | 11+ million full-text records |
| Open Access Subset | ~3.4 million articles |
| Commercial Use Allowed | ~1.5 million (CC0, CC BY, CC BY-SA, CC BY-ND) |
| Non-Commercial Only | ~1.5 million (CC BY-NC variants) |
| Other Licenses | ~400,000 |
| Date Range | Late 1700s to present |
| Annual Growth | Doubled from 5.2M (2018) to 11M (2025) |

### API Access

**Base URL**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

#### Available APIs

| API | URL | Purpose |
|-----|-----|---------|
| E-utilities | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/` | Search and fetch |
| OA Web Service | `https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi` | Discover OA articles |
| BioC API | `https://www.ncbi.nlm.nih.gov/research/bionlp/APIs/BioC-PMC/` | Text mining format |
| OAI-PMH | `https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi` | Metadata harvesting |
| ID Converter | `https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/` | PMID/PMCID/DOI mapping |

#### Rate Limits

| Tier | Rate | Requirements |
|------|------|--------------|
| Anonymous | 3 requests/second | None |
| API Key | 10 requests/second | Free NCBI account |
| Enhanced | >10 requests/second | Written approval |

#### Example API Calls

```bash
# Search PMC for genetics articles
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=SNP+genetics&retmax=100&api_key=YOUR_KEY

# Fetch full-text XML by PMC ID
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMC1234567&rettype=xml

# OA Web Service - find downloadable OA articles
https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC1234567

# BioC API - get article in text mining format
https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC1234567/unicode
```

### Bulk Download

**FTP Base**: `https://ftp.ncbi.nlm.nih.gov/pub/pmc/`

**AWS S3**: `s3://pmc-oa-opendata` (no login required)

| Package Type | Location | Content |
|--------------|----------|---------|
| Commercial Use | `/oa_bulk/oa_comm/` | CC0, CC BY, CC BY-SA, CC BY-ND |
| Non-Commercial | `/oa_bulk/oa_noncomm/` | CC BY-NC variants |
| Other | `/oa_bulk/oa_other/` | Custom/no license |

### Data Formats

- **Primary**: JATS XML (Journal Article Tag Suite)
- **Alternative**: Plain text, PDF (for some articles)
- **Text Mining**: BioC XML/JSON

### Licensing

| License Category | Commercial Use | Text Mining | Examples |
|------------------|----------------|-------------|----------|
| CC0 | Yes | Yes | Public domain |
| CC BY | Yes | Yes | Attribution only |
| CC BY-NC | No | Research only | Non-commercial |
| Other | Varies | Check terms | Custom licenses |

### Size Estimates

- Compressed XML bulk: ~50-100 GB
- Uncompressed with media: ~500 GB - 1 TB
- Genetics subset: ~10-15 GB compressed

### Full-Text Availability

**100%** - PMC is a full-text archive by definition. All articles include complete body text, figures, tables, and references.

---

## 2. Europe PMC

### Overview

**Description**: Open science platform providing access to worldwide life science publications, preprints, patents, and clinical guidelines. Mirrors and extends PubMed/PMC with European sources.

**Maintainer**: European Bioinformatics Institute (EMBL-EBI), funded by 34 European research funders

**Website**: https://europepmc.org/

### Content Statistics

| Metric | Value |
|--------|-------|
| Total Publications | 33+ million |
| Full-Text Articles | 10.2 million |
| Open Access | 6.5 million |
| Reference Lists | 19.4+ million publications |
| Preprints | 880,000+ (from 35 servers) |
| PhD Theses | 176,000+ (UK EThOS) |
| Patents | European Patent Office (EPO) |

### API Access

**Base URL**: `https://www.ebi.ac.uk/europepmc/webservices/rest/`

**Authentication**: None required

#### Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/search` | GET/POST | Search publications |
| `/fields` | GET | Available search fields |
| `/{source}/{id}/citations` | GET | Get citing articles |
| `/{source}/{id}/references` | GET | Get reference list |
| `/{source}/{id}/databaseLinks` | GET | Cross-references |
| `/{source}/{id}/textMinedTerms` | GET | Text-mined annotations |
| `/{id}/fullTextXML` | GET | Full-text XML |
| `/{id}/supplementaryFiles` | GET | Supplementary materials |

#### Rate Limits

- Standard: ~10 requests/second
- No API key required
- Contact epmc-webservices@ebi.ac.uk for higher limits

#### Query Parameters

```
?query=SEARCH_TERM         # Search query (max 1500 chars)
?format=json|xml|dc        # Response format
?resultType=idlist|lite|core  # Detail level
?cursorMark=*              # Pagination cursor
?pageSize=25               # Results per page (max 1000)
?sort=relevance|date|cited # Sort order
```

#### Example API Calls

```bash
# Search for genetics articles with JSON response
https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=SNP%20genetics&format=json&resultType=core

# Get full-text XML
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1234567/fullTextXML

# Get text-mined terms (genes, diseases, chemicals)
https://www.ebi.ac.uk/europepmc/webservices/rest/MED/12345678/textMinedTerms?format=json

# Get citations
https://www.ebi.ac.uk/europepmc/webservices/rest/MED/12345678/citations?format=json
```

### Annotations API

**URL**: `https://europepmc.org/AnnotationsApi`

Provides programmatic access to text-mined annotations:
- Gene/protein names
- Organisms
- Diseases
- Gene Ontology terms
- Chemicals
- Accession numbers
- GeneRIF (Gene Reference into Function)
- Protein-protein interactions

### Bulk Download

**FTP Base**: `https://europepmc.org/ftp/`

| Dataset | Format | Update Frequency |
|---------|--------|------------------|
| Open Access Subset | XML, PDF | Weekly |
| Author Manuscripts | XML, Text | Weekly |
| Preprints | XML | Weekly |
| All Metadata | XML | Weekly |
| ID Mappings | CSV | Monthly |

### Data Formats

| Format | Description |
|--------|-------------|
| XML | Full metadata and content |
| JSON | Structured API responses |
| Dublin Core (DC) | RDF/XML for linked data |

### Licensing

- API: Free, no registration
- Open Access content: CC licenses
- Non-OA: Publisher copyright applies
- Must accept EBI Privacy Notice

### Size Estimates

- Full database metadata: ~50 GB compressed
- Full-text XML: ~200-300 GB compressed
- PDFs: ~500 GB+

### Full-Text Availability

| Content Type | Count | Full-Text |
|--------------|-------|-----------|
| Abstracts | 33M+ | No |
| Full-Text | 10.2M | Yes |
| Open Access | 6.5M | Yes |
| Text-minable | 2.5M+ | Yes |

---

## 3. Semantic Scholar API

### Overview

**Description**: AI-powered academic search engine providing access to scientific publications with advanced features like TLDRs, citation classification, and semantic search.

**Maintainer**: Allen Institute for AI (AI2)

**Website**: https://www.semanticscholar.org/

**API Documentation**: https://api.semanticscholar.org/api-docs/

### Content Statistics

| Metric | Value |
|--------|-------|
| Total Papers | 214+ million |
| All Disciplines | Yes |
| Citation Network | Billions of citations |
| SPECTER2 Embeddings | Available via API |

### API Access

**Base URL**: `https://api.semanticscholar.org/`

#### API Services

| Service | Base Path | Description |
|---------|-----------|-------------|
| Academic Graph | `/graph/v1/` | Papers, authors, citations |
| Recommendations | `/recommendations/v1/` | Similar papers |
| Datasets | `/datasets/v1/` | Bulk data downloads |

#### Rate Limits

| Tier | Rate | Requirements |
|------|------|--------------|
| Unauthenticated | 1000 req/sec (shared) | None |
| API Key (intro) | 1 req/sec | Free registration |
| API Key (standard) | Higher limits | Upon request |

#### Key Endpoints

```
GET /graph/v1/paper/{paper_id}           # Paper details
GET /graph/v1/paper/search               # Search papers
GET /graph/v1/paper/{paper_id}/citations # Paper citations
GET /graph/v1/paper/{paper_id}/references # Paper references
GET /graph/v1/author/{author_id}         # Author details
GET /graph/v1/author/search              # Search authors
```

#### Example API Calls

```bash
# Search for papers
https://api.semanticscholar.org/graph/v1/paper/search?query=CRISPR+genetics&fields=title,authors,abstract,year,citationCount

# Get paper by DOI
https://api.semanticscholar.org/graph/v1/paper/DOI:10.1038/nature12373?fields=title,abstract,authors,citations

# Get paper recommendations
https://api.semanticscholar.org/recommendations/v1/papers/forpaper/CorpusId:12345678
```

#### Available Fields

- `paperId`, `externalIds` (DOI, ArXiv, PubMed, etc.)
- `title`, `abstract`, `venue`, `year`
- `authors`, `fieldsOfStudy`
- `citationCount`, `referenceCount`
- `isOpenAccess`, `openAccessPdf`
- `tldr` (AI-generated summary)
- `embedding` (SPECTER2 vector)

### Bulk Datasets

**URL**: https://api.semanticscholar.org/api-docs/datasets

Available for download:
- Papers dataset
- Authors dataset
- Citations dataset
- Abstracts dataset

### Data Formats

- **API**: JSON
- **Bulk**: Compressed JSON lines

### Licensing

- API: Free with attribution
- Academic use: Permitted
- Commercial use: Contact AI2
- No full-text hosting (links only)

### Size Estimates

- API responses: JSON
- Bulk datasets: Multiple TB

### Full-Text Availability

**Links only** - Semantic Scholar provides `openAccessPdf.url` field pointing to external sources. Does not host full-text content directly.

---

## 4. OpenAlex

### Overview

**Description**: Fully open catalog of the global research system, indexing scholarly works, authors, sources, institutions, topics, and funders. Created as successor to Microsoft Academic Graph.

**Maintainer**: OurResearch (nonprofit)

**Website**: https://openalex.org/

**Documentation**: https://docs.openalex.org/

### Content Statistics

| Metric | Value |
|--------|-------|
| Total Works | 271.3 million |
| Expansion Pack (xpac) | +192 million additional |
| Authors | 90+ million |
| Sources | 250,000+ |
| Institutions | 100,000+ |
| Topics | 4,500+ |
| Updated | Daily |

### API Access

**Base URL**: `https://api.openalex.org/`

**Authentication**: API key recommended (free)

#### Entity Types

| Entity | Endpoint | Description |
|--------|----------|-------------|
| Works | `/works` | Scholarly papers |
| Authors | `/authors` | Researchers |
| Sources | `/sources` | Journals, repositories |
| Institutions | `/institutions` | Universities, labs |
| Topics | `/topics` | Subject areas |
| Publishers | `/publishers` | Publishing organizations |
| Funders | `/funders` | Funding bodies |

#### Rate Limits

| Tier | Limit | Requirements |
|------|-------|--------------|
| Standard | 100,000 calls/day | Email in polite pool |
| Premium | Higher limits | Paid subscription |

#### Query Parameters

```
?filter=field:value              # Filter results
?search=query                    # Full-text search
?sort=field:asc|desc            # Sort results
?per_page=25                     # Results per page (max 200)
?cursor=*                        # Pagination cursor
?select=field1,field2           # Select specific fields
?sample=100                      # Random sample
```

#### Example API Calls

```bash
# Search for genetics works
https://api.openalex.org/works?search=CRISPR%20gene%20editing&per_page=50

# Filter by open access
https://api.openalex.org/works?filter=open_access.is_oa:true,publication_year:2024

# Get work by DOI
https://api.openalex.org/works/doi:10.1038/nature12373

# Get author with works
https://api.openalex.org/authors/A1234567890?select=id,display_name,works_count
```

### Bulk Download

**AWS S3**: `s3://openalex` (Registry of Open Data on AWS)

**Format**: JSON lines, partitioned by entity type

**Update**: Daily snapshots

### Data Formats

- **API**: JSON
- **Bulk**: Compressed JSON lines (.gz)

### Licensing

- **Data**: CC0 (public domain)
- **API**: Free, attribution appreciated
- **Commercial use**: Allowed

### Size Estimates

- Full dataset: ~300 GB compressed
- Works only: ~200 GB compressed

### Full-Text Availability

**Links only** - OpenAlex tracks open access status via Unpaywall integration. Provides:
- `open_access.is_oa`: Boolean
- `open_access.oa_status`: gold, green, hybrid, bronze, closed
- `open_access.oa_url`: Link to OA version

---

## 5. CORE (COnnecting REpositories)

### Overview

**Description**: World's largest aggregator of open access research papers, collecting from institutional repositories, subject repositories, and open access journals globally.

**Maintainer**: Knowledge Media Institute, The Open University (UK)

**Website**: https://core.ac.uk/

**API Documentation**: https://core.ac.uk/services/api

### Content Statistics

| Metric | Value |
|--------|-------|
| Metadata Records | 391+ million |
| Full-Text Links | 323 million (free to read) |
| Full-Text Hosted | 38-46 million |
| Data Providers | 14,000+ |
| Countries | 150+ |

### API Access

**Base URL**: `https://api.core.ac.uk/v3/`

**Documentation**: https://api.core.ac.uk/docs/v3

#### Rate Limits

| Tier | Rate | Requirements |
|------|------|--------------|
| Unauthenticated | 1 batch or 5 single requests per 10 sec | None |
| Registered | Faster rates | Free API key |
| Enterprise | Custom rates | Contact CORE |

#### Key Endpoints

| Endpoint | Description |
|----------|-------------|
| `/search/works` | Search papers |
| `/works/{id}` | Get paper by ID |
| `/search/outputs` | Search research outputs |
| `/data-providers` | List data providers |
| `/discover` | Discovery service |

#### Example API Calls

```bash
# Search for papers
https://api.core.ac.uk/v3/search/works?q=machine%20learning&limit=10&api_key=YOUR_KEY

# Get paper by CORE ID
https://api.core.ac.uk/v3/works/12345678?api_key=YOUR_KEY

# Search with filters
https://api.core.ac.uk/v3/search/works?q=genetics&filter=year:2024&api_key=YOUR_KEY
```

### Services

| Service | Description |
|---------|-------------|
| CORE Discovery | Find OA copies of papers |
| CORE Recommender | Related paper suggestions |
| CORE FastSync | Real-time content sync |
| Repository Dashboard | Analytics for repositories |

### Bulk Download

Contact CORE for bulk data access. Available through:
- CORE Dataset (research purposes)
- CORE FastSync (commercial partners)

### Data Formats

- **API**: JSON
- **Bulk**: JSON, XML

### Licensing

- **Research use**: Free
- **Commercial use**: Paid agreements
- **Content**: Depends on source repository licenses

### Size Estimates

- Full metadata: ~100+ GB
- Full-text corpus: Multiple TB

### Full-Text Availability

| Type | Count |
|------|-------|
| Metadata only | ~350M |
| Full-text links | 323M |
| Full-text hosted | 38-46M |

CORE provides direct access to full-text PDFs and text for millions of papers, unlike aggregators that only provide links.

---

## 6. Unpaywall

### Overview

**Description**: Open database for finding legal free versions of scholarly articles. Looks up DOIs and returns links to open access copies from repositories and publisher sites.

**Maintainer**: OurResearch (nonprofit)

**Website**: https://unpaywall.org/

### Content Statistics

| Metric | Value |
|--------|-------|
| DOIs Indexed | 95+ million (Crossref) |
| Open Access Found | 24-27 million articles |
| Publishers | 50,000+ |
| Repositories | 50,000+ |

### API Access

**Base URL**: `https://api.unpaywall.org/v2/`

**Authentication**: Email address required as parameter

#### Rate Limits

| Usage | Limit |
|-------|-------|
| Standard | 100,000 calls/day |
| Bulk | Use database snapshot |

#### Endpoints

| Endpoint | Description |
|----------|-------------|
| `/{doi}` | Get OA status for DOI |

#### API Parameters

```
?email=your@email.com    # Required - your email address
```

#### Response Fields

| Field | Description |
|-------|-------------|
| `is_oa` | Boolean - is open access? |
| `oa_status` | gold, green, hybrid, bronze, closed |
| `best_oa_location` | Best available OA copy |
| `oa_locations` | All available OA copies |

#### OA Location Object

```json
{
  "url": "https://...",
  "url_for_pdf": "https://...",
  "url_for_landing_page": "https://...",
  "evidence": "oa repository",
  "host_type": "repository",
  "is_best": true,
  "license": "cc-by",
  "version": "publishedVersion"
}
```

#### Example API Calls

```bash
# Check OA status for a DOI
https://api.unpaywall.org/v2/10.1038/nature12373?email=your@email.com

# Response includes best_oa_location with PDF URL
```

### Bulk Download

**Snapshot**: Available for download (updated regularly)

**Format**: JSON lines

**Access**: https://unpaywall.org/products/snapshot

### Data Sources

Unpaywall aggregates from:
- Crossref (DOI registration)
- DOAJ (Open Access journals)
- OAI-PMH repositories
- PubMed Central
- Institutional repositories
- Publisher websites

### Licensing

- **API**: Free for all uses
- **Data**: CC0 (public domain)
- **Attribution**: Appreciated but not required

### Size Estimates

- Database snapshot: ~20 GB compressed
- API responses: Small JSON objects

### Full-Text Availability

**Links only** - Unpaywall does not host content. It provides:
- `url_for_pdf`: Direct PDF link
- `url_for_landing_page`: Article page
- Multiple locations ranked by quality

---

## 7. Preprint Servers (bioRxiv / medRxiv)

### Overview

**Description**: Preprint repositories for biology (bioRxiv) and health sciences (medRxiv). Operated by Cold Spring Harbor Laboratory (transferred to openRxiv nonprofit in 2025).

**Maintainer**: openRxiv (formerly Cold Spring Harbor Laboratory)

**Websites**:
- https://www.biorxiv.org/
- https://www.medrxiv.org/

**API Documentation**: https://api.biorxiv.org/

### Content Statistics

| Server | Preprints | As Of |
|--------|-----------|-------|
| bioRxiv | ~180,000+ | 2022 baseline, growing |
| medRxiv | 61,000+ | December 2024 |
| Combined | 250,000+ | Estimated 2025 |
| COVID-19 | 31,198 | Combined subset |

**Publication Rate**:
- bioRxiv: ~36,000-40,000 per year
- ~67% of bioRxiv papers later published in peer-reviewed journals

### API Access

**Base URLs**:
- `https://api.biorxiv.org/`
- `https://api.medrxiv.org/`

**Authentication**: None required

**Rate Limits**: No specified limits (1-second courtesy delay recommended)

#### Endpoints

| Endpoint | Format | Description |
|----------|--------|-------------|
| `/details/[server]/[interval]/[cursor]/[format]` | JSON, XML, HTML | Article metadata |
| `/pubs/[server]/[interval]/[cursor]` | JSON | Published versions |
| `/pub/[interval]/[cursor]/[format]` | JSON, XML, HTML | bioRxiv only |
| `/publisher/[prefix]/[interval]/[cursor]` | JSON | Publisher-specific |
| `/funder/[server]/[interval]/[funder]/[cursor]/[format]` | JSON | Funder filtering |
| `/sum/[interval]/[format]` | JSON, XML, CSV | Statistics |
| `/usage/[interval]/[server]/[format]` | JSON, XML, CSV | Download stats |

#### Interval Formats

| Format | Example | Description |
|--------|---------|-------------|
| Date range | `2024-01-01/2024-12-31` | YYYY-MM-DD/YYYY-MM-DD |
| Recent N | `100` | N most recent articles |
| Recent days | `30d` | Last N days |

#### Example API Calls

```bash
# Get recent bioRxiv preprints (last 30 days)
https://api.biorxiv.org/details/biorxiv/30d/0/json

# Get medRxiv preprints by date range
https://api.medrxiv.org/details/medrxiv/2024-01-01/2024-06-30/0/json

# Get bioRxiv by category
https://api.biorxiv.org/details/biorxiv/2024-01-01/2024-12-31/0/json?category=genetics

# Get statistics
https://api.biorxiv.org/sum/2024-01-01/2024-12-31/json

# Get usage data
https://api.biorxiv.org/usage/2024-01-01/2024-12-31/biorxiv/json
```

#### Response Pagination

- Results paginated at 100 papers per call
- Use cursor parameter to iterate
- Response includes `cursor` field for next page

### Data Fields

| Field | Description |
|-------|-------------|
| `doi` | Preprint DOI |
| `title` | Article title |
| `authors` | Author list |
| `author_corresponding` | Corresponding author |
| `author_corresponding_institution` | Affiliation |
| `date` | Posting date |
| `version` | Version number |
| `type` | Article type |
| `license` | License type |
| `category` | Subject category |
| `jatsxml` | Link to JATS XML |
| `abstract` | Abstract text |
| `published` | Published DOI (if applicable) |

### Bulk Download

- API supports bulk retrieval via date ranges
- Full repository download: ~1 hour for bioRxiv
- OAI-PMH available for metadata harvesting

### Data Formats

| Format | Description |
|--------|-------------|
| JSON | Default API format |
| XML | OAI-PMH XML |
| HTML | Human-readable |
| CSV | Statistics endpoints |
| JATS XML | Full article content |

### Licensing

| License | Description |
|---------|-------------|
| CC BY | Most common |
| CC BY-NC | Non-commercial |
| CC BY-ND | No derivatives |
| CC0 | Public domain |

All preprints are open access by nature.

### Size Estimates

- bioRxiv: ~50-100 GB full text
- medRxiv: ~20-30 GB full text
- API responses: JSON

### Full-Text Availability

**100%** - All preprints include full text. Available via:
- PDF download
- JATS XML (via `jatsxml` field)
- HTML on website

---

## Comparison Matrix

### Coverage by Discipline

| Source | Life Sciences | Medicine | All Sciences | Preprints |
|--------|---------------|----------|--------------|-----------|
| PMC | Strong | Strong | Limited | Yes |
| Europe PMC | Strong | Strong | Limited | Yes |
| Semantic Scholar | Moderate | Moderate | Strong | Yes |
| OpenAlex | Moderate | Moderate | Strong | Yes |
| CORE | Moderate | Moderate | Strong | Yes |
| Unpaywall | Moderate | Moderate | Strong | Limited |
| bioRxiv/medRxiv | Strong | Strong | No | Yes |

### API Comparison

| Feature | PMC | Europe PMC | Semantic Scholar | OpenAlex | CORE | Unpaywall | bioRxiv |
|---------|-----|------------|------------------|----------|------|-----------|---------|
| Auth Required | No | No | No | No | Optional | Email | No |
| JSON Support | Partial | Full | Full | Full | Full | Full | Full |
| Rate Limit | 3-10/sec | ~10/sec | 1-100/sec | 100K/day | Tiered | 100K/day | None |
| Full-Text | Yes | Yes | Links | Links | Yes | Links | Yes |
| Annotations | BioC | Yes | TLDRs | Topics | No | No | No |
| Bulk Download | Yes | Yes | Yes | Yes | Contact | Yes | Yes |

### Recommended Use Cases

| Use Case | Primary | Secondary | Notes |
|----------|---------|-----------|-------|
| RAG for Genetics | Europe PMC | PMC | Text-mined annotations |
| Citation Network | Semantic Scholar | OpenAlex | Citation classification |
| OA Detection | Unpaywall | OpenAlex | DOI lookup |
| Preprint Tracking | bioRxiv/medRxiv | Europe PMC | Direct full-text |
| Bulk Analysis | OpenAlex | CORE | Large-scale datasets |
| Knowledge Base | PMC OA Subset | Europe PMC | JATS XML parsing |

---

## Implementation Recommendations

### Phase 1: Real-Time APIs

1. **Europe PMC** for search with annotations
2. **Unpaywall** for OA detection by DOI
3. **Semantic Scholar** for recommendations and embeddings

### Phase 2: Bulk Data

1. **PMC Open Access Subset** for full-text corpus
2. **OpenAlex** for comprehensive metadata
3. **bioRxiv/medRxiv** for preprint coverage

### Phase 3: Enrichment

1. **Europe PMC Annotations** for text mining
2. **Semantic Scholar embeddings** for semantic search
3. **OpenAlex topics** for classification

---

## References

- [PMC Open Access Subset](https://pmc.ncbi.nlm.nih.gov/tools/openftlist/)
- [PMC OA Web Service](https://pmc.ncbi.nlm.nih.gov/tools/oa-service/)
- [PMC BioC API](https://www.ncbi.nlm.nih.gov/research/bionlp/APIs/BioC-PMC/)
- [Europe PMC REST API](https://europepmc.org/RestfulWebService)
- [Europe PMC Annotations API](https://europepmc.org/AnnotationsApi)
- [Semantic Scholar API](https://api.semanticscholar.org/api-docs/)
- [OpenAlex Documentation](https://docs.openalex.org/)
- [OpenAlex on AWS](https://registry.opendata.aws/openalex/)
- [CORE API](https://core.ac.uk/services/api)
- [CORE Scientific Data Paper](https://www.nature.com/articles/s41597-023-02208-w)
- [Unpaywall API](https://unpaywall.org/products/api)
- [bioRxiv API](https://api.biorxiv.org/)
- [bioRxiv Wikipedia](https://en.wikipedia.org/wiki/BioRxiv)
- [medRxiv Wikipedia](https://en.wikipedia.org/wiki/MedRxiv)
