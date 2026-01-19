# Alternative/Shadow Library Sources for Research Papers

> **Research Date**: January 2026
> **Purpose**: Evaluate data sources for genetics/health knowledge base platform
> **Constraint**: Platform can only use publicly available data - no commercial licenses

---

## CRITICAL LEGAL WARNING FOR STARTUPS

**Before evaluating these sources, understand the legal landscape:**

As of 2025-2026, shadow libraries are at the center of major copyright litigation. Multiple AI companies (Anthropic, OpenAI, Meta, Apple, Nvidia) face lawsuits specifically for downloading content from shadow libraries like LibGen and Anna's Archive.

**Key Legal Precedents:**
- **Bartz v. Anthropic (2025)**: Judge ruled that downloading books from LibGen constituted willful copyright infringement. Settlement: $1.5 billion.
- **Meta Litigation (2025)**: Unsealed emails revealed Meta downloaded 81+ TB from Anna's Archive torrents.
- **OpenAI (2024-2025)**: Court allowed plaintiffs to pursue shadow library downloading as a separate infringement claim.

**Bottom Line**: For a bootstrapped startup, using shadow library content for commercial purposes carries **extreme legal risk**. Even downloading for "research" creates liability. These sources are documented here for completeness but are **NOT RECOMMENDED** for commercial use.

---

## 1. Sci-Hub

### Current Status (2025-2026)
- **Operational**: Yes, but with limitations
- **Collection Size**: ~88-90 million documents (as of 2022 data)
- **Coverage**: 90.4% of paywalled articles; 92.3% of closed-access content
- **Major Limitation**: Fresh content ingestion halted since late 2020/2021 due to litigation

### Recent Developments
- **August 2025**: Blocked in India following lawsuit by Elsevier, Wiley, and American Chemical Society
- **2025**: Launch of "Sci-Net" - invite-only social network for paper requests (by Alexandra Elbakyan)
- Domain changes continue due to legal pressure

### Access Methods
| Method | Details |
|--------|---------|
| Direct Mirrors | sci-hub.se, sci-hub.st, sci-hub.ru (change frequently) |
| DOI Lookup | Enter DOI directly on homepage |
| TOR | Available via .onion addresses |
| VPN | Recommended in blocked regions |

### Bulk Download Feasibility
- **Torrents**: Full database available (~100 TB total)
- **Source**: Provided via Library Genesis project (not Sci-Hub directly)
- **SQL Metadata**: Not publicly available from Sci-Hub; LibGen provides "scimag" database
- **Latest SQL Dump**: 2020 (outdated, ~82 million records)

### Data Format
- PDFs organized by DOI
- Metadata via LibGen scimag database (MySQL/MariaDB dumps)

### Legal Risk: **EXTREME**
- Subject to ongoing publisher lawsuits globally
- Downloading constitutes copyright infringement in most jurisdictions
- Commercial use would invite immediate legal action

---

## 2. Anna's Archive (annas-archive.org/.li/.pm)

### Current Status (2025-2026)
- **Operational**: Yes, with domain changes
- **Collection**: 61.6 million books + 95.7 million papers
- **Total Size**: ~1.1 petabytes across all torrents
- **Daily Downloads**: ~650,000 (March 2025)

### Recent Developments
- **January 2026**: Primary .org domain suspended; .pm and .in mirrors created
- **December 2025**: Scraped 300TB from Spotify (256M tracks)
- **October 2025**: Blocked in Germany
- **July 2025**: Blocked in Belgium (500,000 euro fines for noncompliance)
- Domains also blocked in Italy, Netherlands, UK

### Access Methods
| Method | Details |
|--------|---------|
| Direct | annas-archive.li, annas-archive.pm, annas-archive.se |
| Search | By title, author, DOI, ISBN, or MD5 |
| Torrents | Full collection available for bulk download |
| SFTP | High-speed access for LLM training ($100K+ or data contribution) |

### Torrent Structure
Anna's Archive provides the "ultimate unified list" of releases from:
- Sci-Hub (scientific papers)
- Library Genesis (books + scimag)
- Z-Library collections
- Internet Archive digitized materials
- DuXiu, MagzDB, Nexus/STC, HathiTrust

**Seeding Status**: 69% of 1.1PB copied in 4+ locations; only 8% in 10+ locations

### Metadata & API
- **Search**: DOI, ISBN, title, author, MD5
- **Database**: Can generate/download as ElasticSearch and MariaDB
- **Third-party API**: Available on RapidAPI
- **MCP Server**: Available for programmatic access

### Bulk Download
- GitHub tool: `cparthiv/annas-torrents` - download by specifying TB amount
- Meta reportedly downloaded 81.7 TB through their torrents
- 30+ companies (mostly Chinese) pay for SFTP access

### Legal Risk: **EXTREME**
- Multiple blocking orders across EU
- Founders/operators anonymous
- Commercial use = guaranteed litigation target

---

## 3. Library Genesis (LibGen)

### Current Status (2025-2026)
- **Operational**: Yes, but under heavy legal pressure
- **Collection**: 2.4M+ non-fiction books, 80M+ scientific articles
- **Status**: Majority of domains seized in December 2024 (Pearson Education lawsuit)

### Recent Developments
- **December 2024**: Domain seizures by publishers led by Pearson Education
- **December 2024**: Germany-wide blocking order by CUII
- **2025**: libgen.fun domain seized

### Working Mirrors (as of January 2026)
| Mirror | Notes |
|--------|-------|
| libgen.rs | Fast, frequently updated |
| libgen.is | Full mirror with Sci-Hub integration |
| libgen.st | Often unblocked globally |
| libgen.li | Stable interface |
| libgen.gs | Backup mirror |

**Direct IPs**: 93.174.95.27, 185.39.10.101

### Databases
LibGen maintains separate databases:
- **libgen** (main): Non-fiction books
- **fiction**: Fiction books
- **scimag**: Scientific articles (Sci-Hub content)
- **comics**: Comic books
- **standards**: Technical standards

### Metadata & Bulk Download
- **SQL Dumps**: Available on Internet Archive
  - Latest: 2023-05-13 snapshot
  - Format: MySQL/MariaDB SQL dumps
- **Figshare**: 2017 scimag dump available as TSV
- **GitHub**: `greenelab/scihub` - processing scripts

### Data Format
- SQL database with DOI, title, author, year, journal, MD5
- Files accessible via MD5 hash

### Legal Risk: **EXTREME**
- Active Pearson Education lawsuit
- Elsevier lawsuit (2015) ongoing
- Domains actively seized by law enforcement

---

## 4. Z-Library

### Current Status (2025-2026)
- **Operational**: Yes, after 2022 shutdown
- **Collection**: 13.35M+ books, 84.8M+ articles (Feb 2023 data)
- **Official Domain**: z-lib.id (confirmed via official channels)

### Legal History
- **November 2022**: FBI seized domains, arrested Russian founders
- **July 2024**: Founders escaped house arrest in Argentina; Interpol warrant issued
- **2025**: Returned with new domain system

### Access Methods
| Method | Details |
|--------|---------|
| Official Domain | z-lib.id only |
| TOR | zlibrary24tuxziyiyfr7zd46ytefdqbqd2axkmxm4o5374ptpc52fad.onion |
| Email Login | Send blank email to official inbox for personal link |
| Desktop App | Version 3.0.0 available |
| Mobile APK | Available from official mirrors |

### Download Limits
- Registered users: 50 downloads/day
- Anonymous: 5 downloads/day

### Bulk Download Feasibility
- **Not directly supported**
- No public API or torrent system
- Content accessible via Anna's Archive torrents

### Legal Risk: **EXTREME**
- Criminal charges against founders (copyright infringement, wire fraud, money laundering)
- US DOJ actively pursuing case
- Using Z-Library = potential criminal liability in US

---

## 5. Internet Archive Scholar

### Current Status (2025-2026)
- **Operational**: Yes
- **Collection**: 35+ million research articles with full text
- **Legal Status**: LEGITIMATE (mostly)

### Recent Developments
- **July 2025**: Designated as Federal Depository Library by US Senate
- **September 2025**: New European headquarters opened
- **2025**: Total collection exceeds 47 million texts
- **Ongoing**: Legal battles over Controlled Digital Lending

### Content Sources
1. **Web Archive**: Open access papers from Wayback Machine
2. **Digitized Materials**: Paper and microform collections
3. **Partner Collections**: Archive.org collaborations

### Access Methods
| Method | Details |
|--------|---------|
| Web Search | scholar.archive.org |
| API | Open REST API via Fatcat |
| Bulk Downloads | Metadata dumps available |

### API & Bulk Access
- **Fatcat**: Open catalog with read/write API, CLI tool
- **Metadata**: Bulk dumps available at fatcat.wiki
- **IAS3 API**: Upload/download/search functionality
- **R Package**: `internetarchive` for programmatic access

### Data Coverage
- Focus on "at-risk" open access content
- Archived faculty manuscripts
- Vanished OA publisher content
- Digitized microfilm of older publications

### Legal Risk: **LOW-MODERATE**
- Legitimate non-profit organization
- Focus on Open Access and public domain
- Controlled Digital Lending faces legal challenges
- Metadata is generally safe to use

---

## LEGAL ALTERNATIVES (RECOMMENDED)

### Crossref
- **What**: DOI registration agency with metadata API
- **Collection**: 156+ million metadata records (2024), 197GB total
- **API**: Free, no authentication required
- **License**: CC0 (public domain) for most metadata
- **Bulk Download**: Annual data file available
- **Legal Risk**: **NONE** - fully legal metadata source
- **URL**: api.crossref.org

### OpenAlex
- **What**: Open bibliographic catalog (Microsoft Academic successor)
- **Collection**: 240+ million works, growing 50K/day
- **API**: Free, 100K requests/day, no auth required
- **License**: CC0
- **Sources**: Crossref, PubMed, arXiv, institutional repositories
- **Legal Risk**: **NONE** - fully legal
- **Funding**: $7.5M grant from Arcadia (2024)
- **URL**: openalex.org

### Unpaywall
- **What**: Database of legal OA article locations
- **Collection**: 20+ million free, legal PDFs indexed
- **Success Rate**: ~52% of papers have free legal version
- **API**: REST API available
- **License**: Legal - indexes OA locations only
- **Legal Risk**: **NONE**
- **URL**: unpaywall.org

### PubMed Central (PMC)
- **What**: NIH's free full-text archive
- **Collection**: 8+ million articles
- **API**: E-utilities API
- **Bulk**: FTP download available
- **Legal Risk**: **NONE** - government-mandated open access
- **URL**: ncbi.nlm.nih.gov/pmc

### arXiv
- **What**: Preprint server for physics, math, CS, biology
- **Collection**: 2+ million preprints
- **API**: Available
- **Bulk**: Bulk download available with agreement
- **Legal Risk**: **NONE** - author-submitted preprints
- **URL**: arxiv.org

### bioRxiv/medRxiv
- **What**: Preprint servers for biology/medicine
- **API**: Available
- **Legal Risk**: **NONE**

---

## COMPARISON MATRIX

| Source | Papers | Bulk DL | API | Metadata | Legal Risk | Commercial Use |
|--------|--------|---------|-----|----------|------------|----------------|
| Sci-Hub | 88M+ | Torrent | No | Via LibGen | EXTREME | NO |
| Anna's Archive | 95M+ | Torrent | Unofficial | Yes | EXTREME | NO |
| LibGen | 80M+ | SQL Dump | No | Yes | EXTREME | NO |
| Z-Library | 84M+ | No | No | Limited | EXTREME | NO |
| IA Scholar | 35M+ | Yes | Yes | Yes | LOW | MAYBE |
| **Crossref** | 156M meta | Yes | Yes | Yes | NONE | YES |
| **OpenAlex** | 240M+ | Yes | Yes | Yes | NONE | YES |
| **Unpaywall** | 20M OA | Yes | Yes | Yes | NONE | YES |
| **PMC** | 8M+ | Yes | Yes | Yes | NONE | YES |

---

## RECOMMENDATIONS FOR BOOTSTRAPPED STARTUP

### DO NOT USE (Legal Exposure Too High):
1. Sci-Hub - Active litigation, criminal implications
2. Anna's Archive - Multiple lawsuits, extreme risk
3. Library Genesis - Domains being seized
4. Z-Library - Criminal charges against operators

### SAFE TO USE:
1. **OpenAlex** - Best overall: large, free, legal, good API
2. **Crossref** - Metadata only, but comprehensive
3. **Unpaywall** - Index of legal OA locations
4. **PubMed Central** - Strong for biomedical content
5. **arXiv/bioRxiv** - Preprints in relevant fields

### HYBRID STRATEGY:
1. Use **OpenAlex** for comprehensive metadata (240M works)
2. Use **Unpaywall** to find legal full-text locations
3. Use **PMC/arXiv/bioRxiv** for direct full-text access
4. Index and link to legal sources rather than hosting content
5. Build value-add layer (analysis, connections, summaries)

### WHAT YOU CAN BUILD LEGALLY:
- Metadata search and discovery platform
- Citation network analysis
- Research trend visualization
- Paper recommendation system
- Links to legal full-text sources
- Summaries and analysis of publicly available papers

### WHAT YOU CANNOT DO LEGALLY:
- Host copyrighted PDFs
- Bulk download from shadow libraries
- Redistribute publisher content
- Scrape and store paywalled content

---

## APPENDIX: Technical Details

### LibGen Scimag SQL Schema (Key Fields)
```sql
- ID (integer)
- DOI (varchar)
- Title (text)
- Author (text)
- Year (integer)
- Journal (varchar)
- Volume, Issue, Pages
- MD5 (file hash)
- Filesize
- TimeAdded
```

### OpenAlex API Example
```bash
# Get work by DOI
curl "https://api.openalex.org/works/doi:10.1038/nature12373"

# Search works
curl "https://api.openalex.org/works?search=genetics&per-page=25"

# Filter by topic
curl "https://api.openalex.org/works?filter=topics.id:T12345"
```

### Crossref API Example
```bash
# Get work by DOI
curl "https://api.crossref.org/works/10.1038/nature12373"

# Search
curl "https://api.crossref.org/works?query=CRISPR+genetics"
```

### Anna's Archive Torrent Categories
- aa_derived_mirror_metadata (metadata databases)
- scihub (scientific papers from Sci-Hub)
- libgen_rs_nf (LibGen non-fiction)
- libgen_rs_fic (LibGen fiction)
- zlib (Z-Library mirrors)
- ia (Internet Archive collections)

---

## SOURCES

### Primary Research
- [Sci-Hub Wikipedia](https://en.wikipedia.org/wiki/Sci-Hub)
- [Anna's Archive Wikipedia](https://en.wikipedia.org/wiki/Anna's_Archive)
- [Library Genesis Wikipedia](https://en.wikipedia.org/wiki/Library_Genesis)
- [Z-Library Wikipedia](https://en.wikipedia.org/wiki/Z-Library)
- [Internet Archive Scholar](https://scholar.archive.org/)

### Legal Coverage
- [Shadow Libraries AI Lawsuits](https://chatgptiseatingtheworld.com/2025/10/26/status-of-all-56-copyright-lawsuits-v-ai-oct-26-2025-apple-hit-with-3-lawsuits/)
- [Anthropic Settlement](https://roninlegalconsulting.com/the-anthropic-settlement-that-wasnt/)
- [Meta Anna's Archive Downloads](https://torrentfreak.com/meta-torrented-over-81-tb-of-data-through-annas-archive-despite-few-seeders-250206/)

### Technical Resources
- [OpenAlex Documentation](https://docs.openalex.org/)
- [Crossref API](https://www.crossref.org/documentation/retrieve-metadata/)
- [Unpaywall](https://unpaywall.org/)
- [LibGen Metadata Snapshots](https://archive.org/details/libgen-meta-20230513)

---

*Last Updated: January 2026*
*This document is for research purposes only. Consult legal counsel before using any data source commercially.*
