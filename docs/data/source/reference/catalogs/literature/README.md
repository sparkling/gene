---
title: "Literature Databases"
parent: ../_index.md
category: shared
last_updated: 2026-01-22
status: draft
---

# Literature Databases

Scientific literature sources for research articles, evidence synthesis, and knowledge extraction.

## Database Catalog

| Database | Type | Tier | Coverage | Access Method | Size |
|----------|------|------|----------|---------------|------|
| **PubMed** | Biomedical Literature | 1 | 36+ million citations | E-utilities API, Entrez | MEDLINE index |
| **PubMed Central (PMC)** | Full-Text Archive | 2 | 10+ million full-text | E-utilities API, FTP | Open access |
| **OpenAlex** | Open Scholarly Graph | 2 | 250+ million works | REST API, Snapshot | Comprehensive |
| **Europe PMC** | European Archive | 2 | 40+ million abstracts | REST API | European focus |
| **Semantic Scholar** | AI-Powered Search | 2 | 200+ million papers | REST API | Citation graph |
| **CORE** | Open Access Aggregator | 3 | 200+ million articles | REST API | OA-focused |
| **arXiv** | Preprint Archive | 3 | 2+ million preprints | OAI-PMH, API | Physics, CS, Math |
| **bioRxiv** | Biology Preprints | 3 | 200,000+ preprints | API | Life sciences |
| **ClinicalTrials.gov** | Clinical Trials | 2 | 450,000+ studies | API, Download | Trial registry |

## Primary Use Cases

### Tier 1 (MVP) Focus

#### PubMed
- **Purpose**: Biomedical and life sciences literature index
- **Content**:
  - 36+ million citations
  - MEDLINE, PMC, BookShelf
  - Abstracts, MeSH terms, metadata
- **Use**:
  - Literature search
  - Evidence extraction
  - Citation mapping
  - Gene-disease-compound associations
- **Access**: https://pubmed.ncbi.nlm.nih.gov/
  - **E-utilities API**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`
  - No API key required (rate-limited)
  - API key recommended for high volume (free)
- **Update**: Real-time (multiple updates daily)

## E-utilities (Entrez Programming Utilities)

### Core E-utility Tools

| Tool | Purpose | Example |
|------|---------|---------|
| `esearch` | Search database | Find PMIDs for "APOE Alzheimer" |
| `efetch` | Retrieve records | Get abstract for PMID 12345678 |
| `esummary` | Get document summaries | Get metadata for PMIDs |
| `elink` | Find related records | Link genes to articles |
| `einfo` | Get database info | List available databases |
| `epost` | Upload ID list | Store search results |
| `egquery` | Global query | Search all NCBI databases |

## Access Methods

### PubMed E-utilities

#### 1. Search (esearch)
```bash
# Basic search
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=APOE+Alzheimer&retmode=json"

# Search with date range
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=APOE+Alzheimer&mindate=2020&maxdate=2024&retmode=json"

# Complex query
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(APOE[Gene]+AND+Alzheimer[Disease/Abnormality])+AND+2023[PDAT]&retmax=100&retmode=json"
```

#### 2. Fetch (efetch)
```bash
# Get abstract (XML format)
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12345678&retmode=xml"

# Get abstract (text format)
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12345678&rettype=abstract&retmode=text"

# Get multiple abstracts
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12345678,23456789&retmode=xml"
```

#### 3. Link (elink)
```bash
# Link gene to articles
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=pubmed&id=348" # APOE gene

# Find related articles
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pubmed&id=12345678&cmd=neighbor"
```

#### 4. Summary (esummary)
```bash
# Get document summaries
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=12345678,23456789&retmode=json"
```

### API Key Usage
```bash
# Register for free API key: https://www.ncbi.nlm.nih.gov/account/

# Use API key (increases rate limit from 3/sec to 10/sec)
API_KEY="your_api_key_here"
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=APOE&api_key=$API_KEY"
```

## Tier 2 Databases

### PubMed Central (PMC)
- **Purpose**: Full-text open access archive
- **Content**: 10+ million full-text articles
- **Use**: Extract detailed methods, results, figures
- **Access**: Same E-utilities as PubMed
  - `db=pmc` instead of `db=pubmed`
- **FTP**: `ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/`

### OpenAlex
- **Purpose**: Open scholarly graph (successor to Microsoft Academic)
- **Content**:
  - 250+ million works
  - Authors, institutions, concepts
  - Citations, topics
- **Use**: Comprehensive literature analysis, citation networks
- **Access**: https://openalex.org/
  - REST API: `https://api.openalex.org/`
  - No authentication required
  - Snapshot downloads available

### Semantic Scholar
- **Purpose**: AI-powered academic search
- **Content**:
  - 200+ million papers
  - Citation context
  - Influential citations
- **Use**: Smart literature discovery, citation analysis
- **Access**: https://www.semanticscholar.org/
  - REST API: `https://api.semanticscholar.org/`
  - API key required (free)

### Europe PMC
- **Purpose**: European biomedical literature archive
- **Content**: 40+ million abstracts, full-text
- **Use**: European literature, preprints
- **Access**: https://europepmc.org/
  - REST API: `https://www.ebi.ac.uk/europepmc/webservices/rest/`

## Common Query Patterns

### Gene-Disease Literature Search
```bash
# Search for gene-disease associations
GENE="APOE"
DISEASE="Alzheimer"

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=pubmed&\
term=${GENE}[Gene]+AND+${DISEASE}[Disease/Abnormality]&\
retmax=100&\
retmode=json"
```

### Compound-Disease Literature
```bash
# Search for compound therapeutic effects
COMPOUND="Curcumin"
DISEASE="Cancer"

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=pubmed&\
term=${COMPOUND}[Title/Abstract]+AND+${DISEASE}[MeSH]&\
retmax=100&\
retmode=json"
```

### Recent High-Impact Articles
```bash
# Search for recent articles in high-impact journals
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=pubmed&\
term=APOE+AND+(Nature[Journal]+OR+Science[Journal]+OR+Cell[Journal])&\
mindate=2023&\
retmax=50&\
retmode=json"
```

### Clinical Trials for Condition
```bash
# Search ClinicalTrials.gov via PubMed
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=pubmed&\
term=Alzheimer[Disease]+AND+Clinical+Trial[Publication+Type]&\
retmax=100&\
retmode=json"
```

## Data Integration Workflow

```
Research Question
    ↓
PubMed Search (esearch)
    ↓
PMID List
    ↓
┌──────────┬──────────┐
│          │          │
efetch     elink      esummary
│          │          │
Full       Related    Metadata
Text       Articles
│          │          │
└──────────┴──────────┘
         ↓
Text Mining / NLP
         ↓
Gene-Disease-Compound Associations
```

## Text Mining Pipeline

### Standard Workflow
```
1. Query PubMed for relevant articles
   ↓
2. Fetch abstracts (or full text from PMC)
   ↓
3. Extract entities:
   - Genes (NER: gene mention detection)
   - Diseases (MeSH terms, UMLS)
   - Compounds (chemical entity recognition)
   - Relationships (dependency parsing)
   ↓
4. Store structured knowledge:
   - Gene-disease associations
   - Compound-target relationships
   - Treatment outcomes
   ↓
5. Rank by evidence strength:
   - Number of articles
   - Journal impact factor
   - Study design (RCT > observational)
```

## Rate Limiting

### PubMed E-utilities
- **Without API key**: 3 requests/second
- **With API key**: 10 requests/second
- **Recommendation**: Always use API key for production

### OpenAlex
- **Rate limit**: 100,000 requests/day (polite pool)
- **Recommendation**: Use `mailto` parameter for higher limits

### Semantic Scholar
- **Rate limit**: 100 requests/5 minutes (with API key)
- **Recommendation**: Batch requests

## Example Use Cases

### Use Case 1: Find Evidence for Gene-Disease Link
```bash
# Search for APOE-Alzheimer's literature
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=pubmed&\
term=APOE[Gene]+AND+Alzheimer+Disease[MeSH]&\
retmax=100&\
retmode=json" > pmids.json

# Fetch abstracts
cat pmids.json | jq -r '.esearchresult.idlist[]' | while read pmid; do
  curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
db=pubmed&\
id=${pmid}&\
rettype=abstract&\
retmode=text" >> abstracts.txt
done
```

### Use Case 2: Find Recent TCM Research
```bash
# Search for recent Traditional Chinese Medicine research
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=pubmed&\
term=Traditional+Chinese+Medicine[Title]+AND+2023[PDAT]&\
retmax=100&\
retmode=json"
```

### Use Case 3: Find Clinical Trials
```bash
# Search ClinicalTrials.gov database
curl "https://clinicaltrials.gov/api/v2/studies?query.cond=Alzheimer&query.intr=APOE&format=json"
```

## Storage Estimates

| Database | Storage Required | Format |
|----------|------------------|--------|
| PubMed abstracts (selective) | 1-10 GB | XML, JSON |
| PMC full-text (selective) | 10-100 GB | XML, PDF |
| OpenAlex snapshot | 200+ GB | JSON |
| Semantic Scholar corpus | 100+ GB | JSON |

## Data Quality Considerations

### PubMed Search Best Practices
- Use MeSH terms for standardized queries
- Combine with free-text for comprehensive coverage
- Filter by publication type (Review, Clinical Trial, etc.)
- Use date ranges to focus on recent literature
- Consider journal quality (impact factor)

### Evidence Ranking Factors
1. **Study design**: RCT > Cohort > Case-control > Case report
2. **Sample size**: Larger is generally better
3. **Journal impact**: Higher impact = more rigorous peer review
4. **Recency**: More recent may be more relevant
5. **Citation count**: Highly cited = influential

## Integration with Gene Platform

### Genetics → Literature
```
Genetic Variant (rs429358 - APOE ε4)
    ↓
PubMed Search: "rs429358 OR APOE ε4"
    ↓
Evidence: Alzheimer's risk, lipid metabolism, cardiovascular disease
```

### Traditional Medicine → Literature
```
TCM Herb (Ginkgo biloba)
    ↓
PubMed Search: "Ginkgo biloba AND clinical trial"
    ↓
Evidence: Cognitive function, neuroprotection
```

### Nutrition → Literature
```
Food Compound (Curcumin)
    ↓
PubMed Search: "Curcumin AND anti-inflammatory"
    ↓
Evidence: NF-κB inhibition, COX-2 suppression
```

## Update Strategy

### Real-Time (Tier 1)
- **PubMed**: Multiple updates daily, query fresh each time

### Periodic Refresh (Tier 2)
- **PMC**: Check weekly for new full-text articles
- **OpenAlex**: Monthly snapshot updates
- **Semantic Scholar**: Quarterly updates for citation graphs

## Navigation

- **Parent**: [Database Sources](../_index.md)
- **Related**: [Genetics](../genetics/_index.md), [Traditional Medicine](../traditional/_index.md), [Nutrition](../nutrition/_index.md), [Pathways](../pathways/_index.md)
