# 8.1 Scientific Literature - Data Dictionary

**Subcategory ID:** 8.1
**Subcategory Name:** Scientific Literature
**Data Sources:** PubMed, PubMed Central (PMC), Europe PMC, OpenAlex, Semantic Scholar

## Overview

This subcategory integrates data from major scientific literature databases, providing unified access to publication metadata, abstracts, citations, and full-text content where available.

---

## Unified Fields

### Core Identifiers

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `article_id` | string | Yes | Unique identifier for the article/publication record | `12345678` |
| `doi` | string | No | Digital Object Identifier for the publication | `10.1038/s41586-019-1666-5` |
| `pmid` | string | No | PubMed identifier | `31534227` |
| `pmcid` | string | No | PubMed Central identifier | `PMC6868014` |

**Source Mappings - Core Identifiers:**

| Field | PubMed | PMC | Europe PMC | OpenAlex | Semantic Scholar |
|-------|--------|-----|------------|----------|------------------|
| `article_id` | PMID | article-id[@pub-id-type='pmc'] | id | id | paperId |
| `doi` | ArticleId[@IdType='doi'] | article-id[@pub-id-type='doi'] | doi | doi | externalIds.DOI |
| `pmid` | PMID | article-id[@pub-id-type='pmid'] | pmid | - | externalIds.PubMed |
| `pmcid` | ArticleId[@IdType='pmc'] | article-id[@pub-id-type='pmc'] | pmcid | - | externalIds.PubMedCentral |

---

### Bibliographic Information

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `title` | string | Yes | Title of the article/publication | `Deep learning for genomic analysis` |
| `abstract` | string | No | Abstract or summary of the article | `Full abstract text here...` |
| `publication_date` | date | No | Date of publication (YYYY-MM-DD) | `2024-01-15` |
| `publication_year` | integer | No | Year of publication | `2024` |
| `journal_title` | string | No | Name of the journal or publication venue | `Nature` |
| `issn` | array[string] | No | International Standard Serial Number(s) | `["0028-0836", "1476-4687"]` |
| `volume` | string | No | Journal volume number | `615` |
| `issue` | string | No | Journal issue number | `7950` |
| `pages` | string | No | Page range in the publication | `100-110` |
| `language` | string | No | Language of the publication (ISO code) | `eng` |

**Source Mappings - Bibliographic Information:**

| Field | PubMed | PMC | Europe PMC | OpenAlex | Semantic Scholar |
|-------|--------|-----|------------|----------|------------------|
| `title` | ArticleTitle | article-title | title | title | title |
| `abstract` | AbstractText | abstract | abstractText | abstract_inverted_index | abstract |
| `publication_date` | PubDate | pub-date | firstPublicationDate | publication_date | publicationDate |
| `publication_year` | PubDate/Year | pub-date/year | pubYear | publication_year | year |
| `journal_title` | Journal/Title | journal-title | journalTitle | primary_location.source.display_name | venue |
| `issn` | ISSN | issn | - | primary_location.source.issn | publicationVenue.issn |
| `volume` | JournalIssue/Volume | - | - | - | journal.volume |
| `issue` | JournalIssue/Issue | - | - | - | - |
| `pages` | MedlinePgn | - | - | - | journal.pages |
| `language` | Language | - | language | - | - |

---

### Authors and Affiliations

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `authors` | array[object] | No | List of authors who contributed to the work | See structure below |
| `affiliations` | array[object] | No | Author institutional affiliations | See structure below |

**Author Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `last_name` | string | Author's family/surname |
| `first_name` | string | Author's given name |
| `full_name` | string | Complete author name |
| `orcid` | string | ORCID identifier (0000-0000-0000-0000) |
| `position` | string | Position in author list (first, middle, last) |
| `is_corresponding` | boolean | Whether author is corresponding author |

**Affiliation Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `name` | string | Institution name |
| `ror_id` | string | Research Organization Registry ID |
| `country` | string | Country of institution |

**Source Mappings - Authors:**

| Field | PubMed | PMC | Europe PMC | OpenAlex | Semantic Scholar |
|-------|--------|-----|------------|----------|------------------|
| `authors` | AuthorList/Author | contrib-group/contrib | authorList | authorships | authors |
| `affiliations` | AffiliationInfo/Affiliation | aff | affiliation | authorships.institutions | affiliations |

---

### Classification and Indexing

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `publication_types` | array[string] | No | Type or category of the publication | `["Journal Article", "Review"]` |
| `mesh_terms` | array[object] | No | Medical Subject Headings (MeSH) indexing terms | See structure below |
| `keywords` | array[string] | No | Author-supplied or controlled keywords | `["TP53", "cancer", "gene mutation"]` |
| `concepts` | array[object] | No | Topic concepts or fields of study | See structure below |

**MeSH Term Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `descriptor_name` | string | MeSH descriptor name |
| `descriptor_ui` | string | MeSH descriptor unique identifier |
| `qualifier_name` | string | MeSH qualifier/subheading |
| `is_major_topic` | boolean | Whether this is a major topic |

**Concept Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `id` | string | Concept identifier |
| `name` | string | Concept name |
| `level` | integer | Hierarchical level (0=broad, 5=specific) |
| `score` | number | Relevance score (0-1) |
| `wikidata_qid` | string | Wikidata entity ID |

**Source Mappings - Classification:**

| Field | PubMed | PMC | Europe PMC | OpenAlex | Semantic Scholar |
|-------|--------|-----|------------|----------|------------------|
| `publication_types` | PublicationType | article-type | pubTypeList | type | publicationTypes |
| `mesh_terms` | MeshHeadingList | - | meshHeadingList | - | - |
| `keywords` | KeywordList | kwd-group | keywordList | - | - |
| `concepts` | - | - | - | concepts | fieldsOfStudy |

---

### Citations and References

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `citation_count` | integer | No | Number of articles citing this work | `42` |
| `references` | array[object] | No | Works cited by this article | See structure below |

**Reference Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `pmid` | string | PubMed ID of referenced work |
| `doi` | string | DOI of referenced work |
| `citation` | string | Formatted citation string |

**Source Mappings - Citations:**

| Field | PubMed | PMC | Europe PMC | OpenAlex | Semantic Scholar |
|-------|--------|-----|------------|----------|------------------|
| `citation_count` | - | - | citedByCount | cited_by_count | citations (count) |
| `references` | ReferenceList | ref-list | - | referenced_works | references |

---

### Open Access Information

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `is_open_access` | boolean | No | Whether the article is openly accessible | `true` |
| `open_access_status` | string | No | Type of open access (gold, green, hybrid, bronze) | `gold` |
| `full_text_url` | string | No | URL to full text content | `https://www.nature.com/articles/....pdf` |
| `license` | string | No | Content license type | `cc-by` |

**Source Mappings - Open Access:**

| Field | PubMed | PMC | Europe PMC | OpenAlex | Semantic Scholar |
|-------|--------|-----|------------|----------|------------------|
| `is_open_access` | - | - | isOpenAccess | is_oa | isOpenAccess |
| `open_access_status` | - | - | - | open_access.oa_status | openAccessPdf.status |
| `full_text_url` | - | - | - | primary_location.pdf_url | openAccessPdf.url |
| `license` | - | license | - | primary_location.license | - |

---

### Funding Information

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `grant_information` | array[object] | No | Funding and grant details | See structure below |

**Grant Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `grant_id` | string | Grant/award number |
| `agency` | string | Funding agency name |
| `country` | string | Country of funding agency |

**Source Mappings - Funding:**

| Field | PubMed | PMC | Europe PMC | OpenAlex | Semantic Scholar |
|-------|--------|-----|------------|----------|------------------|
| `grant_information` | GrantList | - | grantsList | - | - |

---

## Source-Specific Fields

### PubMed-Specific Fields

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `citation_status` | string | MEDLINE indexing status | `MEDLINE`, `PubMed-not-MEDLINE`, `In-Process` |
| `date_completed` | date | Date record was completed | `2024-01-20` |
| `date_revised` | date | Date record was last revised | `2024-02-15` |
| `nlm_unique_id` | string | NLM internal journal identifier | `0410462` |
| `structured_abstract` | object | Abstract with labeled sections | `{background, methods, results, conclusions}` |
| `publication_history` | array | Dates for received, accepted, pubmed | `[{date, status}]` |

### PubMed Central (PMC)-Specific Fields

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `manuscript_id` | string | NIH Manuscript Submission ID | `NIHMS123456` |
| `publisher_id` | string | Publisher-assigned identifier | `pone.0123456` |
| `article_sections` | array | Structured body sections | `["Introduction", "Results", "Methods"]` |
| `figures` | array[object] | Figure elements with label, caption, href | `[{label, caption, href}]` |
| `tables` | array[object] | Table elements with structured data | `[{label, caption, data}]` |
| `supplementary_materials` | array[object] | Supplementary files and data | `[{title, url, type}]` |
| `named_content` | array | Tagged entities (gene, species, compound) | `["TP53", "Homo sapiens"]` |
| `acknowledgements` | string | Acknowledgements section content | `We thank...` |

### Europe PMC-Specific Fields

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `source_type` | string | Source code | `MED`, `PMC`, `PAT`, `AGR`, `PPR`, `ETH`, `CBA` |
| `annotations` | array[object] | SciLite annotations | `[{type, text, prefix, postfix, section}]` |
| `chemical_list` | array | Chemical substances mentioned | `["metformin", "insulin"]` |

**Annotation Types:** Gene_Proteins, Diseases, Chemicals, Organisms, GO_Terms, Accession_Numbers

### OpenAlex-Specific Fields

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `openalex_id` | string | OpenAlex Work ID | `https://openalex.org/W2741809807` |
| `mag_id` | string | Microsoft Academic Graph ID | `2741809807` |
| `is_retracted` | boolean | Whether article has been retracted | `false` |
| `is_paratext` | boolean | Whether content is paratext | `false` |
| `abstract_inverted_index` | object | Compact word-position representation | `{word: [positions]}` |
| `related_works` | array[string] | Similar/related work IDs | `["W123...", "W456..."]` |
| `institution_ror` | string | Research Organization Registry ID | `https://ror.org/03vek6s52` |
| `concept_level` | integer | Hierarchical level of concept | `0` (broad) to `5` (specific) |
| `wikidata_qid` | string | Wikidata entity ID for concepts | `Q11190` |
| `author_position` | string | Position in author list | `first`, `middle`, `last` |
| `is_corresponding` | boolean | Whether author is corresponding | `true` |
| `h_index` | integer | Author h-index metric | `45` |
| `i10_index` | integer | Author i10-index metric | `120` |

### Semantic Scholar-Specific Fields

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `corpus_id` | integer | Semantic Scholar corpus identifier | `123456789` |
| `arxiv_id` | string | ArXiv preprint identifier | `2001.01234` |
| `dblp_id` | string | DBLP bibliography identifier | `conf/nips/SmithJ20` |
| `mag_author_id` | string | Microsoft Academic Graph author ID | `2157025439` |
| `influential_citation_count` | integer | Count of influential citations | `15` |
| `is_influential` | boolean | Whether citation is influential | `true` |
| `citation_contexts` | array[string] | Text passages where citation appears | `["As shown by Smith et al..."]` |
| `citation_intents` | array[string] | Purpose of citation | `["methodology", "background", "result"]` |
| `tldr` | string | AI-generated one-sentence summary | `This paper presents...` |
| `tldr_model` | string | Version of TLDR model used | `tldr@v2.0.0` |
| `embedding` | array[number] | 768-dimensional SPECTER vector | `[0.123, -0.456, ...]` |
| `s2_fields_of_study` | array[object] | Semantic Scholar classified fields | `[{category, source}]` |
| `match_score` | number | Search relevance score | `0.95` |

---

## Data Source Metadata

| Field | Data Type | Description |
|-------|-----------|-------------|
| `_source.primary_source` | string | Name of the primary data source |
| `_source.source_id` | string | Original ID in the source system |
| `_source.extraction_date` | date | Date when data was extracted |
| `_source.source_version` | string | Version of the source data/API |
