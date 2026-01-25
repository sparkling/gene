---
id: schema-openalex
title: "OpenAlex Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, database, literature, scholarly]
---

# OpenAlex - Schema Documentation

## TL;DR

OpenAlex provides a comprehensive scholarly knowledge graph with five core entity types: Works, Authors, Sources (venues), Institutions, and Concepts. All data is available via REST API and bulk snapshot downloads.

## Core Entity Types

### Work (Scholarly Publication)

```json
{
  "id": "https://openalex.org/W2741809807",
  "doi": "https://doi.org/10.1038/s41586-019-1666-5",
  "title": "Article title here",
  "display_name": "Article title here",
  "publication_year": 2019,
  "publication_date": "2019-10-30",
  "type": "journal-article",
  "cited_by_count": 1234,
  "is_retracted": false,
  "is_paratext": false,
  "is_oa": true,
  "primary_location": {
    "source": {
      "id": "https://openalex.org/S137773608",
      "display_name": "Nature"
    },
    "landing_page_url": "https://www.nature.com/articles/...",
    "pdf_url": "https://www.nature.com/articles/....pdf",
    "is_oa": true,
    "version": "publishedVersion",
    "license": "cc-by"
  },
  "authorships": [...],
  "concepts": [...],
  "referenced_works": [...],
  "related_works": [...],
  "abstract_inverted_index": {...}
}
```

### Author

```json
{
  "id": "https://openalex.org/A5023888391",
  "orcid": "https://orcid.org/0000-0002-1825-0097",
  "display_name": "Jane Doe",
  "display_name_alternatives": ["J Doe", "Jane A. Doe"],
  "works_count": 156,
  "cited_by_count": 12345,
  "most_cited_work": "https://openalex.org/W2741809807",
  "summary_stats": {
    "2yr_mean_citedness": 4.5,
    "h_index": 25,
    "i10_index": 45
  },
  "affiliations": [...],
  "last_known_institution": {...},
  "x_concepts": [...]
}
```

### Source (Venue/Journal)

```json
{
  "id": "https://openalex.org/S137773608",
  "issn_l": "0028-0836",
  "issn": ["0028-0836", "1476-4687"],
  "display_name": "Nature",
  "publisher": "Springer Nature",
  "type": "journal",
  "works_count": 412345,
  "cited_by_count": 23456789,
  "is_oa": false,
  "is_in_doaj": false,
  "homepage_url": "https://www.nature.com/",
  "x_concepts": [...]
}
```

### Institution

```json
{
  "id": "https://openalex.org/I136199984",
  "ror": "https://ror.org/03vek6s52",
  "display_name": "Harvard University",
  "country_code": "US",
  "type": "education",
  "homepage_url": "https://www.harvard.edu/",
  "works_count": 567890,
  "cited_by_count": 45678901,
  "associated_institutions": [...],
  "x_concepts": [...]
}
```

### Concept (Topic)

```json
{
  "id": "https://openalex.org/C71924100",
  "wikidata": "https://www.wikidata.org/entity/Q11190",
  "display_name": "Medicine",
  "level": 0,
  "description": "field dealing with health and healing",
  "works_count": 12345678,
  "cited_by_count": 234567890,
  "ancestors": [...],
  "related_concepts": [...]
}
```

## Authorship Object

```json
{
  "author_position": "first",
  "author": {
    "id": "https://openalex.org/A5023888391",
    "display_name": "Jane Doe",
    "orcid": "https://orcid.org/0000-0002-1825-0097"
  },
  "institutions": [
    {
      "id": "https://openalex.org/I136199984",
      "display_name": "Harvard University",
      "ror": "https://ror.org/03vek6s52",
      "country_code": "US",
      "type": "education"
    }
  ],
  "raw_affiliation_string": "Department of Medicine, Harvard University",
  "is_corresponding": true
}
```

## Concept Tagging

```json
{
  "id": "https://openalex.org/C71924100",
  "wikidata": "https://www.wikidata.org/entity/Q11190",
  "display_name": "Medicine",
  "level": 0,
  "score": 0.89
}
```

### Concept Levels

| Level | Scope | Example |
|-------|-------|---------|
| 0 | Broad discipline | Medicine, Biology |
| 1 | Major field | Genetics, Neuroscience |
| 2 | Subfield | Molecular biology |
| 3-5 | Specific topics | DNA repair, p53 |

## Abstract Inverted Index

Compact representation of abstract text:

```json
{
  "abstract_inverted_index": {
    "The": [0, 45],
    "study": [1, 67],
    "demonstrates": [2],
    "that": [3, 89],
    "p53": [4, 12, 78]
  }
}
```

## Work Types

| Type | Description |
|------|-------------|
| journal-article | Peer-reviewed article |
| book-chapter | Book chapter |
| proceedings-article | Conference paper |
| book | Complete book |
| dataset | Data publication |
| dissertation | Thesis/dissertation |
| preprint | Non-peer-reviewed |
| review | Review article |

## Open Access Status

```json
{
  "open_access": {
    "is_oa": true,
    "oa_status": "gold",
    "oa_url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC..."
  }
}
```

### OA Status Values

| Status | Description |
|--------|-------------|
| gold | Published OA in OA venue |
| green | Repository version |
| hybrid | OA in subscription journal |
| bronze | Free to read (no license) |
| closed | Not OA |

## Identifiers Cross-Reference

| Entity | External IDs |
|--------|--------------|
| Work | DOI, PMID, PMCID, MAG ID |
| Author | ORCID, MAG Author ID |
| Source | ISSN, ISSN-L |
| Institution | ROR, GRID, Wikidata |
| Concept | Wikidata QID |

## See Also

- [Download Documentation](./download.md)
- [OpenAlex API Documentation](https://docs.openalex.org/)
