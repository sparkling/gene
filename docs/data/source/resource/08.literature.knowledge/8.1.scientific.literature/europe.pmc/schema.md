---
id: schema-europe-pmc
title: "Europe PMC Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, database, literature, biomedical]
---

# Europe PMC - Schema Documentation

## TL;DR

Europe PMC provides access to biomedical literature with rich metadata via REST API and OAI-PMH. The API returns XML or JSON responses containing article metadata, annotations, citations, and full-text links.

## Primary Data Structures

### Article Record

```json
{
  "id": "12345678",
  "source": "MED",
  "pmid": "12345678",
  "pmcid": "PMC1234567",
  "doi": "10.1000/example",
  "title": "Article title",
  "authorString": "Author A, Author B",
  "journalTitle": "Journal Name",
  "pubYear": "2024",
  "abstractText": "Full abstract text...",
  "isOpenAccess": "Y",
  "citedByCount": 42,
  "firstPublicationDate": "2024-01-15"
}
```

### Author Object

```json
{
  "firstName": "John",
  "lastName": "Smith",
  "initials": "JS",
  "authorId": {
    "type": "ORCID",
    "value": "0000-0002-1825-0097"
  },
  "affiliation": "University Name, Department"
}
```

### Grant Information

```json
{
  "grantId": "BB/X00001X/1",
  "agency": "BBSRC",
  "acronym": "BB",
  "orderIn": 1
}
```

## Core Fields

| Field | Type | Description |
|-------|------|-------------|
| id | string | Internal Europe PMC ID |
| source | string | MED, PMC, PAT, AGR, etc. |
| pmid | string | PubMed identifier |
| pmcid | string | PMC identifier (PMC prefix) |
| doi | string | Digital Object Identifier |
| title | string | Article title |
| authorString | string | Formatted author list |
| authorList | array | Structured author objects |
| journalInfo | object | Journal metadata |
| pubYear | string | Publication year |
| abstractText | string | Abstract content |
| affiliation | string | First author affiliation |
| language | string | ISO language code |
| pubTypeList | array | Publication types (MeSH) |
| meshHeadingList | array | MeSH terms |
| keywordList | array | Author keywords |
| chemicalList | array | Chemical substances |
| grantsList | array | Funding information |

## Source Types

| Source | Code | Description |
|--------|------|-------------|
| MEDLINE | MED | PubMed abstracts |
| PMC | PMC | Full-text articles |
| Patents | PAT | Patent documents |
| Agricultural | AGR | AGRIS/AGRICOLA |
| Preprints | PPR | bioRxiv, medRxiv |
| Theses | ETH | PhD theses |
| CBA | CBA | Chinese Biological Abstracts |

## Citation Data

```json
{
  "citation": {
    "id": "12345679",
    "source": "MED",
    "citationType": "REFERENCES",
    "title": "Cited article title",
    "pubYear": "2020"
  }
}
```

## Annotations (SciLite)

```json
{
  "annotation": {
    "type": "Gene_Proteins",
    "exact": "TP53",
    "prefix": "the tumor suppressor ",
    "postfix": " is frequently mutated",
    "section": "Abstract",
    "tags": [
      {
        "name": "TP53",
        "uri": "https://identifiers.org/ncbigene:7157"
      }
    ]
  }
}
```

### Annotation Types

| Type | Description |
|------|-------------|
| Gene_Proteins | Gene and protein mentions |
| Diseases | Disease terms (EFO) |
| Chemicals | Chemical compounds |
| Organisms | Species mentions |
| GO_Terms | Gene Ontology terms |
| Accession_Numbers | Database accessions |

## Full-Text XML Structure (JATS)

```xml
<article>
  <front>
    <article-meta>
      <article-id pub-id-type="pmid">12345678</article-id>
      <title-group>
        <article-title>Title</article-title>
      </title-group>
      <abstract>
        <p>Abstract text...</p>
      </abstract>
    </article-meta>
  </front>
  <body>
    <sec>
      <title>Introduction</title>
      <p>Section content...</p>
    </sec>
  </body>
  <back>
    <ref-list>
      <ref id="R1">...</ref>
    </ref-list>
  </back>
</article>
```

## API Response Wrapper

```json
{
  "version": "6.9",
  "hitCount": 1234,
  "request": {
    "queryString": "TP53",
    "resultType": "lite",
    "cursorMark": "*",
    "pageSize": 25,
    "sort": ""
  },
  "resultList": {
    "result": [...]
  }
}
```

## Result Types

| Type | Description |
|------|-------------|
| idlist | IDs only |
| lite | Core fields |
| core | Full metadata |

## See Also

- [Download Documentation](./download.md)
- [Europe PMC API](https://www.ebi.ac.uk/europepmc/webservices/rest)
