---
id: schema-wikipedia
title: "Wikipedia Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: final
tags: [schema, database, knowledge-base, encyclopedia]
---

# Wikipedia - Schema Documentation

## TL;DR

Wikipedia content is available through the MediaWiki API and REST API, providing article text, metadata, revisions, and links. Data can also be accessed through structured extractions like DBpedia.

## Article Page Object

```json
{
  "pageid": 1234567,
  "ns": 0,
  "title": "TP53",
  "contentmodel": "wikitext",
  "pagelanguage": "en",
  "touched": "2024-01-15T10:30:00Z",
  "lastrevid": 1234567890,
  "length": 45678,
  "fullurl": "https://en.wikipedia.org/wiki/TP53",
  "editurl": "https://en.wikipedia.org/w/index.php?title=TP53&action=edit"
}
```

## REST API Response

### Page Summary

```json
{
  "type": "standard",
  "title": "TP53",
  "displaytitle": "TP53",
  "namespace": {
    "id": 0,
    "text": ""
  },
  "wikibase_item": "Q21173022",
  "titles": {
    "canonical": "TP53",
    "normalized": "TP53",
    "display": "TP53"
  },
  "pageid": 1234567,
  "thumbnail": {
    "source": "https://upload.wikimedia.org/...",
    "width": 220,
    "height": 200
  },
  "originalimage": {
    "source": "https://upload.wikimedia.org/...",
    "width": 800,
    "height": 727
  },
  "lang": "en",
  "dir": "ltr",
  "revision": "1234567890",
  "tid": "abc123-def456",
  "timestamp": "2024-01-15T10:30:00Z",
  "description": "Human gene and protein",
  "description_source": "local",
  "content_urls": {
    "desktop": {
      "page": "https://en.wikipedia.org/wiki/TP53",
      "revisions": "https://en.wikipedia.org/wiki/TP53?action=history",
      "edit": "https://en.wikipedia.org/wiki/TP53?action=edit",
      "talk": "https://en.wikipedia.org/wiki/Talk:TP53"
    },
    "mobile": {
      "page": "https://en.m.wikipedia.org/wiki/TP53"
    }
  },
  "extract": "TP53 (tumor protein P53) is a gene that codes for a protein...",
  "extract_html": "<p><b>TP53</b> (tumor protein P53) is a gene...</p>"
}
```

## MediaWiki API Responses

### Query Page Content

```json
{
  "query": {
    "pages": {
      "1234567": {
        "pageid": 1234567,
        "ns": 0,
        "title": "TP53",
        "revisions": [
          {
            "revid": 1234567890,
            "parentid": 1234567889,
            "user": "Username",
            "timestamp": "2024-01-15T10:30:00Z",
            "comment": "Edit summary",
            "slots": {
              "main": {
                "contentmodel": "wikitext",
                "contentformat": "text/x-wiki",
                "content": "{{Infobox gene}}\n'''TP53'''..."
              }
            }
          }
        ]
      }
    }
  }
}
```

### Page Categories

```json
{
  "query": {
    "pages": {
      "1234567": {
        "categories": [
          {"ns": 14, "title": "Category:Tumor suppressor genes"},
          {"ns": 14, "title": "Category:Genes on human chromosome 17"},
          {"ns": 14, "title": "Category:Transcription factors"}
        ]
      }
    }
  }
}
```

### Page Links

```json
{
  "query": {
    "pages": {
      "1234567": {
        "links": [
          {"ns": 0, "title": "Cancer"},
          {"ns": 0, "title": "DNA repair"},
          {"ns": 0, "title": "Apoptosis"}
        ]
      }
    }
  }
}
```

### External Links

```json
{
  "query": {
    "pages": {
      "1234567": {
        "extlinks": [
          {"url": "https://www.ncbi.nlm.nih.gov/gene/7157"},
          {"url": "https://www.uniprot.org/uniprot/P04637"},
          {"url": "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000141510"}
        ]
      }
    }
  }
}
```

## Namespace Codes

| NS | Name | Description |
|----|------|-------------|
| 0 | Main | Articles |
| 1 | Talk | Article discussion |
| 2 | User | User pages |
| 6 | File | Media files |
| 10 | Template | Templates |
| 14 | Category | Categories |

## Wikitext Infobox (Gene Example)

```wikitext
{{Infobox gene
| Name = TP53
| HGNCid = 11998
| Symbol = TP53
| AltSymbols = p53, LFS1
| EntrezGene = 7157
| OMIM = 191170
| RefSeq = NM_000546
| UniProt = P04637
| ECnumber =
| Chromosome = 17
| Arm = p
| Band = 13.1
| LocusSupplementaryData =
}}
```

## Parsed HTML Response

```json
{
  "parse": {
    "title": "TP53",
    "pageid": 1234567,
    "text": {
      "*": "<div class=\"mw-parser-output\">..."
    },
    "sections": [
      {"toclevel": 1, "level": "2", "line": "Function", "number": "1"},
      {"toclevel": 1, "level": "2", "line": "Structure", "number": "2"},
      {"toclevel": 1, "level": "2", "line": "Clinical significance", "number": "3"}
    ],
    "categories": [...],
    "links": [...],
    "templates": [...],
    "images": [...],
    "externallinks": [...]
  }
}
```

## DBpedia Structured Data

```turtle
@prefix dbr: <http://dbpedia.org/resource/> .
@prefix dbo: <http://dbpedia.org/ontology/> .

dbr:P53 a dbo:Gene ;
  rdfs:label "p53"@en ;
  dbo:abstract "p53, also known as..."@en ;
  dbo:chromosome "17" ;
  dbo:entrezgene "7157" ;
  dbo:hgncid "11998" ;
  dbo:omim "191170" ;
  dbo:refseq "NM_000546" ;
  dbo:symbol "TP53" ;
  dbo:uniprot "P04637" ;
  owl:sameAs <http://www.wikidata.org/entity/Q21173022> .
```

## Wikidata QID Mapping

Each Wikipedia article links to a Wikidata item:

| Wikipedia | Wikidata QID |
|-----------|--------------|
| TP53 | Q21173022 |
| Cancer | Q12078 |
| BRCA1 | Q17853272 |

## Dump File Formats

| Format | Description |
|--------|-------------|
| XML | Full page content with revisions |
| SQL | Database tables |
| JSON | Entity dumps (Wikidata) |

## See Also

- [Download Documentation](./download.md)
- [MediaWiki API](https://www.mediawiki.org/wiki/API:Main_page)
