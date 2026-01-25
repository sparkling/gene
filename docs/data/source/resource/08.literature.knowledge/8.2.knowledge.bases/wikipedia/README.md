---
id: wikipedia
title: "Wikipedia"
type: data-source
category: literature
subcategory: knowledge.bases
parent: ../README.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [encyclopedia, text, descriptions, multilingual, open-content]
---

# Wikipedia

**Category:** [Literature](../../../README.md) > [Knowledge Bases](../README.md)

## Overview

Wikipedia is the free encyclopedia that anyone can edit, containing millions of articles across hundreds of languages. For biomedical applications, Wikipedia provides accessible descriptions of genes, proteins, diseases, drugs, and biological concepts.

While not a primary scientific resource, Wikipedia articles often provide excellent contextual information and are maintained by domain experts. Many bioinformatics resources use Wikipedia as a source for entity descriptions and definitions.

Wikipedia is valuable for generating human-readable descriptions, training NLP models, and providing context for technical database entries.

## Key Statistics

| Metric | Value |
|--------|-------|
| English Articles | 6.8M+ |
| Total Articles (all languages) | 63M+ |
| Languages | 300+ |
| Active Editors | 125,000+ |
| Last Update | Continuous |

## Primary Use Cases

1. Entity description retrieval
2. NLP training data
3. Definition lookup
4. Cross-language content
5. Educational resource

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Page Title | URL-encoded | TP53 |
| Page ID | Numeric | 1234567 |
| Wikidata QID | `Q[0-9]+` | Q21173022 |
| Wikipedia URL | Full URL | en.wikipedia.org/wiki/TP53 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.wikipedia.org | Browse and read |
| REST API | https://en.wikipedia.org/api/rest_v1 | Content retrieval |
| MediaWiki API | https://en.wikipedia.org/w/api.php | Full API |
| Dumps | https://dumps.wikimedia.org | Full database |
| DBpedia | https://dbpedia.org | Structured extraction |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY-SA 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [Wikidata](../wikidata/README.md) - Structured data
