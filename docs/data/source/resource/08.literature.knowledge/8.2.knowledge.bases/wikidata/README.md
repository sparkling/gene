---
id: wikidata
title: "Wikidata"
type: data-source
category: literature
subcategory: knowledge.bases
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [knowledge-graph, linked-data, semantic-web, identifiers, open-data]
---

# Wikidata

**Category:** [Literature](../../../_index.md) > [Knowledge Bases](../_index.md)

## Overview

Wikidata is a free, collaborative, multilingual knowledge base operated by the Wikimedia Foundation. It serves as central storage for structured data of Wikipedia and other Wikimedia projects, but is also widely used as a general-purpose knowledge graph.

Wikidata contains items (entities) connected by statements (claims) with references and qualifiers. It has extensive coverage of biomedical entities including genes, proteins, diseases, drugs, and publications, with cross-references to specialized databases.

Wikidata is essential for entity resolution, knowledge graph construction, and linking identifiers across databases.

## Key Statistics

| Metric | Value |
|--------|-------|
| Items | 110M+ |
| Statements | 1.7B+ |
| Properties | 12,000+ |
| Languages | 300+ |
| Last Update | Continuous |

## Primary Use Cases

1. Entity resolution and linking
2. Knowledge graph construction
3. Cross-database identifier mapping
4. Multilingual label retrieval
5. Semantic web applications

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| QID (Item) | `Q[0-9]+` | Q21173022 (TP53) |
| PID (Property) | `P[0-9]+` | P352 (UniProt ID) |
| External IDs | Via properties | UniProt: P04637 |
| URI | https://www.wikidata.org/entity/ | Q21173022 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.wikidata.org | Browse and edit |
| SPARQL | https://query.wikidata.org | RDF queries |
| REST API | https://www.wikidata.org/w/api.php | MediaWiki API |
| Dumps | https://dumps.wikimedia.org/wikidatawiki | Full database |
| Linked Data | Content negotiation | RDF/JSON-LD |

## License

| Aspect | Value |
|--------|-------|
| License | CC0 (Public Domain) |
| Commercial Use | Yes |
| Attribution | Not required (appreciated) |

## See Also

- [Schema Documentation](./schema.md)
- [Wikipedia](../wikipedia/_index.md) - Textual content
- [UniProt ID Mapping](../../8.3.identifier.mapping/uniprot.id.mapping/_index.md) - Protein IDs
