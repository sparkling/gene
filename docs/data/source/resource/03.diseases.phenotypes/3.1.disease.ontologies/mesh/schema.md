---
id: schema-mesh
title: "Medical Subject Headings (MeSH) Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, vocabulary, medical-terms, nlm, pubmed, indexing, diseases]
---

# Medical Subject Headings (MeSH) Schema Documentation

**Document ID:** SCHEMA-MESH

---

## TL;DR

MeSH is the NLM controlled vocabulary thesaurus for indexing PubMed articles and NLM catalog books. It contains over 30,000 descriptors in a hierarchical tree structure with 16 major categories, providing standardized disease terminology for consistent literature searches.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Total Descriptors | ~30,000 | NLM |
| Disease Descriptors | ~5,000 | NLM |
| Supplementary Concepts | ~280,000 | NLM |
| Tree Categories | 16 | NLM |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | XML, ASCII, RDF/N-Triples, JSON-LD |
| API | Yes (SPARQL endpoint) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| MeSH UI | `D[0-9]{6}` | D003920 (Diabetes Mellitus) |
| Tree Number | `C[0-9]{2}(\.[0-9]{3})+` | C18.452.394.750 |
| MeSH Term | Free text | Diabetes Mellitus, Type 2 |
| Entry Term | Free text | Type 2 Diabetes |

---

## References

See [Overview](./_index.md) for full details.
