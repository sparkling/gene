---
id: schema-omim
title: "Online Mendelian Inheritance in Man (OMIM) Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, mendelian, genetics, rare-diseases, gene-disease, inheritance]
---

# Online Mendelian Inheritance in Man (OMIM) Schema Documentation

**Document ID:** SCHEMA-OMIM

---

## TL;DR

OMIM is the authoritative compendium of human genes and genetic phenotypes, providing expert-authored entries on all known Mendelian disorders and over 16,000 genes. It uses a structured MIM number system with prefixes indicating entry type (gene, phenotype, or combined).

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Gene Entries | 16,000+ | OMIM |
| Phenotype Entries | 9,000+ | OMIM |
| Gene-Phenotype Relationships | 7,200+ | OMIM |
| Allelic Variants | 29,000+ | OMIM |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | Text files, JSON (API) |
| API | Yes (requires registration) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| MIM Number (Gene) | `*[0-9]{6}` | *154700 |
| MIM Number (Phenotype) | `#[0-9]{6}` | #154700 |
| MIM Number (Combined) | `+[0-9]{6}` | +134610 |
| MIM Number (Unknown) | `%[0-9]{6}` | %176300 |

---

## References

See [Overview](./_index.md) for full details.
