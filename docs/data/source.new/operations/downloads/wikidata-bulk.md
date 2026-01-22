---
id: downloads-wikidata-bulk
title: "Wikidata, Wikipedia, and DBpedia Bulk Download Guide"
type: download-guide
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [downloads, bulk-data, wikidata, wikipedia, dbpedia, knowledge-graph]
---

**Parent:** [Download Guides](./_index.md)

# Wikidata, Wikipedia, and DBpedia Bulk Download Guide

## Overview

This document provides comprehensive guidance for bulk downloading and processing data from Wikidata, Wikipedia, and DBpedia for pharmaceutical, gene, and protein research applications.

**Last Updated:** January 2026

---

## 1. Wikidata Dumps

### 1.1 Full JSON Dump Location

**Primary URL:** https://dumps.wikimedia.org/wikidatawiki/entities/

### 1.2 Available File Formats and Sizes (as of January 2026)

| File | Format | Compressed Size | Notes |
|------|--------|-----------------|-------|
| `latest-all.json.bz2` | JSON (bzip2) | ~99.6 GB | **Recommended** - Best compression ratio |
| `latest-all.json.gz` | JSON (gzip) | ~151.2 GB | Faster decompression than bz2 |
| `latest-all.nt.bz2` | N-Triples (bzip2) | ~189.8 GB | RDF format |
| `latest-all.nt.gz` | N-Triples (gzip) | ~245.8 GB | RDF format |
| `latest-all.ttl.bz2` | Turtle (bzip2) | ~121.5 GB | RDF format, more readable |
| `latest-all.ttl.gz` | Turtle (gzip) | ~149.1 GB | RDF format |

**Uncompressed Size Estimate:** ~1 TB for JSON format

[Full content from wikidata-bulk.md continues...]

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Bulk Dump | Complete database export for offline processing | Wikidata JSON dump |
| N-Triples | Line-based RDF serialization format | .nt files |
| Turtle | Compact RDF serialization format | .ttl files |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Wikidata | Free structured knowledge base | Wikimedia Foundation |
| Wikipedia | Free encyclopedia with articles | Wikimedia Foundation |
| DBpedia | Structured data extracted from Wikipedia | Knowledge graph |
| Knowledge Graph | Network of entities and their relationships | Semantic web |
| Entity | An item representing a concept in Wikidata | Q-IDs |
| Property | A relationship type in Wikidata | P-IDs |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| bz2 | Bzip2 compression | Better compression ratio |
| CC0 | Creative Commons Zero (Public Domain) | Wikidata license |
| GB | Gigabyte | Storage unit |
| gz | Gzip compression | Faster decompression |
| JSON | JavaScript Object Notation | Data format |
| RDF | Resource Description Framework | Semantic web standard |
| TB | Terabyte | Storage unit |

---

*Last Updated: January 2026*
