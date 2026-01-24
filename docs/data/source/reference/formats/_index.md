---
id: reference-formats
title: Data Format Specifications
description: Technical specifications for pathway databases, literature sources, and data extraction
category: reference
parent: ../_index.md
weight: 10
last_updated: 2026-01-23
status: active
tags: [formats, specifications, pathway, literature, xml, json]
---

# Data Format Specifications

Technical reference documentation for data formats used across the knowledge base.

## Contents

### [Pathway Formats](./pathway-formats.md)
Comprehensive technical specifications for pathway database formats including:
- Reactome Neo4j graph schema
- KEGG KGML format
- WikiPathways GPML format
- BioPAX (Biological Pathway Exchange)
- SBML (Systems Biology Markup Language)
- PSI-MI (Molecular Interactions)
- CX Format (NDEx)

### [Literature Formats](./literature-formats.md)
Data structures for research paper processing:
- PubMed XML format and DTD structure
- PMC JATS XML schema
- TypeScript interfaces for paper entities
- RuVector storage schemas for articles and chunks
- Embedding considerations for RAG systems

### [Text Extraction](./text-extraction.md)
Analysis of abstracts vs full-text approaches:
- Availability and coverage comparison
- Storage requirements analysis
- RAG quality metrics
- Processing complexity comparison
- Legal and licensing considerations
- Hybrid strategy recommendations

---

## Format Quick Reference

| Format | Primary Use | File Type | Parser |
|--------|-------------|-----------|--------|
| PubMed XML | Abstracts, metadata | XML | lxml |
| JATS XML | Full-text articles | XML | lxml |
| KGML | KEGG pathways | XML | BioPython |
| GPML | WikiPathways | XML | Custom |
| BioPAX | Pathway exchange | OWL/RDF | Paxtools |
| SBML | Systems biology | XML | libSBML |
| PSI-MI | Interactions | XML/MITAB | PSICQUIC |
| CX | Network exchange | JSON | ndex2 |

---

## Related Documentation

- [Schemas Index](../../operations/schemas/_index.md) - Database-specific schema documentation
- [Integration Guide](../../operations/integration/_index.md) - Data integration procedures
- [Downloads](../../operations/downloads/_index.md) - Data acquisition methods
