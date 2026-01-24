---
id: reference
title: Reference Documentation
description: Technical specifications, format definitions, and architectural documentation
category: reference
parent: ../_index.md
weight: 40
last_updated: 2026-01-23
status: active
tags: [reference, documentation, formats, architecture, specifications]
---

# Reference Documentation

This section contains technical reference documentation including format specifications, architectural designs, and schema definitions.

## Contents

### [Formats](./formats/_index.md)
Data format specifications and technical references for pathway databases, literature sources, and text extraction methods.

- **Pathway Formats** - Reactome, KEGG, WikiPathways, BioPAX, SBML technical specifications
- **Literature Formats** - PubMed XML, JATS, data structures for research papers
- **Text Extraction** - Abstracts vs full-text comparison and extraction strategies

### [Architecture](./architecture/_index.md)
System architecture documentation including schema designs, data models, and integration patterns.

- **Three Worlds Schema** - RuVector schema design for Genetics, Traditional Medicine, and Nutrition domains
- **Unified Schema Analysis** - Cross-database entity mapping and relationship patterns
- **Sample Data** - Example API responses and data structures from major databases

### [Database Catalogs](./catalogs/_index.md)
Comprehensive database catalogs with detailed specifications for 100+ databases across all categories.

- **Genetics Catalog** - dbSNP, ClinVar, gnomAD, PharmGKB and 20+ genetics databases
- **Compounds Catalog** - ChEMBL, DrugBank, PubChem, natural products databases
- **Traditional Medicine Catalog** - TCM, Ayurveda, Kampo, Western Herbal databases
- **Literature Catalog** - PubMed, PMC, OpenAlex, and literature pipeline design

---

## Quick Links

| Document | Description | Use Case |
|----------|-------------|----------|
| [Pathway Formats](./formats/pathway-formats.md) | Technical specs for pathway databases | Implementing pathway parsers |
| [Literature Formats](./formats/literature-formats.md) | PubMed/PMC data structures | Building literature pipeline |
| [Text Extraction](./formats/text-extraction.md) | Abstracts vs full-text analysis | RAG system design |
| [Three Worlds Schema](./architecture/three-worlds-schema.md) | RuVector database design | Database implementation |
| [Unified Schema](./architecture/unified-schema-analysis.md) | Entity mapping analysis | Data integration |
| [Sample Data](./architecture/sample-data.md) | API response examples | Testing and validation |

---

## Related Sections

- [Resource](../resource/_index.md) - Data source documentation by category
- [Guides](../guides/_index.md) - Integration and workflow guides
- [Operations](../operations/_index.md) - Data processing methodology
- [Research](../research/_index.md) - Research priorities and gap analysis
