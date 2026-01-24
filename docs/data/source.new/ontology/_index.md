---
id: ontology
title: "Ontology"
type: ontology
parent: ../_index.md
last_updated: 2026-01-23
status: active
tags: [ontology, semantic-web, owl, shacl, skos, rdf]
---

# Ontology

**Parent:** [Data Sources](../_index.md)

Formal ontology definitions for the Gene Platform data catalog using semantic web standards. Enables machine-readable data integration, validation, and interoperability.

## Contents

| Folder | Format | Purpose |
|--------|--------|---------|
| [owl/](./owl/) | OWL 2 | Class hierarchies, properties, and formal axioms |
| [shacl/](./shacl/) | SHACL | Data validation shapes and constraints |
| [skos/](./skos/) | SKOS | Concept schemes and taxonomic vocabularies |

## Core Files

| File | Description |
|------|-------------|
| [ontology.yaml](./ontology.yaml) | Human-readable ontology definition (Three Worlds model) |

## Three Worlds Model

The ontology defines three primary domains:

| World | Focus | Color |
|-------|-------|-------|
| **Genetics** | Modern genomic and molecular biology | #4A90E2 |
| **Traditional Medicine** | Time-tested healing knowledge | #7ED321 |
| **Nutrition** | Evidence-based nutritional science | #F5A623 |

## Semantic Web Standards

- **OWL 2** (Web Ontology Language) - Formal class definitions and reasoning
- **SHACL** (Shapes Constraint Language) - Data validation and quality
- **SKOS** (Simple Knowledge Organization System) - Taxonomies and thesauri

## Related

- [taxonomy/](../taxonomy/) - Classification system documentation
- [reference/architecture/](../reference/architecture/) - Schema designs
