---
id: ontology-owl
title: "OWL Ontologies"
type: ontology
parent: ../_index.md
last_updated: 2026-01-23
status: placeholder
tags: [owl, ontology, semantic-web, reasoning]
---

# OWL Ontologies

**Parent:** [Ontology](../_index.md)

OWL 2 (Web Ontology Language) definitions for formal class hierarchies, object properties, and reasoning axioms.

## Purpose

- Define formal class hierarchies (Gene, Variant, Compound, Disease, Pathway)
- Specify object and data properties with domains/ranges
- Enable automated reasoning and inference
- Support interoperability with external ontologies (GO, ChEBI, MONDO)

## Planned Files

| File | Description |
|------|-------------|
| `gene-platform.owl` | Core ontology for Gene Platform entities |
| `data-sources.owl` | Ontology for data source metadata |
| `imports/` | External ontology imports (GO, ChEBI, etc.) |

## External Ontologies

| Ontology | URI | Use |
|----------|-----|-----|
| Gene Ontology | http://purl.obolibrary.org/obo/go.owl | Biological processes, functions |
| ChEBI | http://purl.obolibrary.org/obo/chebi.owl | Chemical entities |
| MONDO | http://purl.obolibrary.org/obo/mondo.owl | Disease classification |
| SO | http://purl.obolibrary.org/obo/so.owl | Sequence features |
