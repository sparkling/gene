---
id: ontology-skos
title: "SKOS Vocabularies"
type: ontology
parent: ../_index.md
last_updated: 2026-01-23
status: placeholder
tags: [skos, vocabulary, taxonomy, thesaurus, concepts]
---

# SKOS Vocabularies

**Parent:** [Ontology](../_index.md)

SKOS (Simple Knowledge Organization System) concept schemes for taxonomies, controlled vocabularies, and thesauri.

## Purpose

- Define hierarchical concept schemes (broader/narrower)
- Manage preferred and alternative labels
- Support multilingual terminology
- Enable faceted classification

## Planned Files

| File | Description |
|------|-------------|
| `categories.skos.ttl` | 9-category data source classification |
| `data-types.skos.ttl` | Data type vocabulary (API, bulk, SPARQL) |
| `licenses.skos.ttl` | License classification scheme |
| `domains.skos.ttl` | Health domain concepts |

## Concept Scheme Structure

```turtle
# Example: Category concept scheme
ex:CategoryScheme a skos:ConceptScheme ;
    skos:prefLabel "Data Source Categories"@en ;
    skos:hasTopConcept ex:Genetics, ex:Compounds, ex:Diseases .

ex:Genetics a skos:Concept ;
    skos:prefLabel "Genetics & Genomics"@en ;
    skos:narrower ex:VariantDatabases, ex:ExpressionDatabases ;
    skos:inScheme ex:CategoryScheme .
```

## Alignment with Taxonomy

Maps directly to [taxonomy/](../../taxonomy/) classification:

| SKOS Concept | Taxonomy Category |
|--------------|-------------------|
| ex:Genetics | 01.genetics.genomics |
| ex:Compounds | 02.compounds.molecules |
| ex:Diseases | 03.diseases.phenotypes |
| ex:Pathways | 04.pathways.networks |
