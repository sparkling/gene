---
id: ontology-owl
title: "OWL Ontologies"
type: ontology
parent: ../_index.md
last_updated: 2026-01-25
status: active
tags: [owl, ontology, semantic-web, reasoning]
---

# OWL Ontologies

**Parent:** [Ontology](../_index.md)

OWL 2 (Web Ontology Language) definitions for formal class hierarchies, object properties, and reasoning axioms.

## Purpose

- Define formal class hierarchies for data sources and metadata
- Specify object and data properties with domains/ranges
- Enable automated reasoning and inference
- Support interoperability with SKOS taxonomy

## Files

| File | Description |
|------|-------------|
| [datasource-ontology.ttl](./datasource-ontology.ttl) | OWL 2 ontology for data source metadata |

## Namespace

```turtle
@prefix ds: <https://gene.example.org/ontology/datasource#> .
```

## Classes

### Core Hierarchy

| Class | Description |
|-------|-------------|
| `ds:Category` | Top-level domain (e.g., Genetics and Genomics) |
| `ds:Subcategory` | Second-level classification |
| `ds:DataSource` | Individual data source (e.g., ClinVar, UniProt) |

### Metadata

| Class | Description |
|-------|-------------|
| `ds:Version` | Release version with date, size, checksum |
| `ds:License` | Usage terms and restrictions |
| `ds:AccessMethod` | API, download, or query endpoint |
| `ds:Format` | Data format (JSON, XML, VCF, etc.) |
| `ds:Maintainer` | Organization maintaining the source |

### Schema

| Class | Description |
|-------|-------------|
| `ds:Schema` | Data structure definition |
| `ds:Field` | Individual field with type and constraints |
| `ds:FieldGroup` | Logical grouping of fields |

### Cross-References

| Class | Description |
|-------|-------------|
| `ds:CrossReference` | Link between sources |
| `ds:Identifier` | Identifier type (e.g., VCV, UniProt accession) |
| `ds:IdentifierMapping` | Mapping between identifier types |

## Key Properties

| Property | Type | Domain | Range |
|----------|------|--------|-------|
| `ds:belongsToCategory` | Object | DataSource/Subcategory | Category |
| `ds:hasVersion` | Object | DataSource | Version |
| `ds:hasLicense` | Object | DataSource | License |
| `ds:tier` | Object | DataSource | skos:Concept |
| `ds:status` | Object | HierarchyLevel | skos:Concept |

## Defined Classes (Inference)

| Class | Definition |
|-------|------------|
| `ds:Tier1Source` | DataSource with tier = tax:Tier1 |
| `ds:ActiveSource` | DataSource with status = tax:Active |
| `ds:OpenAccessSource` | DataSource with commercialUseAllowed = true |
| `ds:ApiAccessibleSource` | DataSource with REST API access |

## External Ontologies

| Ontology | URI | Use |
|----------|-----|-----|
| SKOS | http://www.w3.org/2004/02/skos/core# | Taxonomy concepts |
| FOAF | http://xmlns.com/foaf/0.1/ | Organization class |
| DCT | http://purl.org/dc/terms/ | Dublin Core metadata |
