---
id: ontology
title: "Ontology"
type: ontology
parent: ../_index.md
last_updated: 2026-01-25
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
| [skos/datasource-taxonomy.ttl](./skos/datasource-taxonomy.ttl) | SKOS taxonomy with 6 concept schemes |
| [owl/datasource-ontology.ttl](./owl/datasource-ontology.ttl) | OWL 2 classes and properties for data sources |
| [shacl/datasource-shapes.ttl](./shacl/datasource-shapes.ttl) | SHACL validation shapes |

## Namespaces

| Prefix | URI | Description |
|--------|-----|-------------|
| `tax:` | `https://gene.example.org/taxonomy/` | SKOS taxonomy concepts |
| `ds:` | `https://gene.example.org/ontology/datasource#` | OWL classes and properties |
| `shapes:` | `https://gene.example.org/shapes/datasource#` | SHACL validation shapes |

## Data Source Ontology

### SKOS Taxonomy (6 Concept Schemes)

- **DataSourceTaxonomy** - 9 categories, ~45 subcategories (3-level hierarchy)
- **TierScheme** - Tier1, Tier2, Tier3 (quality classification)
- **StatusScheme** - Draft, Active, Deprecated (lifecycle)
- **AccessMethodScheme** - REST API, SPARQL, FTP, S3, etc.
- **FormatScheme** - JSON, XML, CSV, VCF, FASTA, etc.
- **LicenseTypeScheme** - Public Domain, CC-BY, Academic Only, etc.

### OWL Classes

- **ds:Category** - Top-level domain (e.g., Genetics and Genomics)
- **ds:Subcategory** - Second-level classification
- **ds:DataSource** - Individual data source (e.g., ClinVar)
- **ds:Version** - Release with date, size, checksum
- **ds:License** - Usage terms and restrictions
- **ds:AccessMethod** - API or download endpoint
- **ds:Schema** - Data structure definition
- **ds:Field** - Schema field with type and constraints
- **ds:CrossReference** - Link between sources

### SHACL Shapes

- **DataSourceShape** - Required: id, title, category, subcategory, status
- **VersionShape** - Required: versionNumber, releaseDate
- **LicenseShape** - Required: licenseName

## Semantic Web Standards

- **OWL 2** (Web Ontology Language) - Formal class definitions and reasoning
- **SHACL** (Shapes Constraint Language) - Data validation and quality
- **SKOS** (Simple Knowledge Organization System) - Taxonomies and thesauri

## Related

- [taxonomy/](../taxonomy/) - Classification system documentation
- [reference/architecture/](../reference/architecture/) - Schema designs
- [templates/](../templates/) - Data source documentation templates
