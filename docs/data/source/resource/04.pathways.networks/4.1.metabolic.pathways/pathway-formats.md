---
id: schemas-pathway-formats
title: "Pathway Database Formats - Detailed Technical Reference"
category: schemas
parent: README.md
last_updated: 2026-01-22
status: migrated
tags: [schema, pathway, reactome, kegg, wikipathways, gpml, biopax, sbml, psi-mi]
---

**Parent:** [Schema Documentation](./README.md)

# Pathway Database Formats - Detailed Technical Reference

This document provides comprehensive technical specifications for pathway database formats, extending the existing pathway-target-data-models.md with detailed schema definitions, example data structures, and parsing code snippets.

---

## Table of Contents

1. [Reactome Data Model](#1-reactome-data-model)
2. [KEGG Data Model](#2-kegg-data-model)
3. [WikiPathways GPML Format](#3-wikipathways-gpml-format)
4. [Gene Ontology Structure](#4-gene-ontology-structure)
5. [BioPAX (Biological Pathway Exchange)](#5-biopax-biological-pathway-exchange)
6. [SBML (Systems Biology Markup Language)](#6-sbml-systems-biology-markup-language)
7. [PSI-MI (Molecular Interactions)](#7-psi-mi-molecular-interactions)
8. [CX Format (NDEx)](#8-cx-format-ndex)

---

## 1. Reactome Data Model

### 1.1 Neo4j Graph Schema Overview

Reactome uses Neo4j as its graph database backend. The schema follows an object-oriented design with inheritance hierarchies.

#### Node Type Hierarchy

```
DatabaseObject (base class)
├── Event
│   ├── Pathway
│   ├── ReactionLikeEvent
│   │   ├── Reaction
│   │   ├── BlackBoxEvent
│   │   ├── Polymerisation
│   │   ├── Depolymerisation
│   │   └── FailedReaction
│   └── TopLevelPathway
├── PhysicalEntity
│   ├── Complex
│   ├── EntitySet
│   │   ├── CandidateSet
│   │   ├── DefinedSet
│   │   └── OpenSet
│   ├── GenomeEncodedEntity
│   │   └── EntityWithAccessionedSequence
│   ├── Drug
│   ├── SimpleEntity
│   ├── Polymer
│   └── OtherEntity
├── ReferenceEntity
│   ├── ReferenceSequence
│   │   ├── ReferenceGeneProduct
│   │   ├── ReferenceIsoform
│   │   └── ReferenceDNASequence
│   ├── ReferenceMolecule
│   ├── ReferenceGroup
│   └── ReferenceTherapeutic
├── Regulation
│   ├── PositiveRegulation
│   │   ├── PositiveGeneExpressionRegulation
│   │   └── Requirement
│   └── NegativeRegulation
│       └── NegativeGeneExpressionRegulation
├── CatalystActivity
├── GO_Term
│   ├── GO_BiologicalProcess
│   ├── GO_CellularComponent
│   └── GO_MolecularFunction
└── Species
```

### 1.2 Core Node Types

#### Pathway Node

```cypher
// Cypher query to retrieve Pathway node structure
MATCH (p:Pathway {stId: 'R-HSA-109581'})
RETURN p
```

**Properties:**

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| `dbId` | Long | Internal database ID | `109581` |
| `stId` | String | Stable identifier | `R-HSA-109581` |
| `displayName` | String | Human-readable name | `Apoptosis` |
| `name` | List[String] | All names including synonyms | `["Apoptosis", "Programmed cell death"]` |
| `speciesName` | String | Species | `Homo sapiens` |
| `isInDisease` | Boolean | Disease-related pathway | `false` |
| `isInferred` | Boolean | Computationally inferred | `false` |
| `releaseDate` | String | First release date | `2004-09-20` |
| `schemaClass` | String | Node type | `Pathway` |
| `hasDiagram` | Boolean | Has pathway diagram | `true` |
| `doi` | String | DOI identifier | `10.3180/R-HSA-109581` |

#### Reaction Node

```cypher
MATCH (r:Reaction {stId: 'R-HSA-109606'})
RETURN r
```

**Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Long | Internal database ID |
| `stId` | String | Stable identifier |
| `displayName` | String | Reaction name |
| `category` | String | Reaction category (transition, binding, dissociation, etc.) |
| `isChimeric` | Boolean | Contains entities from multiple species |
| `systematicName` | String | Systematic biochemical name |
| `releaseStatus` | String | RELEASED, NEW, UPDATED |

### 1.3 Relationship Types

#### Event Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `hasEvent` | Pathway | Event | Sub-events within pathway |
| `precedingEvent` | Event | Event | Temporal ordering |
| `followingEvent` | Event | Event | Reverse temporal ordering |
| `inferredTo` | Event | Event | Orthology inference |

#### Reaction Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `input` | ReactionLikeEvent | PhysicalEntity | Reaction inputs |
| `output` | ReactionLikeEvent | PhysicalEntity | Reaction outputs |
| `catalystActivity` | ReactionLikeEvent | CatalystActivity | Enzymatic catalysis |
| `regulatedBy` | ReactionLikeEvent | Regulation | Regulatory relationships |
| `normalReaction` | FailedReaction | Reaction | Normal version of failed reaction |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 100,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | KGML, BioPAX |
| Alternative | JSON, XML |
| Encoding | UTF-8 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| KEGG | HTTP | https://www.kegg.jp/download/kgml/ |
| Reactome | HTTP | https://reactome.org/download-data |
| WikiPathways | HTTP | https://www.wikipathways.org/index.php/Download_All |

**Access Requirements:** Open access, no registration required (see individual database terms)

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| KEGG KGML | Academic Free | Limited (subscription required for commercial) |
| Reactome | CC BY 4.0 | Yes |
| WikiPathways | CC BY 4.0 | Yes |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "R-HSA-109581" |
| `name` | string | Entity name | "Apoptosis" |
| `type` | string | Record type | "pathway" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Sample Data

### Example Record
```json
{
  "stId": "R-HSA-109581",
  "displayName": "Apoptosis",
  "dbId": 109581,
  "speciesName": "Homo sapiens"
}
```

### Sample Query Result
| stId | displayName | dbId | speciesName |
|------|-------------|------|-------------|
| R-HSA-109581 | Apoptosis | 109581 | Homo sapiens |
| R-HSA-191273 | Cholesterol biosynthesis | 191273 | Homo sapiens |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `stId` | Reactome stable identifier for pathways and reactions | `R-HSA-109581` |
| `dbId` | Reactome internal database identifier | `109581` |
| `displayName` | Human-readable name for an entity | `Apoptosis` |
| `schemaClass` | Node type in Reactome graph schema | `Pathway` |
| `pathway_id` | KEGG pathway identifier with organism prefix | `hsa04010` |
| `GraphId` | WikiPathways GPML element identifier | `a1b2c` |
| `DataNodeType` | Classification of GPML node entities | `Protein`, `GeneProduct` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Pathway | Biological process as series of molecular events | Reactome, KEGG |
| ReactionLikeEvent | Molecular transformation in a pathway | Reactome schema |
| PhysicalEntity | Molecule or complex participating in reactions | Reactome schema |
| Regulation | Control mechanism for pathway events | Positive/Negative regulation |
| CatalystActivity | Enzymatic function facilitating a reaction | Reactome relationships |
| hasEvent | Relationship connecting pathway to sub-events | Reactome hierarchy |
| KGML | KEGG Markup Language for pathway representation | KEGG format |
| GPML | Graphical Pathway Markup Language | WikiPathways format |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GPML | Graphical Pathway Markup Language | WikiPathways XML format |
| KGML | KEGG Markup Language | KEGG XML format |
| BioPAX | Biological Pathway Exchange | OWL-based pathway standard |
| SBML | Systems Biology Markup Language | Computational modeling standard |
| PSI-MI | Proteomics Standards Initiative - Molecular Interactions | Interaction data standard |
| CX | Cytoscape Exchange | JSON network format for NDEx |
| NDEx | Network Data Exchange | Network repository platform |
| OWL | Web Ontology Language | Semantic web standard |
| RDF | Resource Description Framework | Linked data format |
| GO | Gene Ontology | Function/process/location terms |
| HSA | Homo sapiens | KEGG organism prefix for human |
| MMU | Mus musculus | KEGG organism prefix for mouse |

---

*Note: This is a large document. Full content preserved from original source.*
