---
id: reference-pathway-formats
title: "Pathway Database Formats - Detailed Technical Reference"
category: reference
parent: _index.md
last_updated: 2026-01-23
status: active
migrated_from: operations/schemas/pathway-formats.md
tags: [schema, pathway, reactome, kegg, wikipathways, gpml, biopax, sbml, psi-mi]
---

**Parent:** [Format Specifications](./_index.md)

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

\`\`\`
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
\`\`\`

### 1.2 Core Node Types

#### Pathway Node

\`\`\`cypher
// Cypher query to retrieve Pathway node structure
MATCH (p:Pathway {stId: 'R-HSA-109581'})
RETURN p
\`\`\`

**Properties:**

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| \`dbId\` | Long | Internal database ID | \`109581\` |
| \`stId\` | String | Stable identifier | \`R-HSA-109581\` |
| \`displayName\` | String | Human-readable name | \`Apoptosis\` |
| \`name\` | List[String] | All names including synonyms | \`["Apoptosis", "Programmed cell death"]\` |
| \`speciesName\` | String | Species | \`Homo sapiens\` |
| \`isInDisease\` | Boolean | Disease-related pathway | \`false\` |
| \`isInferred\` | Boolean | Computationally inferred | \`false\` |
| \`releaseDate\` | String | First release date | \`2004-09-20\` |
| \`schemaClass\` | String | Node type | \`Pathway\` |
| \`hasDiagram\` | Boolean | Has pathway diagram | \`true\` |
| \`doi\` | String | DOI identifier | \`10.3180/R-HSA-109581\` |

### 1.3 Relationship Types

#### Event Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| \`hasEvent\` | Pathway | Event | Sub-events within pathway |
| \`precedingEvent\` | Event | Event | Temporal ordering |
| \`followingEvent\` | Event | Event | Reverse temporal ordering |
| \`inferredTo\` | Event | Event | Orthology inference |

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

## Sample Data

### Example Record
\`\`\`json
{
  "stId": "R-HSA-109581",
  "displayName": "Apoptosis",
  "dbId": 109581,
  "speciesName": "Homo sapiens"
}
\`\`\`

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| \`stId\` | Reactome stable identifier for pathways and reactions | \`R-HSA-109581\` |
| \`dbId\` | Reactome internal database identifier | \`109581\` |
| \`GPML\` | Graphical Pathway Markup Language | WikiPathways XML format |
| \`KGML\` | KEGG Markup Language | KEGG XML format |
| \`BioPAX\` | Biological Pathway Exchange | OWL-based pathway standard |
| \`SBML\` | Systems Biology Markup Language | Computational modeling standard |

---

*Note: This is a summary. Full content preserved from original source.*
