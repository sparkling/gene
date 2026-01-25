# Pathway Commons - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pathwaycommons |
| **Name** | Pathway Commons |
| **Total Fields** | 28 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| uri | String | Yes | Unique identifier URI for the entity | http://identifiers.org/reactome/R-HSA-109581 |
| displayName | String | Yes | Human-readable name | Apoptosis |
| name | Array | No | Alternative names for the entity | ["Programmed Cell Death"] |
| biopaxClass | Enum | No | BioPAX Level 3 entity class | Pathway, BiochemicalReaction, Protein |
| organism.displayName | String | No | Species name | Homo sapiens |
| organism.taxonId | String | No | NCBI Taxonomy ID | 9606 |
| dataSource | Array | No | Contributing database(s) | ["Reactome", "KEGG"] |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| URI | http://identifiers.org/* | Standard identifier URI | http://identifiers.org/reactome/R-HSA-109581 |
| Taxon ID | Numeric string | NCBI Taxonomy identifier | 9606 |
| GO ID | GO:####### | Gene Ontology term | GO:0005737 |
| UniProt | Alphanumeric | Protein accession | P04637 |

---

## Enumerations

### BioPAX Entity Classes

| Value | Description |
|-------|-------------|
| Pathway | Biological pathway |
| BiochemicalReaction | Chemical transformation |
| Control | Regulatory interaction |
| Catalysis | Enzymatic catalysis |
| Protein | Protein entity |
| SmallMolecule | Chemical compound |
| Complex | Molecular complex |
| Dna | DNA entity |
| Rna | RNA entity |
| PhysicalEntity | Generic physical entity |

### Conversion Direction

| Value | Description |
|-------|-------------|
| LEFT-TO-RIGHT | Forward reaction |
| RIGHT-TO-LEFT | Reverse reaction |
| REVERSIBLE | Bidirectional reaction |

### Control Type

| Value | Description |
|-------|-------------|
| ACTIVATION | Positive regulation |
| INHIBITION | Negative regulation |

### SIF Interaction Types

| Value | Description |
|-------|-------------|
| INTERACTS_WITH | Generic interaction |
| IN_COMPLEX_WITH | Complex membership |
| CONTROLS-STATE-CHANGE-OF | State modification |
| CONTROLS-TRANSPORT-OF | Localization control |
| CONTROLS-PHOSPHORYLATION-OF | Phosphorylation regulation |
| CONTROLS-EXPRESSION-OF | Expression regulation |
| CATALYSIS-PRECEDES | Sequential catalysis |
| NEIGHBOR_OF | Network neighbor |
| CHEMICAL-AFFECTS | Chemical effect |

---

## Entity Relationships

### Pathway to Component
- **Cardinality:** 1:N
- **Description:** Pathway contains reactions and controls
- **Key Fields:** uri, pathwayComponent

### Reaction to Participant
- **Cardinality:** N:M
- **Description:** Reactions have substrates (left) and products (right)
- **Key Fields:** left, right

### Control to Target
- **Cardinality:** N:1
- **Description:** Controller entities regulate controlled processes
- **Key Fields:** controller, controlled

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| BioPAX | Biological Pathway Exchange | Standard pathway format |
| SIF | Simple Interaction Format | Binary interaction format |
| PC | Pathway Commons | Database abbreviation |
| URI | Uniform Resource Identifier | Entity identifier |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
