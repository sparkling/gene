# Reactome - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | reactome |
| **Name** | Reactome |
| **Total Fields** | 45 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| dbId | Integer | Yes | Internal database identifier | 1430728 |
| stId | String | Yes | Stable identifier (R-{species}-{number}) | R-HSA-1430728 |
| displayName | String | Yes | Human-readable name | Metabolism |
| speciesName | String | Yes | Species name | Homo sapiens |
| schemaClass | String | Yes | Entity type | TopLevelPathway, Reaction |
| isInDisease | Boolean | No | Disease pathway flag | true, false |
| isInferred | Boolean | No | Inferred from orthology | true, false |
| hasDiagram | Boolean | No | Has pathway diagram | true, false |
| hasEHLD | Boolean | No | Has Enhanced High-Level Diagram | true, false |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Reactome stId | R-{species}-{number} | Stable identifier | R-HSA-1430728 |
| UniProt | [A-Z][0-9]{5} | Protein accession | Q9UBM7 |
| ChEBI | CHEBI:##### | Small molecule identifier | CHEBI:30616 |
| GO | GO:####### | Gene Ontology term | GO:0005654 |
| Ensembl | ENSG########### | Gene identifier | ENSG00000113161 |
| PubMed | ######## | Literature reference | 10583946 |

---

## Enumerations

### Schema Classes (Event)

| Value | Description |
|-------|-------------|
| TopLevelPathway | Root-level pathway entries |
| Pathway | Collection of related events |
| Reaction | Direct input-to-output conversion |
| BlackBoxEvent | Incompletely characterized reaction |
| Polymerisation | Polymer assembly |
| Depolymerisation | Polymer breakdown |
| FailedReaction | Known failed reaction |

### Schema Classes (Physical Entity)

| Value | Description |
|-------|-------------|
| EntityWithAccessionedSequence | Protein/nucleic acid with accession |
| SimpleEntity | Small molecule |
| Complex | Multi-component assembly |
| DefinedSet | Specific set of entities |
| CandidateSet | Candidate entities |
| Polymer | Polymeric molecule |
| Drug | Therapeutic compound |

### Species Codes

| Code | Species | Tax ID |
|------|---------|--------|
| HSA | Homo sapiens | 9606 |
| MMU | Mus musculus | 10090 |
| RNO | Rattus norvegicus | 10116 |
| DRE | Danio rerio | 7955 |
| DME | Drosophila melanogaster | 7227 |
| CEL | Caenorhabditis elegans | 6239 |
| SCE | Saccharomyces cerevisiae | 4932 |

---

## Entity Relationships

### Pathway to Event
- **Cardinality:** 1:N
- **Description:** Pathways contain sub-events
- **Key Fields:** hasEvent relationship

### Reaction to Physical Entity
- **Cardinality:** N:M
- **Description:** Reactions have inputs and outputs
- **Key Fields:** input, output relationships

### Physical Entity to Reference Entity
- **Cardinality:** N:1
- **Description:** Entities reference canonical forms
- **Key Fields:** referenceEntity relationship

### Entity to Compartment
- **Cardinality:** N:1
- **Description:** Entities localized to compartments
- **Key Fields:** compartment relationship

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| stId | Stable Identifier | R-HSA-###### format |
| HSA | Homo sapiens | Human species code |
| GO | Gene Ontology | Compartment/function annotations |
| BioPAX | Biological Pathway Exchange | OWL-based pathway format |
| SBML | Systems Biology Markup Language | Simulation format |
| SBGN | Systems Biology Graphical Notation | Visualization standard |
| EHLD | Enhanced High-Level Diagram | Reactome illustration |
| ChEBI | Chemical Entities of Biological Interest | Small molecule database |
| PSI-MITAB | Proteomics Standards Initiative MITAB | Interaction format |
| OWL | Web Ontology Language | RDF-based ontology format |
| Neo4j | Graph Database | Reactome storage platform |

---

## See Also

- [schema.md](./schema.md) - Full schema documentation
- [sample.json](./sample.json) - Example records
