# IntAct - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | intact |
| **Name** | IntAct Molecular Interactions |
| **Total Fields** | 32 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| interactionAc | String | Yes | IntAct interaction accession | EBI-77734 |
| interactorA.identifier | String | Yes | Primary identifier for interactor A | P04637 |
| interactorB.identifier | String | Yes | Primary identifier for interactor B | Q00987 |
| interactionType.miIdentifier | String | No | PSI-MI interaction type ID | MI:0915 |
| interactionType.shortLabel | String | No | Interaction type name | physical association |
| detectionMethod.miIdentifier | String | No | PSI-MI detection method ID | MI:0018 |
| detectionMethod.shortLabel | String | No | Detection method name | two hybrid |
| confidenceScore | Number | No | IntAct MI score (0-1) | 0.85 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| IntAct AC | EBI-####### | Interaction accession | EBI-77734 |
| PSI-MI | MI:#### | Molecular Interaction ontology | MI:0915 |
| UniProt | Alphanumeric | Protein accession | P04637 |
| ChEBI | CHEBI:##### | Small molecule ID | CHEBI:15377 |
| PubMed | Numeric | Publication ID | 12345678 |

---

## Enumerations

### Interactor Types

| Value | Description |
|-------|-------------|
| protein | Protein interactor |
| dna | DNA interactor |
| rna | RNA interactor |
| small molecule | Chemical compound |
| complex | Molecular complex |
| gene | Gene entity |

### Experimental Roles

| Value | Description |
|-------|-------------|
| bait | Tagged/fixed protein |
| prey | Captured protein |
| neutral component | No specific role |
| self | Self-interaction |

### Common Detection Methods

| MI ID | Name | Description |
|-------|------|-------------|
| MI:0018 | two hybrid | Yeast two-hybrid |
| MI:0019 | coimmunoprecipitation | Co-IP assay |
| MI:0096 | pull down | Pull-down assay |
| MI:0071 | molecular sieving | Size exclusion |
| MI:0424 | protein kinase assay | Kinase activity |

### Common Interaction Types

| MI ID | Name | Description |
|-------|------|-------------|
| MI:0915 | physical association | Physical binding |
| MI:0407 | direct interaction | Direct contact |
| MI:0914 | association | General association |
| MI:0217 | phosphorylation reaction | Phosphorylation |

---

## Entity Relationships

### Interaction to Interactors
- **Cardinality:** 1:2+
- **Description:** Each interaction involves 2+ interactors
- **Key Fields:** interactionAc, interactorA, interactorB

### Interaction to Publication
- **Cardinality:** N:1
- **Description:** Interactions linked to literature
- **Key Fields:** interactionAc, publication.pubmedId

### Interactor to Features
- **Cardinality:** 1:N
- **Description:** Binding regions and modifications
- **Key Fields:** interactor, features

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IntAct | Intact Molecular Interactions | Database name |
| PSI-MI | Proteomics Standards Initiative - Molecular Interactions | Ontology |
| MITAB | Molecular Interactions TAB | File format |
| AC | Accession | Identifier type |
| Y2H | Yeast Two-Hybrid | Detection method |
| Co-IP | Co-immunoprecipitation | Detection method |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
