---
id: chemical.ontology.classification
title: "Chemical Ontology and Classification"
type: subcategory
parent: ../README.md
last_updated: 2026-01-23
status: active
tags: [ontology, classification, chebi, pubchem, taxonomy]
---

# Chemical Ontology and Classification

**Parent:** [Compounds & Molecules](../README.md)

## Overview

Chemical ontology and classification databases provide standardized chemical identifiers, hierarchical classifications, and structural annotations. These resources enable consistent compound identification and semantic integration across databases.

Key resources include ChEBI (chemical ontology), PubChem (comprehensive chemical database), ClassyFire (automated classification), and NPClassifier (natural product classification). Together they provide the foundation for chemical data integration.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [ChEBI](./chebi/README.md) | 1 | Chemical Entities of Biological Interest |
| [PubChem](./pubchem/README.md) | 1 | NIH chemical database |
| [ClassyFire](./classyfire/README.md) | 2 | Automated chemical classification |
| [NPClassifier](./npclassifier/README.md) | 2 | Natural product classifier |

## Integration Notes

PubChem CIDs provide universal compound identifiers. ChEBI offers semantic ontology relationships. ClassyFire assigns hierarchical chemical classes automatically from structures. NPClassifier specializes in biosynthetic pathway-based classification. Use InChIKey for structure-based deduplication.
