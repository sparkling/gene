---
id: syngo
title: "SynGO - Synaptic Gene Ontologies"
type: source
parent: ../README.md
tier: 2
status: active
category: diseases.phenotypes
subcategory: mental.health.neurological
tags:
  - synapse
  - neuroscience
  - gene-ontology
  - synaptic-function
  - brain-disorders
---

# SynGO - Synaptic Gene Ontologies

## Overview

SynGO (Synaptic Gene Ontologies) is an expert-curated resource for synaptic function gene annotations, providing a knowledge base for synaptic protein localization and function. The resource uses a structured ontology framework to describe synaptic biology at both cellular component (where) and biological process (what function) levels.

Developed by a consortium of synapse biology experts, SynGO addresses the gap in Gene Ontology coverage for synaptic processes. Each annotation is supported by published experimental evidence, with over 3,000 genes annotated across pre-synaptic and post-synaptic compartments and their specific functional roles in synaptic transmission.

SynGO is particularly valuable for interpreting genetic studies of neuropsychiatric and neurodevelopmental disorders, which show strong enrichment for synaptic genes. The resource enables researchers to test whether disease-associated genes converge on specific synaptic processes or compartments.

## Key Statistics

| Metric | Value |
|--------|-------|
| Annotated Genes | 1,200+ |
| Annotations | 5,800+ |
| Supporting Papers | 2,800+ |
| Ontology Terms | 180+ |
| Expert Contributors | 60+ |

## Primary Use Cases

1. Gene set enrichment analysis for synaptic genes
2. Neurological disease gene interpretation
3. Synaptic protein localization lookup
4. Synaptic function annotation
5. Drug target prioritization for brain disorders

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| SynGO Term | `SYNGO:[0-9]+` | SYNGO:synapse |
| Gene Symbol | HGNC | SYN1 |
| UniProt ID | `[A-Z0-9]{6,10}` | P17600 |
| GO Term | `GO:[0-9]{7}` | GO:0045202 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| SynGO Portal | https://syngoportal.org/ | Web interface |
| Annotation Browser | https://syngoportal.org/browse | Query genes |
| Downloads | https://syngoportal.org/download | Bulk data |
| GitHub | https://github.com/syngoportal/SynGO-data | Source files |

## Ontology Structure

| Branch | Description | Terms |
|--------|-------------|-------|
| Cellular Component | Synaptic localization | 90+ |
| Biological Process | Synaptic function | 90+ |
| Pre-synaptic | Axon terminal | 50+ |
| Post-synaptic | Dendritic spine | 60+ |

## Data Formats

| Format | File | Description |
|--------|------|-------------|
| OBO | syngo.obo | Ontology structure |
| GAF | syngo_annotations.gaf | Gene annotations |
| TSV | syngo_genes.tsv | Gene list |
| JSON | API responses | Web service |

## Evidence Codes

| Code | Meaning |
|------|---------|
| ECO:0005593 | Biological assay |
| ECO:0005589 | Immunolocalization |
| ECO:0005644 | Electron microscopy |
| ECO:0006063 | Biochemical assay |

## Limitations

- Focused on synapse biology; non-synaptic brain genes not covered
- Annotation depth varies by gene; well-studied genes better annotated
- Evidence limited to published experimental data
- Mouse data may not fully translate to human

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [Allen Brain Atlas](../allen.brain.atlas/README.md) - Brain expression data
- [PGC](../pgc/README.md) - Psychiatric genetics
