---
id: gene.function.ontology
title: "Gene Function and Ontology"
type: subcategory
parent: ../README.md
last_updated: 2026-01-23
status: active
tags: [gene-ontology, go, function, annotations, gene-sets]
---

# Gene Function and Ontology

**Parent:** [Pathways & Networks](../README.md)

## Overview

Gene function databases provide standardized functional annotations and curated gene sets for enrichment analysis. These resources enable systematic characterization of gene functions and identification of enriched biological processes.

Key resources include Gene Ontology (GO) for functional annotations and MSigDB (Molecular Signatures Database) for curated gene sets. Together they support gene set enrichment analysis and functional interpretation.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [Gene Ontology](./gene.ontology/README.md) | 1 | Functional annotation ontology |
| [MSigDB](./msigdb/README.md) | 1 | Molecular Signatures Database |

## Integration Notes

GO provides three ontologies: Molecular Function, Biological Process, and Cellular Component. MSigDB offers hallmark, positional, and curated gene sets. Use GO for individual gene annotation; use MSigDB for enrichment analysis. Both support GSEA and similar methods.
