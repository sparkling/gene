---
id: disease.ontologies
title: "Disease Ontologies"
type: subcategory
parent: ../_index.md
last_updated: 2026-01-23
status: active
tags: [ontology, classification, icd, mondo, mesh, terminology]
---

# Disease Ontologies

**Parent:** [Diseases & Phenotypes](../_index.md)

## Overview

Disease ontologies provide standardized vocabularies and hierarchical classifications for diseases and medical conditions. These resources enable consistent disease annotation across databases and support semantic integration of biomedical data.

Key resources include MONDO (unified disease ontology), EFO (Experimental Factor Ontology), ICD (International Classification of Diseases), and MeSH (Medical Subject Headings). Together they provide comprehensive coverage for research and clinical applications.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [MONDO](./mondo/_index.md) | 1 | Unified disease ontology |
| [EFO](./efo/_index.md) | 1 | Experimental Factor Ontology |
| [ICD](./icd/_index.md) | 1 | International Classification of Diseases |
| [MeSH](./mesh/_index.md) | 1 | Medical Subject Headings |

## Integration Notes

MONDO provides mappings to multiple disease vocabularies. EFO is used by GWAS Catalog and Open Targets. ICD codes are standard in clinical settings. MeSH terms link to PubMed literature. Use MONDO for cross-database integration; use ICD for clinical data mapping.
