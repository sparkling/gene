---
title: "Methodology"
parent: ../_index.md
last_updated: 2026-01-23
status: draft
---

# Methodology

This section documents data acquisition methodologies, pipeline designs, coverage analysis approaches, and source evaluation frameworks.

## Contents

| Document | Description |
|----------|-------------|
| [Literature Pipeline](literature-pipeline.md) | ETL pipeline design for research papers from PubMed/PMC |
| [Coverage Analysis](coverage-analysis.md) | Literature coverage gap analysis across databases |
| [Public Sources](public-sources.md) | Guide to public research paper sources (PubMed, PMC, Europe PMC) |
| [Literature Sources](literature-sources.md) | Comprehensive literature and research paper source reference |
| [Curation Framework](curation-framework.md) | Data source evaluation and pruning methodology |

## Overview

The methodology section provides:

- **Pipeline Design**: ETL workflows for data extraction and processing
- **Coverage Analysis**: Assessment of database coverage and gaps
- **Source Documentation**: Detailed guides for accessing public data sources
- **Curation Standards**: Frameworks for evaluating and selecting data sources

## Key Methodologies

### Literature Pipeline

The research papers ETL pipeline processes 1-10M genetics-relevant papers from PubMed/PMC:
- Download, filter, and parse XML records
- Extract entities (genes, SNPs, compounds)
- Generate embeddings for semantic search
- Store in RuVector with graph relationships

### Coverage Analysis

Systematic evaluation of database coverage to identify:
- High-coverage areas (90-95% for SNP/GWAS studies)
- Coverage gaps (alternative medicine, non-English journals, preprints)
- Multi-source strategies for gap-filling

### Public Sources

Documentation of freely accessible research literature:
- PubMed (39M+ citations)
- PubMed Central (10M+ full-text articles)
- Europe PMC (43M+ abstracts with text mining)
- OpenAlex (250M+ works)

## Navigation

- [Parent: Operations](../_index.md)
- [Downloads](../downloads/_index.md)
- [Integration](../integration/_index.md)
- [Governance](../governance/_index.md)
- [Planning](../planning/_index.md)
- [Schemas](../schemas/_index.md)
