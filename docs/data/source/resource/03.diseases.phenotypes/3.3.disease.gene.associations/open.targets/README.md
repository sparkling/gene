---
id: open.targets
title: "Open Targets Platform"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: disease.gene.associations
tags:
  - drug-targets
  - evidence
  - associations
  - graphql
  - target-validation
---

# Open Targets Platform

## Overview

The Open Targets Platform is a comprehensive resource for systematic identification and prioritization of drug targets. It integrates evidence from genetics, genomics, transcriptomics, drugs, animal models, and scientific literature to build target-disease associations with transparent evidence scoring.

Developed through a public-private partnership between the Wellcome Sanger Institute, EMBL-EBI, GSK, and other pharmaceutical partners, Open Targets provides both a web platform and programmatic access via GraphQL API. The platform enables researchers to evaluate the therapeutic potential of targets with evidence from over 20 data sources.

Each target-disease association includes detailed evidence breakdown by data type (genetic associations, somatic mutations, known drugs, pathways, RNA expression, text mining, animal models). The platform also provides tractability assessments, safety information, and known drug mechanisms of action, supporting the full target validation workflow.

## Key Statistics

| Metric | Value |
|--------|-------|
| Target-Disease Associations | 14M+ |
| Targets (Genes) | 62,000+ |
| Diseases | 20,000+ |
| Evidence Strings | 17M+ |
| Data Sources | 20+ |

## Primary Use Cases

1. Drug target identification and prioritization
2. Evidence integration for target validation
3. Known drug and clinical trial exploration
4. Genetic association analysis
5. Safety and tractability assessment

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Ensembl Gene | `ENSG[0-9]{11}` | ENSG00000146648 |
| EFO Disease | `EFO_[0-9]{7}` | EFO_0000685 |
| ChEMBL Drug | `CHEMBL[0-9]+` | CHEMBL941 |
| UniProt | `[A-Z0-9]{6,10}` | P00533 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Platform Web | https://platform.opentargets.org/ | Web interface |
| GraphQL API | https://api.platform.opentargets.org/api/v4/graphql | API access |
| Data Downloads | https://platform.opentargets.org/downloads | Bulk files |
| GitHub | https://github.com/opentargets | Source code |

## Data Formats

| Format | File | Notes |
|--------|------|-------|
| JSON | GraphQL responses | API access |
| Parquet | Bulk downloads | Big data analysis |
| TSV | Selected exports | Spreadsheet-friendly |

## Evidence Types

| Data Type | Sources |
|-----------|---------|
| Genetic Associations | GWAS Catalog, Gene2Phenotype, ClinGen |
| Somatic Mutations | Cancer Gene Census, IntOGen |
| Known Drugs | ChEMBL |
| Pathways | Reactome, SLAPenrich |
| RNA Expression | Expression Atlas |
| Animal Models | IMPC, MGI |
| Text Mining | Europe PMC |

## Limitations

- Evidence scores are aggregate; individual source quality varies
- Tractability assessments are predictions, not experimental validation
- Quarterly releases mean some data lags source databases
- EFO disease ontology may not map to all clinical terminologies

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [DisGeNET](../disgenet/README.md) - Complementary gene-disease resource
