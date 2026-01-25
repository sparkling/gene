---
id: panelapp
title: "PanelApp - Genomics England Gene Panels"
type: source
parent: ../README.md
tier: 2
status: active
category: diseases.phenotypes
subcategory: rare.orphan.diseases
tags:
  - gene-panels
  - clinical-genomics
  - diagnostic-genes
  - rare-diseases
  - nhs
---

# PanelApp - Genomics England Gene Panels

## Overview

PanelApp is a crowdsourced knowledge base for gene-disease relationships developed by Genomics England for use in the NHS National Genomic Medicine Service. The platform provides curated gene panels for clinical diagnostic use, with expert-reviewed evidence supporting gene-disease associations.

Each gene panel in PanelApp represents a clinical indication (disease or group of diseases) and contains genes categorized by their diagnostic evidence level: Green (high evidence for diagnostic use), Amber (moderate evidence), and Red (low/no evidence). The evidence review process involves clinical scientists, geneticists, and disease specialists contributing to a consensus.

PanelApp panels are used in the 100,000 Genomes Project and NHS clinical genomics services to determine which genes should be analyzed for specific clinical presentations. The platform is also used internationally by clinical laboratories for panel design and variant interpretation.

## Key Statistics

| Metric | Value |
|--------|-------|
| Gene Panels | 350+ |
| Green Genes | 4,500+ |
| Expert Reviews | 50,000+ |
| Countries Using | 20+ |
| Updates | Continuous |

## Primary Use Cases

1. Clinical gene panel design for diagnostics
2. Variant interpretation for rare diseases
3. Evidence assessment for gene-disease links
4. Test development for clinical laboratories
5. Research into gene-disease associations

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Panel ID | `[0-9]+` | 245 |
| Panel Name | Text | Intellectual disability |
| Gene Symbol | HGNC symbol | MECP2 |
| HGNC ID | `HGNC:[0-9]+` | HGNC:6990 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| PanelApp Web | https://panelapp.genomicsengland.co.uk/ | UK version |
| PanelApp Australia | https://panelapp.agha.umccr.org/ | Australian mirror |
| REST API | https://panelapp.genomicsengland.co.uk/api/ | Programmatic |
| Downloads | https://panelapp.genomicsengland.co.uk/panels/ | Panel exports |

## Evidence Levels

| Level | Color | Criteria |
|-------|-------|----------|
| 3 | Green | Diagnostic grade, high evidence |
| 2 | Amber | Moderate evidence, emerging |
| 1 | Red | Low evidence, not for diagnosis |
| 0 | Gray | No evidence yet |

## Data Formats

| Format | Endpoint | Notes |
|--------|----------|-------|
| JSON | API responses | Full panel data |
| TSV | Panel exports | Gene lists |
| BED | Gene coordinates | For bioinformatics |

## Limitations

- Panels reflect UK clinical practice; may differ internationally
- Evidence levels are crowd-sourced expert opinion
- Continuous updates mean panel versions may change
- Not all rare diseases have curated panels

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [Orphanet](../orphanet/README.md) - Rare disease database
