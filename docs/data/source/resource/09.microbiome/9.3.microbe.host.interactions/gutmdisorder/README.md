---
id: gutmdisorder
title: "gutMDisorder"
type: data-source
category: microbiome
subcategory: microbe.host.interactions
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [microbiome, disease, dysbiosis, associations, curated]
---

# gutMDisorder

**Category:** [Microbiome](../../../_index.md) > [Microbe-Host Interactions](../_index.md)

## Overview

gutMDisorder is a manually curated database of gut microbiota-disease associations. It compiles evidence from published studies linking specific microbial taxa to human diseases, noting whether species are increased, decreased, or altered in disease states.

The database covers a wide range of conditions including metabolic disorders, inflammatory diseases, neurological conditions, and cancers. Each association is annotated with the original publication, sample information, and evidence strength.

gutMDisorder is essential for understanding dysbiosis patterns across diseases and identifying microbial biomarkers.

## Key Statistics

| Metric | Value |
|--------|-------|
| Microbe-Disease Pairs | 10,000+ |
| Diseases | 200+ |
| Microbial Taxa | 2,000+ |
| Publications | 2,000+ |
| Last Update | 2023 |

## Primary Use Cases

1. Disease-microbiome association lookup
2. Dysbiosis pattern analysis
3. Microbial biomarker identification
4. Systematic review support
5. Research hypothesis generation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Disease Name | Text | Type 2 Diabetes |
| ICD Code | ICD-10 | E11 |
| NCBI Taxon ID | Numeric | 816 |
| PMID | Numeric | 12345678 |
| Association Type | Enum | Increased, Decreased |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://bio-annotation.cn/gutMDisorder | Search and browse |
| Downloads | Available | TSV format |
| API | N/A | No public API |

## License

| Aspect | Value |
|--------|-------|
| License | Free for academic use |
| Commercial Use | Contact maintainers |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [GMrepo](../../9.1.gut.microbiome/gmrepo/_index.md) - Sample data
- [VMH](../vmh/_index.md) - Metabolic models
