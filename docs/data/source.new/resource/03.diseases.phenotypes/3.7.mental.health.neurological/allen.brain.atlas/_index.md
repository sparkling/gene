---
id: allen.brain.atlas
title: "Allen Brain Atlas"
type: data-source
category: diseases
subcategory: mental.health.neurological
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [brain, neuroanatomy, gene-expression, transcriptomics, neuroscience]
---

# Allen Brain Atlas

**Category:** [Diseases & Phenotypes](../../_index.md) > [Mental Health & Neurological](../_index.md)

## Overview

The Allen Brain Atlas is a comprehensive collection of publicly available online resources integrating gene expression and neuroanatomical data across the mammalian brain. Developed by the Allen Institute for Brain Science, these atlases provide foundational data for understanding brain structure, function, and dysfunction in neurological and psychiatric disorders.

The Allen Human Brain Atlas contains genome-wide microarray-based gene expression data from ~900 anatomically-defined brain structures across multiple adult human brains. Complementary resources include the BrainSpan developmental transcriptome atlas, the Allen Cell Types Database with single-cell RNA-seq data, and connectivity atlases mapping neural circuits.

These resources enable researchers to identify brain region-specific gene expression patterns relevant to neurological diseases, map candidate disease genes to specific cell types and circuits, and understand how genetic variants associated with psychiatric disorders affect brain biology.

## Key Statistics

| Metric | Value |
|--------|-------|
| Human Brain Samples | 6 donors |
| Brain Structures | ~900 |
| Genes Profiled | 20,000+ |
| Single-Cell Profiles | 1M+ |
| Mouse Brain Sections | 85,000+ |

## Primary Use Cases

1. Neurological disease gene expression mapping
2. Cell type-specific gene expression analysis
3. Brain region prioritization for disease genes
4. Neurodevelopmental gene expression patterns
5. Cross-species brain evolution studies

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Structure ID | `[0-9]+` | 4020 (hippocampus) |
| Structure Acronym | `[A-Za-z]+` | HIP |
| Gene Symbol | HGNC | BDNF |
| Donor ID | `H[0-9]+` | H0351.2001 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Allen Portal | https://portal.brain-map.org/ | Main portal |
| Human Brain Atlas | https://human.brain-map.org/ | Human data |
| BrainSpan | https://www.brainspan.org/ | Development |
| API | https://api.brain-map.org/ | Programmatic |

## Data Resources

| Atlas | Species | Data Types |
|-------|---------|------------|
| Human Brain Atlas | Human | Microarray, ISH |
| BrainSpan | Human | RNA-seq (development) |
| Mouse Brain Atlas | Mouse | ISH, connectivity |
| Cell Types Database | Human/Mouse | Single-cell RNA-seq |
| Mouse Connectivity | Mouse | Projection mapping |

## Data Formats

| Format | Description | Access |
|--------|-------------|--------|
| CSV | Expression matrices | API download |
| JSON | API responses | REST API |
| NIfTI | Brain images | Download portal |
| HDF5 | Large datasets | Bulk download |

## License

| Aspect | Value |
|--------|-------|
| License | Allen Institute Terms of Use |
| Commercial Use | Yes with attribution |
| Attribution | Required |
| Citation | Allen Institute for Brain Science |

## See Also

- [Schema Documentation](./schema.md)
- [GTEx](../../../01.genetics.genomics/1.5.expression.regulation/gtex/_index.md) - Multi-tissue expression
- [PGC](../pgc/_index.md) - Psychiatric genetics
