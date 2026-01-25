---
id: allen.brain.atlas
title: "Allen Brain Atlas"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: mental.health.neurological
tags:
  - brain
  - neuroanatomy
  - gene-expression
  - transcriptomics
  - neuroscience
---

# Allen Brain Atlas

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

## Limitations

- Human brain data from limited number of donors (6)
- Postmortem tissue may not reflect living brain biology
- Microarray data lower resolution than RNA-seq
- Some brain regions under-sampled

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [PGC](../pgc/README.md) - Psychiatric genetics
