---
id: schema-allen-brain-atlas
title: "Allen Brain Atlas Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, brain, neuroanatomy, gene-expression, transcriptomics, neuroscience]
---

# Allen Brain Atlas Schema Documentation

**Document ID:** SCHEMA-ALLEN-BRAIN-ATLAS

---

## TL;DR

The Allen Brain Atlas is a comprehensive collection of publicly available resources integrating gene expression and neuroanatomical data across the mammalian brain. It contains genome-wide expression data from ~900 anatomically-defined brain structures with complementary single-cell and developmental resources.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Human Brain Samples | 6 donors | Allen Institute |
| Brain Structures | ~900 | Allen Institute |
| Genes Profiled | 20,000+ | Allen Institute |
| Single-Cell Profiles | 1M+ | Allen Institute |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | CSV, JSON, NIfTI, HDF5 |
| API | Yes (REST API) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Structure ID | `[0-9]+` | 4020 (hippocampus) |
| Structure Acronym | `[A-Za-z]+` | HIP |
| Gene Symbol | HGNC | BDNF |
| Donor ID | `H[0-9]+` | H0351.2001 |

---

## Data Resources

| Atlas | Species | Data Types |
|-------|---------|------------|
| Human Brain Atlas | Human | Microarray, ISH |
| BrainSpan | Human | RNA-seq (development) |
| Mouse Brain Atlas | Mouse | ISH, connectivity |
| Cell Types Database | Human/Mouse | Single-cell RNA-seq |

---

## References

See [Overview](./README.md) for full details.
