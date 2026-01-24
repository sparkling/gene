---
id: schema-syngo
title: "SynGO Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, synapse, neuroscience, gene-ontology, synaptic-function, brain-disorders]
---

# SynGO Schema Documentation

**Document ID:** SCHEMA-SYNGO

---

## TL;DR

SynGO is an expert-curated resource for synaptic function gene annotations, providing a knowledge base for synaptic protein localization and function. It uses a structured ontology framework with over 1,200 annotated genes and 5,800 annotations supported by published experimental evidence.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Annotated Genes | 1,200+ | SynGO |
| Annotations | 5,800+ | SynGO |
| Supporting Papers | 2,800+ | SynGO |
| Ontology Terms | 180+ | SynGO |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | OBO, GAF, TSV, JSON |
| API | Yes (Web service) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| SynGO Term | `SYNGO:[0-9]+` | SYNGO:synapse |
| Gene Symbol | HGNC | SYN1 |
| UniProt ID | `[A-Z0-9]{6,10}` | P17600 |
| GO Term | `GO:[0-9]{7}` | GO:0045202 |

---

## Ontology Structure

| Branch | Description | Terms |
|--------|-------------|-------|
| Cellular Component | Synaptic localization | 90+ |
| Biological Process | Synaptic function | 90+ |
| Pre-synaptic | Axon terminal | 50+ |
| Post-synaptic | Dendritic spine | 60+ |

---

## Evidence Codes

| Code | Meaning |
|------|---------|
| ECO:0005593 | Biological assay |
| ECO:0005589 | Immunolocalization |
| ECO:0005644 | Electron microscopy |
| ECO:0006063 | Biochemical assay |

---

## References

See [Overview](./_index.md) for full details.
