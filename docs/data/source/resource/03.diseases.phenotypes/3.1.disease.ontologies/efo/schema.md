---
id: schema-efo
title: "Experimental Factor Ontology (EFO) Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, ontology, traits, diseases, experimental-factors, gwas, ebi]
---

# Experimental Factor Ontology (EFO) Schema Documentation

**Document ID:** SCHEMA-EFO

---

## TL;DR

EFO provides a systematic description of experimental variables for EBI databases, integrating terms from multiple domain ontologies including Disease Ontology, HPO, and ChEBI. It serves as the key ontology for describing diseases, phenotypes, and traits in GWAS Catalog and Open Targets.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Total Terms | ~50,000+ | EBI |
| Disease Terms | ~20,000 | EBI |
| Phenotype Terms | ~10,000 | EBI |
| Release Cycle | Monthly | EBI |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | OWL, OBO, JSON |
| API | Yes (OLS Browser) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| EFO ID | `EFO_[0-9]{7}` | EFO_0000685 |
| Orphanet XREF | `Orphanet:[0-9]+` | Orphanet:558 |
| DOID XREF | `DOID:[0-9]+` | DOID:14323 |
| HP XREF | `HP:[0-9]{7}` | HP:0001250 |

---

## References

See [Overview](./README.md) for full details.
