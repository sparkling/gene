---
id: schema-decipher
title: "DECIPHER Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, rare-diseases, cnv, developmental-disorders, clinical-genetics, genomic-variants]
---

# DECIPHER Schema Documentation

**Document ID:** SCHEMA-DECIPHER

---

## TL;DR

DECIPHER is an interactive database for clinical interpretation of rare CNVs and sequence variants. It collects anonymized patient phenotypes and genomic variants from clinical genetics centers across 35 countries, using HPO terminology for phenotype annotation.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Patient Records | 50,000+ | DECIPHER |
| CNV Records | 40,000+ | DECIPHER |
| SNV Records | 60,000+ | DECIPHER |
| Contributing Centers | 300+ | DECIPHER |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | Web interface, API |
| API | Yes |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| DECIPHER ID | `[0-9]+` | 250000 |
| Patient ID | `DECIPHER:[0-9]+` | DECIPHER:250001 |
| HPO Term | `HP:[0-9]{7}` | HP:0001249 |
| Ensembl Gene | `ENSG[0-9]{11}` | ENSG00000146648 |

---

## Scoring Metrics

| Metric | Description |
|--------|-------------|
| HI Score | Haploinsufficiency prediction |
| pLI | Probability of LoF intolerance |
| LOEUF | LoF observed/expected upper bound |
| Pathogenicity | Variant classification |

---

## References

See [Overview](./README.md) for full details.
