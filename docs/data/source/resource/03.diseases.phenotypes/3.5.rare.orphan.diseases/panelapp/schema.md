---
id: schema-panelapp
title: "PanelApp Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, gene-panels, clinical-genomics, diagnostic-genes, rare-diseases, nhs]
---

# PanelApp Schema Documentation

**Document ID:** SCHEMA-PANELAPP

---

## TL;DR

PanelApp is a crowdsourced knowledge base for gene-disease relationships developed by Genomics England for NHS clinical genomics. It provides curated gene panels with expert-reviewed evidence levels (Green, Amber, Red) for diagnostic use in the 100,000 Genomes Project and NHS services.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Gene Panels | 350+ | Genomics England |
| Green Genes | 4,500+ | Genomics England |
| Expert Reviews | 50,000+ | Genomics England |
| Countries Using | 20+ | Genomics England |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | JSON, TSV, BED |
| API | Yes (REST API - free, no registration) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Panel ID | `[0-9]+` | 245 |
| Panel Name | Text | Intellectual disability |
| Gene Symbol | HGNC symbol | MECP2 |
| HGNC ID | `HGNC:[0-9]+` | HGNC:6990 |

---

## Evidence Levels

| Level | Color | Criteria |
|-------|-------|----------|
| 3 | Green | Diagnostic grade, high evidence |
| 2 | Amber | Moderate evidence, emerging |
| 1 | Red | Low evidence, not for diagnosis |
| 0 | Gray | No evidence yet |

---

## References

See [Overview](./README.md) for full details.
