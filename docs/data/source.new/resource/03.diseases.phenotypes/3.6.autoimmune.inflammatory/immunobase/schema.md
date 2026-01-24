---
id: schema-immunobase
title: "ImmunoBase Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, autoimmune, gwas, genetic-associations, immune-diseases, snps]
---

# ImmunoBase Schema Documentation

**Document ID:** SCHEMA-IMMUNOBASE

---

## TL;DR

ImmunoBase is a curated database of genetic associations for 12 major immune-mediated diseases, integrating data from GWAS, ImmunoChip studies, and fine-mapping analyses. It provides a unified view of genetic variation across autoimmune and inflammatory conditions.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Diseases Covered | 12 | ImmunoBase |
| Associated Regions | 350+ | ImmunoBase |
| SNP Associations | 1,000+ | ImmunoBase |
| Credible Sets | 200+ | ImmunoBase |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | Web interface, bulk downloads |
| API | No |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| dbSNP ID | `rs[0-9]+` | rs2476601 |
| Region ID | Custom | 1p13.2 |
| Gene Symbol | HGNC | PTPN22 |
| Disease Code | Custom | T1D, RA, CEL |

---

## Diseases Covered

| Disease | Abbreviation | GWAS Loci |
|---------|--------------|-----------|
| Type 1 Diabetes | T1D | 60+ |
| Rheumatoid Arthritis | RA | 100+ |
| Celiac Disease | CEL | 40+ |
| Multiple Sclerosis | MS | 200+ |
| Crohn's Disease | CD | 140+ |
| Ulcerative Colitis | UC | 100+ |

---

## References

See [Overview](./_index.md) for full details.
