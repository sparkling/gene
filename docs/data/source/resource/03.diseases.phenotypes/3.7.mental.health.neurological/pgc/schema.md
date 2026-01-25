---
id: schema-pgc
title: "Psychiatric Genomics Consortium (PGC) Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, psychiatric, gwas, schizophrenia, depression, bipolar, genetics]
---

# Psychiatric Genomics Consortium (PGC) Schema Documentation

**Document ID:** SCHEMA-PGC

---

## TL;DR

PGC is the largest international collaboration studying the genetic basis of psychiatric disorders. It conducts mega-analyses combining datasets from 800+ research groups worldwide, identifying hundreds of genetic loci and revealing extensive genetic overlap between different psychiatric conditions.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Participating Groups | 800+ | PGC |
| Total Samples | 4M+ | PGC |
| Published GWAS | 50+ | PGC |
| Genetic Loci Identified | 1,000+ | PGC |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | Summary statistics (TSV), LD scores |
| API | No (direct downloads) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| dbSNP ID | `rs[0-9]+` | rs1006737 |
| GWAS Study ID | Custom | PGC_SCZ3 |
| Gene Symbol | HGNC | CACNA1C |
| Disorder Code | Custom | SCZ, BIP, MDD |

---

## Major Studies

| Disorder | Sample Size | Loci | Reference |
|----------|-------------|------|-----------|
| Schizophrenia (SCZ3) | 320,000 | 287 | Trubetskoy 2022 |
| Major Depression (MDD3) | 1.2M | 243 | Howard 2019 |
| Bipolar Disorder (BIP3) | 413,000 | 64 | Mullins 2021 |
| ADHD | 225,000 | 27 | Demontis 2023 |

---

## References

See [Overview](./README.md) for full details.
