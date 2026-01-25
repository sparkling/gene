---
id: schema-ipd-imgt-hla
title: "IPD-IMGT/HLA Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: draft
tags: [schema, hla, immunogenetics, transplantation, alleles, mhc, histocompatibility]
---

# IPD-IMGT/HLA Schema Documentation

**Document ID:** SCHEMA-IPD-IMGT-HLA

---

## TL;DR

IPD-IMGT/HLA is the international repository of HLA gene sequences, providing a reference for human MHC allele nomenclature. With over 30,000 allele sequences, it is critical for transplantation matching, disease association studies, and pharmacogenomics.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Total Alleles | 30,000+ | EMBL-EBI |
| HLA-A Alleles | 7,000+ | EMBL-EBI |
| HLA-B Alleles | 8,500+ | EMBL-EBI |
| HLA-C Alleles | 7,000+ | EMBL-EBI |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | FASTA, XML, DAT (EMBL format) |
| API | Yes (REST API) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| HLA Allele | `HLA-[A-Z]+[0-9]*\*[0-9:]+` | HLA-A*02:01:01:01 |
| Allele Group | `HLA-[A-Z]+\*[0-9]+` | HLA-A*02 |
| IPD Accession | `HLA[0-9]+` | HLA00001 |
| Gene Symbol | `HLA-[A-Z]+[0-9]*` | HLA-DRB1 |

---

## Naming Convention

| Field | Example | Meaning |
|-------|---------|---------|
| Gene | HLA-A | Gene locus |
| Allele Group | *02 | Serological equivalent |
| Specific Protein | :01 | Protein sequence |
| Synonymous | :01 | Coding synonymous |
| Non-coding | :01 | Non-coding differences |

---

## References

See [Overview](./README.md) for full details.
