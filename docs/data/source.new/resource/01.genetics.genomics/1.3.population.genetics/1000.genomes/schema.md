---
id: schema-1000-genomes
title: "1000 Genomes Project Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, database, population, frequency, diversity]
---

# 1000 Genomes Project Schema Documentation

**Document ID:** SCHEMA-1000-GENOMES
**Version:** 1.0
**Source Version:** Phase 3 / 30x High-Coverage

---

## TL;DR

The 1000 Genomes Project provides population-level allele frequencies for variants discovered in 2,504 individuals across 26 global populations. Data is available in VCF format with per-population and super-population frequency annotations.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Individuals | 2,504 | Phase 3 |
| Populations | 26 | 5 continents |
| Super-populations | 5 | AFR, AMR, EAS, EUR, SAS |
| SNVs | 84,700,000+ | Phase 3 |
| Indels | 3,600,000+ | Phase 3 |

---

## Entity Relationship Overview

```
┌────────────────────────────────────────────────────────┐
│                      Variant                            │
├────────────────────────────────────────────────────────┤
│ Position → Global AF → Super-pop AF → Population AF    │
│ (chr:pos)   (AF)        (AFR_AF, etc)   (YRI_AF, etc) │
└────────────────────────────────────────────────────────┘
```

---

## Core Tables/Entities

### VCF INFO Fields

**Description:** Population frequency annotations

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| AF | float | Yes | Global allele frequency |
| AC | integer | Yes | Allele count |
| AN | integer | Yes | Total alleles |
| AFR_AF | float | No | African super-population AF |
| AMR_AF | float | No | American super-population AF |
| EAS_AF | float | No | East Asian super-population AF |
| EUR_AF | float | No | European super-population AF |
| SAS_AF | float | No | South Asian super-population AF |

### Population Codes

| Super-pop | Population | Description |
|-----------|------------|-------------|
| AFR | YRI | Yoruba in Ibadan, Nigeria |
| AFR | LWK | Luhya in Webuye, Kenya |
| EUR | GBR | British in England and Scotland |
| EUR | FIN | Finnish in Finland |
| EAS | CHB | Han Chinese in Beijing |
| EAS | JPT | Japanese in Tokyo |
| SAS | GIH | Gujarati Indians in Houston |
| AMR | MXL | Mexican Ancestry in LA |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | VCF (Variant Call Format) |
| Alternative | BCF (binary VCF) |
| Encoding | UTF-8 |
| Compression | bgzip with tabix index |

---

## Sample Record

```vcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10177	rs367896724	A	AC	100	PASS	AF=0.425;AC=2130;AN=5008;AFR_AF=0.49;AMR_AF=0.34;EAS_AF=0.40;EUR_AF=0.38;SAS_AF=0.48
```

---

## Glossary

| Term | Definition |
|------|------------|
| AF | Allele Frequency (0-1) |
| AC | Allele Count in called genotypes |
| AN | Total number of alleles in called genotypes |
| Super-population | Continental grouping of populations |

---

## References

1. https://www.internationalgenome.org/
2. 1000 Genomes Consortium (2015) Nature. DOI: 10.1038/nature15393
