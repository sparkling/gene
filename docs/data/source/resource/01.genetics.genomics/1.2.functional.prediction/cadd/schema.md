---
id: schema-cadd
title: "CADD Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, database, deleteriousness, annotation, prediction]
---

# CADD Schema Documentation

**Document ID:** SCHEMA-CADD
**Version:** 1.0
**Source Version:** v1.7 (GRCh38)

---

## TL;DR

CADD provides deleteriousness scores for all possible SNVs and indels in the human genome. Pre-computed files contain genomic coordinates with raw scores and PHRED-scaled scores, where higher values indicate more likely deleteriousness.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| SNV Coverage | ~9 billion | All possible |
| Features Integrated | 60+ | Annotations |
| Current Version | v1.7 | Release |
| Reference Builds | GRCh38, GRCh37 | Available |
| Pre-scored File Size | ~80 GB | GRCh38 |

---

## Entity Relationship Overview

```
┌──────────────────┐
│   Variant        │
├──────────────────┤
│ chromosome       │
│ position         │
│ reference        │
│ alternate        │
│ raw_score        │
│ phred_score      │
└──────────────────┘
```

---

## Core Tables/Entities

### Variant Score Record

**Description:** Single variant with CADD scores

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Chrom | string | Yes | Chromosome |
| Pos | integer | Yes | Position (1-based) |
| Ref | string | Yes | Reference allele |
| Alt | string | Yes | Alternate allele |
| RawScore | float | Yes | Untransformed SVM score |
| PHRED | float | Yes | PHRED-scaled score |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | TSV (tab-separated) |
| Alternative | VCF annotations |
| Encoding | UTF-8 |
| Compression | bgzip with tabix index |

---

## Sample Record

```tsv
#Chrom	Pos	Ref	Alt	RawScore	PHRED
1	10001	T	A	0.234567	8.123
1	10001	T	C	0.456789	12.456
1	10001	T	G	0.345678	10.234
```

---

## Score Interpretation

| PHRED Score | Percentile | Interpretation |
|-------------|------------|----------------|
| >= 30 | Top 0.1% | Very likely deleterious |
| >= 20 | Top 1% | Likely deleterious |
| >= 15 | Top 3% | Possibly deleterious |
| >= 10 | Top 10% | Mildly deleterious |
| < 10 | Bottom 90% | Likely benign |

---

## Glossary

| Term | Definition |
|------|------------|
| RawScore | Untransformed SVM output score |
| PHRED | -10 * log10(rank/total), scaled like quality scores |
| C-score | Alternative name for CADD score |

---

## References

1. https://cadd.gs.washington.edu/
2. Rentzsch et al. (2019) Nucleic Acids Res. DOI: 10.1093/nar/gky1016
