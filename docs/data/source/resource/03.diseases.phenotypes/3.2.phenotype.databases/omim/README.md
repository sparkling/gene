---
id: omim
title: "Online Mendelian Inheritance in Man (OMIM)"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: phenotype.databases
tags:
  - mendelian
  - genetics
  - rare-diseases
  - gene-disease
  - inheritance
---

# Online Mendelian Inheritance in Man (OMIM)

## Overview

OMIM (Online Mendelian Inheritance in Man) is the authoritative compendium of human genes and genetic phenotypes, providing comprehensive expert-authored entries on all known Mendelian disorders and over 16,000 genes. Founded by Victor McKusick in 1966, OMIM has been continuously curated for over 50 years and represents the gold standard for clinical genetics information.

Each OMIM entry contains detailed clinical synopses, molecular genetics information, genotype-phenotype correlations, and extensive references to the scientific literature. Gene entries describe the gene's function, associated disorders, and allelic variants, while phenotype entries detail clinical features, inheritance patterns, and mapping information.

OMIM uses a structured MIM number system where prefixes indicate entry type: asterisk (*) for genes, number sign (#) for phenotypes with known molecular basis, plus (+) for gene/phenotype combinations, and percent (%) for phenotypes with unknown molecular basis. This classification system is widely used in clinical genetics reporting and genetic counseling.

## Key Statistics

| Metric | Value |
|--------|-------|
| Gene Entries | 16,000+ |
| Phenotype Entries | 9,000+ |
| Gene-Phenotype Relationships | 7,200+ |
| Allelic Variants | 29,000+ |
| Literature References | 180,000+ |

## Primary Use Cases

1. Clinical genetics reference for Mendelian disorders
2. Gene-disease relationship curation source
3. Variant interpretation in diagnostic laboratories
4. Phenotype-genotype correlation research
5. Medical genetics education and training

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| MIM Number (Gene) | `*[0-9]{6}` | *154700 |
| MIM Number (Phenotype) | `#[0-9]{6}` | #154700 |
| MIM Number (Combined) | `+[0-9]{6}` | +134610 |
| MIM Number (Unknown) | `%[0-9]{6}` | %176300 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| OMIM Website | https://omim.org/ | Official portal |
| OMIM API | https://api.omim.org/ | Requires registration |
| FTP Download | https://data.omim.org/ | Licensed access |
| Gene Map | https://omim.org/geneMap | Interactive map |

## Data Formats

| Format | File | Notes |
|--------|------|-------|
| Text | mim2gene.txt | Gene mappings |
| Text | genemap2.txt | Gene-phenotype map |
| Text | mimTitles.txt | MIM titles |
| JSON | API responses | Programmatic access |

## Limitations

- Commercial use requires paid license agreement
- API access requires registration and approval
- Focuses on Mendelian; complex disease coverage limited
- Full-text entries not available via bulk download

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [HPO](../hpo/README.md) - Phenotype annotations linked to OMIM
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/README.md) - Rare disease database
