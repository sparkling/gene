---
id: omim
title: "Online Mendelian Inheritance in Man (OMIM)"
type: data-source
category: diseases
subcategory: phenotype.databases
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [mendelian, genetics, rare-diseases, gene-disease, inheritance]
---

# Online Mendelian Inheritance in Man (OMIM)

**Category:** [Diseases & Phenotypes](../../_index.md) > [Phenotype Databases](../_index.md)

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

## License

| Aspect | Value |
|--------|-------|
| License | Custom OMIM License |
| Commercial Use | Requires license |
| Academic Use | Free with registration |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [HPO](../hpo/_index.md) - Phenotype annotations linked to OMIM
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/_index.md) - Rare disease database
