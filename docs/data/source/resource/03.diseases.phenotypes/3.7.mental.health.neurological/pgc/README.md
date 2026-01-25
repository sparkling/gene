---
id: pgc
title: "Psychiatric Genomics Consortium (PGC)"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: mental.health.neurological
tags:
  - psychiatric
  - gwas
  - schizophrenia
  - depression
  - bipolar
  - genetics
---

# Psychiatric Genomics Consortium (PGC)

## Overview

The Psychiatric Genomics Consortium (PGC) is the largest international collaboration studying the genetic basis of psychiatric disorders. By combining datasets from hundreds of research groups worldwide, PGC conducts mega-analyses with unprecedented statistical power to identify genetic variants associated with major mental illnesses.

Founded in 2007, PGC has conducted landmark GWAS for schizophrenia, bipolar disorder, major depressive disorder, ADHD, autism spectrum disorder, anorexia nervosa, and other conditions. The consortium's studies have identified hundreds of genetic loci associated with psychiatric disorders and revealed extensive genetic overlap between different conditions.

PGC makes summary statistics from published studies freely available, enabling researchers to perform downstream analyses including polygenic risk scoring, genetic correlation analysis, and pathway enrichment. The consortium also supports cross-disorder analyses examining shared genetic architecture across mental illnesses.

## Key Statistics

| Metric | Value |
|--------|-------|
| Participating Groups | 800+ |
| Total Samples | 4M+ |
| Countries | 40+ |
| Published GWAS | 50+ |
| Genetic Loci Identified | 1,000+ |

## Primary Use Cases

1. Psychiatric disorder genetic risk discovery
2. Polygenic risk score development
3. Cross-disorder genetic correlation analysis
4. Drug target identification for mental health
5. Biological pathway analysis

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| dbSNP ID | `rs[0-9]+` | rs1006737 |
| GWAS Study ID | Custom | PGC_SCZ3 |
| Gene Symbol | HGNC | CACNA1C |
| Disorder Code | Custom | SCZ, BIP, MDD |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| PGC Website | https://pgc.unc.edu/ | Main portal |
| Summary Stats | https://pgc.unc.edu/for-researchers/download-results/ | GWAS results |
| LD Hub | http://ldsc.broadinstitute.org/ | Analysis platform |
| GWAS Catalog | https://www.ebi.ac.uk/gwas/ | Study records |

## Major Studies

| Disorder | Sample Size | Loci | Reference |
|----------|-------------|------|-----------|
| Schizophrenia (SCZ3) | 320,000 | 287 | Trubetskoy 2022 |
| Major Depression (MDD3) | 1.2M | 243 | Howard 2019 |
| Bipolar Disorder (BIP3) | 413,000 | 64 | Mullins 2021 |
| ADHD | 225,000 | 27 | Demontis 2023 |
| Autism (ASD) | 54,000 | 12 | Grove 2019 |

## Data Types

| Data Type | Description | Access |
|-----------|-------------|--------|
| Summary Statistics | GWAS effect sizes | Open |
| LD Scores | Linkage disequilibrium | Open |
| PRS Weights | Polygenic scores | Open |
| Individual Data | Genotypes, phenotypes | Controlled (dbGaP) |

## Limitations

- Individual-level data requires dbGaP controlled access
- European ancestry overrepresented in most studies
- Phenotype definitions vary across contributing cohorts
- Rare variant data limited in GWAS format

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [Allen Brain Atlas](../allen.brain.atlas/README.md) - Brain expression
