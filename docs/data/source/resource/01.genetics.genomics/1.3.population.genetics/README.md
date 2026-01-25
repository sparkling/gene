---
id: population.genetics
title: "Population Genetics"
type: subcategory
parent: ../README.md
last_updated: 2026-01-23
status: active
tags: [population, allele-frequency, gnomad, biobank, diversity]
---

# Population Genetics

**Parent:** [Genetics & Genomics](../README.md)

## Overview

Population genetics databases aggregate variant frequency data from large-scale sequencing projects across diverse human populations. These resources are essential for filtering common variants in disease studies and understanding human genetic diversity.

Key resources include gnomAD (aggregated exome/genome data), UK Biobank (large-scale biobank), and TOPMed (deep whole-genome sequencing). The 1000 Genomes Project provides foundational population reference data.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [gnomAD](./gnomad/README.md) | 1 | Genome Aggregation Database |
| [1000 Genomes](./1000.genomes/README.md) | 1 | International population reference |
| [UK Biobank](./uk.biobank/README.md) | 1 | Large-scale biobank genetics |
| [TOPMed](./topmed/README.md) | 1 | Trans-Omics for Precision Medicine |

## Integration Notes

gnomAD is the primary resource for population allele frequencies. Filter variants with AF > 0.01 for rare disease analysis. Use population-specific frequencies for ancestry-matched studies. UK Biobank requires approved access. TOPMed provides deep WGS coverage.

---

## Database Selection Guide

### TL;DR

Population genetics databases aggregate variant frequency data from large-scale sequencing projects across diverse human populations. TOPMed BRAVO (868M+ variants, ~60% non-European), ALFA R4 (public domain, 960K+ ClinVar RSIDs), and gnomAD v4.1 (807K individuals) provide comprehensive global coverage. GenomeAsia 100K and H3Africa address critical gaps in Asian and African diversity.

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Primary population frequencies | gnomAD v4.1 | 807K+ individuals; 8 ancestry groups; open access |
| Diverse rare variants | TOPMed BRAVO | 868M+ variants; ~60% non-European; Freeze 10 |
| NCBI-integrated frequencies | ALFA R4 | Public domain; doubled cohort from R3; 960K+ ClinVar coverage |
| Large biobank | UK Biobank AFB | 490K WGS (full); public AFB for summary data |
| Asian diversity | GenomeAsia 100K | 66M+ variants from 219 Asian population groups |
| African diversity | H3Africa/AGVP | 3M+ novel variants; 50 ethnolinguistic groups |

### Comparative Statistics

| Database | Variants | Individuals | Non-European | License |
|----------|----------|-------------|--------------|---------|
| gnomAD v4.1 | 786M+ | 807K | ~40% | Open Access |
| TOPMed BRAVO | 868M+ | 150K | ~60% | Public (summary) |
| All of Us | 1B+ | 245K+ | ~50% | Data Use Agreement |
| ALFA R4 | - | 200K+ | 12 populations | Public Domain |
| UK Biobank AFB | - | 490K | European focus | Public (AFB) |
| GenomeAsia 100K | 66M+ | 1,739 | 100% Asian | Data Access Agreement |
| H3Africa/AGVP | 3M+ novel | 426+ WGS | 100% African | Controlled Access |

### Tier Priorities

| Tier | Databases | Priority Reason |
|------|-----------|-----------------|
| **1 (Critical)** | gnomAD v4.1, ALFA R4 | Broad coverage; NCBI integration; open licenses |
| **2 (High Value)** | TOPMed BRAVO, UK Biobank AFB, All of Us | Deep sequencing; rare variants; diverse populations |
| **3 (Enrichment)** | GenomeAsia 100K, H3Africa/AGVP, jMorp | Population-specific reference panels |

### API Availability

| Database | API Type | Rate Limits | Bulk Download |
|----------|----------|-------------|---------------|
| gnomAD | GraphQL | None stated | Yes |
| TOPMed BRAVO | REST | Unknown | Limited |
| ALFA | E-utilities | 3-10/sec | Yes |
| All of Us | Workbench APIs | Per-project | Registration required |
