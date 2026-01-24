---
id: variant.repositories
title: "Variant Repositories"
type: subcategory
parent: ../_index.md
last_updated: 2026-01-23
status: active
tags: [variants, clinical, snp, structural, ncbi]
---

# Variant Repositories

**Parent:** [Genetics & Genomics](../_index.md)

## Overview

Variant repositories are centralized databases that collect, curate, and distribute information about human genetic variants. These resources serve as foundational references for variant annotation, providing identifiers, genomic coordinates, and clinical interpretations.

The primary variant repositories include ClinVar for clinically-interpreted variants, dbSNP for short sequence variants, and dbVar for structural variants. Together they provide comprehensive coverage of the human variant landscape.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [ClinVar](./clinvar/_index.md) | 1 | Clinical variant interpretations |
| [dbSNP](./dbsnp/_index.md) | 1 | Short nucleotide polymorphisms |
| [dbVar](./dbvar/_index.md) | 1 | Structural variants |

## Integration Notes

ClinVar provides clinical significance; dbSNP provides population frequency data for SNVs; dbVar covers larger structural variants. Use rsIDs from dbSNP to cross-reference with population databases. ClinVar VCV/RCV identifiers link to specific variant-condition assertions.

---

## Database Selection Guide

### TL;DR

Primary genetics databases provide population allele frequencies, functional variant annotations, epigenetic marks, and structural variant catalogs essential for personalized health intelligence. dbNSFP v4.9 (35 prediction scores) and AlphaMissense (71M missense variants, CC BY 4.0) are recommended as Tier 1 priorities for functional annotation. gnomAD-SV v4.1 provides population-scale structural variant frequencies.

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Primary functional annotation | dbNSFP v4.9 | 35 prediction scores in single download; includes AlphaMissense, SpliceAI, CADD, REVEL |
| Missense pathogenicity | AlphaMissense | 71M variants scored; CC BY 4.0 allows commercial use; state-of-art accuracy |
| Population frequencies (diverse) | TOPMed BRAVO | 868M+ variants from 150K genomes; ~60% non-European ancestry |
| Population frequencies (NCBI) | ALFA R4 | Public domain; 200K+ subjects; excellent ClinVar coverage (960K+ RSIDs) |
| Structural variants | gnomAD-SV v4.1 | 1.2M SVs from 527K individuals; GraphQL API; open access |
| Regulatory annotation | ENCODE 4 | 926K human cCREs; comprehensive ChIP-seq, DNase-seq, ATAC-seq |
| Splice prediction | SpliceAI | Deep learning splice predictions; pre-computed scores available |

### Tier Priorities

| Tier | Databases | Priority Reason |
|------|-----------|-----------------|
| **1 (Critical)** | dbNSFP v4.9, AlphaMissense, gnomAD-SV v4.1, ALFA R4 | Core annotation; open/commercial-friendly licenses |
| **2 (High Value)** | TOPMed BRAVO, ENCODE 4, MaveDB, ClinGen Dosage | Diverse populations; regulatory elements; experimental validation |
| **3 (Enrichment)** | GenomeAsia 100K, H3Africa/AGVP, 4D Nucleome, FANTOM5 | Population-specific; 3D genome context |

### Licensing Summary

| License Type | Databases | Commercial Use |
|--------------|-----------|----------------|
| CC BY 4.0 | AlphaMissense, FANTOM5 | Yes (with attribution) |
| Public Domain | ALFA | Yes |
| Open Access | gnomAD-SV, ENCODE, MaveDB | Yes |
| Academic Only | dbNSFP (v4.9a), CADD, REVEL | No |
