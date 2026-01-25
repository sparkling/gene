---
id: refseq
title: "NCBI RefSeq - Reference Sequence Database"
type: data-source
category: proteins
subcategory: protein.sequences.annotations
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [sequences, ncbi, reference, proteins, transcripts, genomes]
---

# NCBI RefSeq - Reference Sequence Database

**Category:** [Proteins](../../../README.md) > [Protein Sequences](../README.md)

## Overview

NCBI Reference Sequence (RefSeq) database provides a comprehensive, integrated, non-redundant, well-annotated set of reference sequences including genomic DNA, transcripts, and proteins. RefSeq sequences are explicitly linked to their source data and include both computationally derived and manually curated records.

RefSeq is organized into multiple categories including complete genomic molecules, model organisms, and reference genomes. For proteins, it provides stable reference records that are used throughout NCBI resources and widely cited in the literature.

RefSeq is essential for genome annotation, sequence comparison, primer design, and as a standard reference for variant interpretation.

## Key Statistics

| Metric | Value |
|--------|-------|
| Protein Records | 400M+ |
| Transcript Records | 200M+ |
| Organisms | 120,000+ |
| Human Proteins | 110,000+ |
| Last Update | Bimonthly |

## Primary Use Cases

1. Reference sequence for variant annotation
2. Genome annotation pipeline
3. Sequence comparison and alignment
4. Primer and probe design
5. Cross-database linking

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| RefSeq Protein | `[ANYXWZ]P_[0-9]+` | NP_000001.1 |
| RefSeq mRNA | `[NX]M_[0-9]+` | NM_000001.3 |
| RefSeq ncRNA | `[NX]R_[0-9]+` | NR_000001.1 |
| RefSeq Genome | `N[CZW]_[0-9]+` | NC_000001.11 |
| GI Number | Numeric (deprecated) | 4557225 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.ncbi.nlm.nih.gov/refseq | Search and browse |
| Entrez API | https://eutils.ncbi.nlm.nih.gov/entrez/eutils | E-utilities |
| FTP | https://ftp.ncbi.nlm.nih.gov/refseq | Full downloads |
| BLAST | https://blast.ncbi.nlm.nih.gov | Sequence search |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [UniProt](../uniprot/README.md) - Universal Protein Resource
- [NCBI Gene](../../../../01.genetics.genomics/1.2.gene.databases/ncbi.gene/README.md) - Gene records
