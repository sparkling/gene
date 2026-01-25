---
id: uniprot
title: "UniProt - Universal Protein Resource"
type: data-source
category: proteins
subcategory: protein.sequences.annotations
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [proteins, sequences, annotations, swiss-prot, trembl, curation]
---

# UniProt - Universal Protein Resource

**Category:** [Proteins](../../../README.md) > [Protein Sequences](../README.md)

## Overview

The Universal Protein Resource (UniProt) is the most comprehensive, high-quality, and freely accessible database of protein sequence and functional annotation. It is a collaboration between the European Bioinformatics Institute (EMBL-EBI), the Swiss Institute of Bioinformatics (SIB), and the Protein Information Resource (PIR).

UniProt consists of two main sections: Swiss-Prot (manually annotated and reviewed) and TrEMBL (automatically annotated and unreviewed). It provides detailed protein annotations including function, domains, post-translational modifications, variants, and disease associations.

UniProt is the gold standard for protein annotation and is extensively cross-referenced by virtually all major biological databases.

## Key Statistics

| Metric | Value |
|--------|-------|
| Swiss-Prot (reviewed) | 570,000+ |
| TrEMBL (unreviewed) | 250M+ |
| Organisms | 20,000+ |
| Human Proteins | 20,400 (Swiss-Prot) |
| Last Update | Monthly |

## Primary Use Cases

1. Protein function annotation
2. Sequence retrieval and analysis
3. Domain and family classification
4. Variant interpretation
5. Cross-database linking

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| UniProt Accession | `[OPQ][0-9][A-Z0-9]{3}[0-9]` | P04637 |
| UniProt ID | `[A-Z]+_[A-Z]+` | P53_HUMAN |
| Entry Name | `[A-Z0-9]+_[A-Z]+` | TP53_HUMAN |
| Gene Name | Text | TP53 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.uniprot.org | Search, browse, BLAST |
| REST API | https://rest.uniprot.org | Programmatic access |
| SPARQL | https://sparql.uniprot.org | RDF queries |
| FTP | https://ftp.uniprot.org | Full downloads |
| ID Mapping | https://www.uniprot.org/id-mapping | Identifier conversion |

## License

| Aspect | Value |
|--------|-------|
| License | Creative Commons Attribution 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Sample Data](./sample.json) - Example records
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - Transform to unified schema
- [RefSeq](../refseq/README.md) - NCBI reference sequences
- [UniProt ID Mapping](../../../08.literature.knowledge/8.3.identifier.mapping/uniprot.id.mapping/README.md) - ID conversion
