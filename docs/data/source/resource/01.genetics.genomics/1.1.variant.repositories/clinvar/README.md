---
id: clinvar
title: "ClinVar"
type: data-source
category: genetics
subcategory: variant.repositories
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [variant, clinical, pathogenicity, ncbi, germline, somatic]
---

# ClinVar

**Category:** [Genetics & Genomics](../../_index.md) > [Variant Repositories](../_index.md)

## Overview

ClinVar is NCBI's public archive of interpretations of clinically relevant variants. It aggregates information about genomic variation and its relationship to human health, collecting submissions from clinical laboratories, research groups, and expert panels worldwide.

The database uses a three-tier accession system: VCV (Variation Archive) for aggregate interpretations across all submitters, RCV (Record) for individual variant-condition pairs, and SCV (Submitted Record) for individual submissions. This structure enables tracking of both consensus interpretations and the underlying evidence from multiple sources.

ClinVar integrates with other NCBI resources including dbSNP, dbVar, and MedGen, providing comprehensive variant annotation including genomic coordinates, HGVS nomenclature, gene associations, and disease relationships with standardized phenotype terms.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Submissions | 10,000,000+ |
| Unique Variants | 1,500,000+ |
| Top Submitter | LabCorp Genetics (1.89M) |
| Reference Genome | GRCh38, GRCh37 |
| Update Frequency | Weekly (Mondays) |

## Primary Use Cases

1. Clinical variant interpretation and pathogenicity assessment
2. Diagnostic laboratory variant classification lookup
3. Research variant prioritization and filtering
4. Gene-disease association discovery

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| VCV | VCV + 9 digits + version | VCV000000123.4 |
| RCV | RCV + 9 digits + version | RCV000012345.6 |
| SCV | SCV + 9 digits + version | SCV000123456.7 |
| AlleleID | Integer | 12345 |
| VariationID | Integer | 67890 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| FTP | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/ | VCF, XML, TSV files |
| Web | https://www.ncbi.nlm.nih.gov/clinvar/ | Interactive search |
| Entrez | esearch/efetch | Programmatic access |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (CC BY) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [dbSNP](../dbsnp/_index.md) - SNP identifiers
- [dbVar](../dbvar/_index.md) - Structural variants
