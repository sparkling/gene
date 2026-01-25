---
id: ipd.imgt.hla
title: "IPD-IMGT/HLA - Immuno Polymorphism Database"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: autoimmune.inflammatory
tags:
  - hla
  - immunogenetics
  - transplantation
  - alleles
  - mhc
  - histocompatibility
---

# IPD-IMGT/HLA - Immuno Polymorphism Database

## Overview

The IPD-IMGT/HLA Database is the international repository of HLA (Human Leukocyte Antigen) gene sequences, providing a reference for human MHC (Major Histocompatibility Complex) allele nomenclature. The database is critical for transplantation matching, disease association studies, pharmacogenomics, and understanding immune system genetics.

Maintained by the EMBL-EBI in collaboration with the WHO Nomenclature Committee for Factors of the HLA System, IPD-IMGT/HLA contains over 30,000 HLA allele sequences across the classical (HLA-A, -B, -C, -DRB1, -DQB1, -DPB1) and non-classical HLA genes. The database follows a structured naming convention that captures allele groups, specific proteins, synonymous variants, and non-coding differences.

HLA genes are the most polymorphic in the human genome and play essential roles in immune recognition. The extensive variation catalogued in IPD-IMGT/HLA supports precision medicine applications including organ transplant compatibility, drug hypersensitivity prediction, and autoimmune disease risk assessment.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Alleles | 30,000+ |
| HLA-A Alleles | 7,000+ |
| HLA-B Alleles | 8,500+ |
| HLA-C Alleles | 7,000+ |
| HLA-DRB1 Alleles | 3,500+ |

## Primary Use Cases

1. Transplant donor-recipient matching
2. HLA typing result interpretation
3. Disease association study reference
4. Pharmacogenomics (drug hypersensitivity)
5. Vaccine and immunotherapy development

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| HLA Allele | `HLA-[A-Z]+[0-9]*\*[0-9:]+` | HLA-A*02:01:01:01 |
| Allele Group | `HLA-[A-Z]+\*[0-9]+` | HLA-A*02 |
| IPD Accession | `HLA[0-9]+` | HLA00001 |
| Gene Symbol | `HLA-[A-Z]+[0-9]*` | HLA-DRB1 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| IPD Website | https://www.ebi.ac.uk/ipd/imgt/hla/ | Web interface |
| Allele Query | https://www.ebi.ac.uk/ipd/imgt/hla/allele.html | Search tool |
| FTP | https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/ | Bulk downloads |
| REST API | https://www.ebi.ac.uk/ipd/imgt/hla/api/ | Programmatic |

## Data Formats

| Format | File | Description |
|--------|------|-------------|
| FASTA | hla_nuc.fasta | Nucleotide sequences |
| FASTA | hla_prot.fasta | Protein sequences |
| XML | hla.xml | Full allele data |
| DAT | hla.dat | EMBL format |

## Naming Convention

| Field | Example | Meaning |
|-------|---------|---------|
| Gene | HLA-A | Gene locus |
| Allele Group | *02 | Serological equivalent |
| Specific Protein | :01 | Protein sequence |
| Synonymous | :01 | Coding synonymous |
| Non-coding | :01 | Non-coding differences |

## Limitations

- HLA naming nomenclature complex for non-specialists
- Not all populations equally represented in allele discovery
- Functional impact of novel alleles often unknown
- Disease associations curated elsewhere (not in IPD)

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [ImmunoBase](../immunobase/README.md) - Autoimmune genetics
