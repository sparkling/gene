---
id: immunobase
title: "ImmunoBase - Autoimmune Disease Genetics Resource"
type: data-source
category: diseases
subcategory: autoimmune.inflammatory
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [autoimmune, gwas, genetic-associations, immune-diseases, snps]
---

# ImmunoBase - Autoimmune Disease Genetics Resource

**Category:** [Diseases & Phenotypes](../../_index.md) > [Autoimmune & Inflammatory](../_index.md)

## Overview

ImmunoBase is a curated database of genetic associations for immune-mediated diseases, integrating data from GWAS, ImmunoChip studies, and fine-mapping analyses. The resource focuses on autoimmune and inflammatory conditions where genetic studies have revealed shared and disease-specific susceptibility variants.

Developed as a collaboration between the Wellcome Trust Sanger Institute and the University of Cambridge, ImmunoBase provides a unified view of genetic variation across 12 major autoimmune diseases including type 1 diabetes, rheumatoid arthritis, celiac disease, multiple sclerosis, and inflammatory bowel disease. The platform emphasizes the interconnected genetic architecture of immune diseases.

ImmunoBase enables researchers to explore associated regions, candidate genes, and their overlap across diseases. The database includes detailed annotation of associated variants, credible sets from fine-mapping studies, and integration with functional genomics data to support causal gene and variant identification.

## Key Statistics

| Metric | Value |
|--------|-------|
| Diseases Covered | 12 |
| Associated Regions | 350+ |
| SNP Associations | 1,000+ |
| Credible Sets | 200+ |
| Studies Integrated | 50+ |

## Primary Use Cases

1. Autoimmune disease genetic risk exploration
2. Cross-disease association analysis
3. Fine-mapping variant prioritization
4. Shared pathway identification
5. Drug target discovery in autoimmunity

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| dbSNP ID | `rs[0-9]+` | rs2476601 |
| Region ID | Custom | 1p13.2 |
| Gene Symbol | HGNC | PTPN22 |
| Disease Code | Custom | T1D, RA, CEL |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| ImmunoBase Web | https://www.immunobase.org/ | Web interface |
| Region Browser | https://www.immunobase.org/region/ | Genomic view |
| Downloads | https://www.immunobase.org/downloads/ | Bulk data |

## Diseases Covered

| Disease | Abbreviation | GWAS Loci |
|---------|--------------|-----------|
| Type 1 Diabetes | T1D | 60+ |
| Rheumatoid Arthritis | RA | 100+ |
| Celiac Disease | CEL | 40+ |
| Multiple Sclerosis | MS | 200+ |
| Crohn's Disease | CD | 140+ |
| Ulcerative Colitis | UC | 100+ |
| Ankylosing Spondylitis | AS | 30+ |
| Psoriasis | PSO | 60+ |

## Data Types

| Data Type | Description |
|-----------|-------------|
| Association Summary | P-values, odds ratios |
| Credible Sets | Fine-mapped variants |
| Gene Annotations | Candidate genes |
| Cross-disease | Shared associations |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [GWAS Catalog](../../../01.genetics.genomics/1.5.expression.regulation/gwas.catalog/_index.md) - Broader GWAS data
- [IPD-IMGT/HLA](../ipd.imgt.hla/_index.md) - HLA immunogenetics
