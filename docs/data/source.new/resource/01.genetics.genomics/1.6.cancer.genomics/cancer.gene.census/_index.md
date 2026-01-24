---
id: cancer-gene-census
title: "Cancer Gene Census"
type: data-source
category: genetics
subcategory: cancer.genomics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [cancer, driver, oncogene, tumor-suppressor, census]
---

# Cancer Gene Census

**Category:** [Genetics & Genomics](../../_index.md) > [Cancer Genomics](../_index.md)

## Overview

The Cancer Gene Census (CGC) is an ongoing effort to catalogue genes whose mutations are causally implicated in cancer. Maintained by the Wellcome Sanger Institute as part of COSMIC, it distinguishes between Tier 1 genes with documented cancer-causing mutations and Tier 2 genes with strong evidence of cancer involvement.

Each CGC entry includes the gene's role (oncogene, tumor suppressor, or fusion partner), associated cancer types, mutation types commonly observed, and supporting evidence from the literature. The census provides a gold-standard reference for cancer gene identification and driver mutation analysis.

The CGC is continuously updated as new cancer genes are identified through research and clinical studies. It serves as a foundation for cancer gene panel design, variant interpretation, and research prioritization.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Genes | 736 |
| Tier 1 Genes | 588 |
| Tier 2 Genes | 148 |
| Oncogenes | 300+ |
| Tumor Suppressors | 300+ |

## Primary Use Cases

1. Cancer panel gene selection
2. Driver mutation identification
3. Oncogene/TSG classification
4. Research prioritization

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Gene | HUGO symbol | TP53 |
| Tier | 1 or 2 | Tier 1 |
| Role | Oncogene/TSG/Fusion | TSG |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web | https://cancer.sanger.ac.uk/census | Interactive list |
| COSMIC | Part of COSMIC | Integrated access |
| Downloads | Registered access | TSV format |

## License

| Aspect | Value |
|--------|-------|
| License | COSMIC license terms |
| Commercial Use | Requires license |
| Academic Use | Free with registration |

## See Also

- [COSMIC](../cosmic/_index.md) - Somatic mutations
- [OncoKB](../oncokb/_index.md) - Actionable genes
- [CIViC](../civic/_index.md) - Clinical interpretations
