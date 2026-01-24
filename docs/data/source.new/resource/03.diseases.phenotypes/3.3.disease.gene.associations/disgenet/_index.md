---
id: disgenet
title: "DisGeNET - Disease Gene Network"
type: data-source
category: diseases
subcategory: disease.gene.associations
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [gene-disease, variant-disease, associations, scoring, text-mining]
---

# DisGeNET - Disease Gene Network

**Category:** [Diseases & Phenotypes](../../_index.md) > [Disease Gene Associations](../_index.md)

## Overview

DisGeNET is the largest publicly available collection of genes and variants associated with human diseases. The platform integrates data from expert-curated repositories, GWAS catalogues, animal models, and the scientific literature to provide a comprehensive resource for understanding genotype-phenotype relationships.

The database contains over 628,000 gene-disease associations (GDAs) involving 17,500 genes and 24,000 diseases, as well as 210,000 variant-disease associations (VDAs). Each association is scored using the DisGeNET scoring system, which considers the number and type of supporting sources, providing a confidence metric for prioritization.

DisGeNET uses UMLS Concept Unique Identifiers (CUIs) as the primary disease vocabulary, enabling broad coverage and integration with clinical terminologies. The platform also provides disease semantic types, MeSH classifications, and cross-references to HPO, DOID, and Orphanet for flexible disease filtering.

## Key Statistics

| Metric | Value |
|--------|-------|
| Gene-Disease Associations | 628,685+ |
| Variant-Disease Associations | 210,498+ |
| Genes | 17,549 |
| Variants | 117,337 |
| Diseases/Phenotypes | 24,166 |

## Primary Use Cases

1. Gene prioritization for disease research
2. Variant interpretation in clinical genomics
3. Drug target identification and validation
4. Biomarker discovery and validation
5. Network-based disease analysis

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| UMLS CUI | `C[0-9]{7}` | C0020443 (Hypercholesterolemia) |
| NCBI Gene ID | `[0-9]+` | 3157 |
| dbSNP ID | `rs[0-9]+` | rs2066843 |
| Gene Symbol | HGNC symbol | HMGCS1 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| DisGeNET Web | https://www.disgenet.org/ | Web interface |
| REST API | https://www.disgenet.org/api | API access |
| Cytoscape App | https://apps.cytoscape.org/apps/disgenetapp | Network analysis |
| R Package | disgenet2r | Bioconductor |

## Data Formats

| Format | File | Description |
|--------|------|-------------|
| TSV | curated_gene_disease_associations.tsv | Expert curated |
| TSV | all_gene_disease_associations.tsv | Complete dataset |
| TSV | variant_disease_associations.tsv | VDA file |
| SQLite | disgenet_2020.db | Database dump |

## Score Metrics

| Score | Range | Interpretation |
|-------|-------|----------------|
| GDA Score | 0-1 | Association confidence |
| Evidence Index (EI) | 0-1 | Evidence consistency |
| DSI | 0-1 | Disease Specificity Index |
| DPI | 0-1 | Disease Pleiotropy Index |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY-NC-SA 4.0 |
| Commercial Use | Requires separate license |
| Academic Use | Free |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [Open Targets](../open.targets/_index.md) - Complementary target platform
- [GWAS Catalog](../../../01.genetics.genomics/1.5.expression.regulation/gwas.catalog/_index.md) - GWAS data source
