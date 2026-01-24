---
id: dgidb
title: "DGIdb - Drug Gene Interaction Database"
type: data-source
category: compounds
subcategory: compound-target-interactions
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [drug-gene-interactions, druggable-genome, targets, pharmacogenomics, graphql]
---

# DGIdb - Drug Gene Interaction Database

**Category:** [Compounds & Molecules](../../_index.md) > [Compound-Target Interactions](../_index.md)

## Overview

DGIdb (Drug Gene Interaction Database) aggregates drug-gene interaction data from multiple trusted sources into a unified, queryable resource. The database provides information about how drugs interact with genes and their protein products, along with categorization of genes into druggable classes.

DGIdb integrates data from over 40 sources including DrugBank, PharmGKB, ChEMBL, TTD, and clinical databases. Each interaction is annotated with interaction type (inhibitor, activator, antagonist, etc.), source databases, and supporting publications. Gene druggability categories help identify potentially tractable therapeutic targets.

The database is particularly valuable for identifying existing drugs targeting genes of interest, exploring drug repurposing opportunities, and understanding the druggable genome landscape.

## Key Statistics

| Metric | Value |
|--------|-------|
| Drug-Gene Interactions | 90,000+ |
| Genes | 13,000+ |
| Drugs | 10,000+ |
| Source Databases | 40+ |
| Gene Categories | 40+ |
| Last Update | 2024 (v5) |

## Primary Use Cases

1. **Target Validation** - Assess druggability of genes of interest
2. **Drug Repurposing** - Find existing drugs for new targets
3. **Interaction Lookup** - Query drug-gene relationships
4. **Druggable Genome** - Explore tractable target space
5. **Clinical Translation** - Connect genomic findings to drugs

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Gene Symbol | HGNC | EGFR |
| Concept ID | Normalized | HGNC:3236 |
| Drug Name | Text | Gefitinib |
| Source | Database name | ChEMBL |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://dgidb.org | Search and browse |
| GraphQL API | https://dgidb.org/api/graphql | Flexible queries |
| Bulk Download | https://dgidb.org/downloads | TSV files |
| R Package | rDGIdb | Bioconductor |

## Interaction Types

| Type | Description |
|------|-------------|
| inhibitor | Drug inhibits gene/protein |
| activator | Drug activates gene/protein |
| antagonist | Drug blocks receptor |
| agonist | Drug activates receptor |
| antibody | Antibody-based therapeutic |
| modulator | Drug modulates activity |

## Gene Categories

| Category | Description |
|----------|-------------|
| KINASE | Protein kinase genes |
| DRUGGABLE GENOME | Potentially druggable |
| CLINICALLY ACTIONABLE | Clinical relevance |
| TUMOR SUPPRESSOR | Cancer suppressor genes |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution Required | Yes |

## Related Resources

- [Open Targets](../../_index.md) - Target-disease associations
- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Bioactivity
- [TTD](../ttd/_index.md) - Therapeutic targets

## See Also

- [Schema Documentation](./schema.md)

## References

1. Freshour SL, et al. (2021) "Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts." Nucleic Acids Res. 49(D1):D1144-D1151.
