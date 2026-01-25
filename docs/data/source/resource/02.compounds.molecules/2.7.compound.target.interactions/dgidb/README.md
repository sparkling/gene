---
id: dgidb
title: "DGIdb - Drug Gene Interaction Database"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: compound.target.interactions
tags:
  - drug-gene-interactions
  - druggable-genome
  - targets
  - pharmacogenomics
  - graphql
---

# DGIdb - Drug Gene Interaction Database

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

## Limitations

- Aggregated data may contain conflicting interaction types
- Druggability categories are predictions, not experimental validation
- Source database versions may lag behind primary sources
- Gene nomenclature normalized but some ambiguities remain

## Related Resources

- [Open Targets](../../_index.md) - Target-disease associations
- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Bioactivity
- [TTD](../ttd/_index.md) - Therapeutic targets

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## References

1. Freshour SL, et al. (2021) "Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts." Nucleic Acids Res. 49(D1):D1144-D1151.
