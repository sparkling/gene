---
id: gtopdb
title: "GtoPdb - Guide to Pharmacology"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: compound.target.interactions
tags:
  - pharmacology
  - iuphar
  - receptors
  - ion-channels
  - drug-targets
  - curated
---

# GtoPdb - Guide to Pharmacology

## Overview

The Guide to Pharmacology (GtoPdb), developed by IUPHAR/BPS, is an expert-curated database of pharmacological targets and the substances that act on them. It provides authoritative, peer-reviewed information on drug targets including GPCRs, ion channels, nuclear hormone receptors, kinases, and other protein families.

Unlike large-scale aggregated databases, GtoPdb emphasizes expert curation and quality over quantity. Each target entry includes detailed pharmacological information, endogenous ligands, selective tool compounds, and approved drugs. Interaction data includes quantitative affinities (pKi, pIC50, pEC50) with standardized conditions.

GtoPdb serves as the official database of the International Union of Basic and Clinical Pharmacology (IUPHAR), providing nomenclature standards and reference information for the pharmacology community.

## Key Statistics

| Metric | Value |
|--------|-------|
| Targets | ~3,000 |
| Target Families | 800+ |
| Ligands | ~13,000 |
| Approved Drugs | ~2,200 |
| Interactions | 170,000+ |
| References | 85,000+ |
| Last Update | Quarterly |

## Primary Use Cases

1. **Pharmacological Reference** - Authoritative target and ligand information
2. **Drug Discovery** - Identify selective ligands for targets
3. **Target Classification** - IUPHAR-standard nomenclature
4. **Tool Compound Selection** - Find selective probe molecules
5. **Teaching Resource** - Educational pharmacology content

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Target ID | Integer | 290 (BRAF) |
| Ligand ID | Integer | 5085 (vemurafenib) |
| HGNC Symbol | Gene symbol | BRAF |
| ChEMBL ID | Cross-reference | CHEMBL1667 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.guidetopharmacology.org | Search and browse |
| REST API | https://www.guidetopharmacology.org/services | Programmatic access |
| Bulk Download | https://www.guidetopharmacology.org/download.jsp | CSV, PostgreSQL |

## Target Families

| Family | Description |
|--------|-------------|
| GPCRs | G protein-coupled receptors |
| Ion Channels | Voltage-gated, ligand-gated |
| NHRs | Nuclear hormone receptors |
| Enzymes | Kinases, proteases, etc. |
| Transporters | SLC, ABC families |

## Affinity Types (pX notation)

| Type | Description |
|------|-------------|
| pKi | -log10(Ki) inhibitor binding |
| pIC50 | -log10(IC50) inhibition |
| pEC50 | -log10(EC50) activation |
| pKd | -log10(Kd) dissociation |

## Limitations

- Coverage prioritizes established drug targets over novel targets
- Quarterly updates mean recently approved drugs may lag
- Expert curation means smaller scale than aggregated databases
- pX values converted from varied experimental conditions

## Related Resources

- [BindingDB](../bindingdb/_index.md) - Binding affinities
- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Bioactivity
- [DrugBank](../../2.2.pharmaceuticals/drugbank/_index.md) - Drug information

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## References

1. Harding SD, et al. (2024) "The IUPHAR/BPS Guide to PHARMACOLOGY in 2024." Nucleic Acids Res. 52(D1):D1438-D1449.
