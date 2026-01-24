---
id: bindingdb
title: "BindingDB - Binding Affinity Database"
type: data-source
category: compounds
subcategory: compound-target-interactions
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [binding-affinity, drug-target, ic50, ki, pharmacology, drug-discovery]
---

# BindingDB - Binding Affinity Database

**Category:** [Compounds & Molecules](../../_index.md) > [Compound-Target Interactions](../_index.md)

## Overview

BindingDB is a public database of measured binding affinities, focusing chiefly on the interactions of proteins considered to be drug targets with small, drug-like molecules. It contains over 2.9 million binding data points for 1.3 million compounds and 9,400 protein targets, making it one of the largest resources for quantitative drug-target interaction data.

The database provides detailed affinity measurements including IC50, Ki, Kd, and EC50 values, along with experimental conditions (pH, temperature), assay descriptions, and literature references. This quantitative data is essential for structure-activity relationship (SAR) analysis, computational model training, and drug candidate prioritization.

BindingDB integrates data from scientific publications, patents, and direct depositions, with extensive cross-references to PubChem, ChEMBL, DrugBank, and UniProt.

## Key Statistics

| Metric | Value |
|--------|-------|
| Binding Measurements | 2,900,000+ |
| Distinct Compounds | 1,300,000+ |
| Protein Targets | ~9,400 |
| IC50 Measurements | 1,800,000+ |
| Ki Measurements | 560,000+ |
| Publications | 40,000+ |
| Patents | 35,000+ |
| Last Update | Weekly |

## Primary Use Cases

1. **Lead Identification** - Find active compounds for target of interest
2. **SAR Analysis** - Explore structure-activity relationships
3. **Model Training** - Develop binding affinity prediction models
4. **Target Deconvolution** - Identify targets for active compounds
5. **Selectivity Profiling** - Compare affinity across target families

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| MonomerID | Integer | 50000001 |
| UniProt ID | Accession | P35354 |
| InChI Key | 27 characters | Standard format |
| PDB ID | 4 characters | 1CX2 |
| ChEMBL ID | CHEMBL + digits | CHEMBL521 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.bindingdb.org | Search and browse |
| REST API | https://bindingdb.org/rest/ | Programmatic access |
| Bulk Download | https://www.bindingdb.org/download | TSV, SDF formats |

## Affinity Types

| Type | Description | Count |
|------|-------------|-------|
| IC50 | Half-maximal inhibitory concentration | 1.8M |
| Ki | Inhibition constant | 560K |
| Kd | Dissociation constant | 100K |
| EC50 | Half-maximal effective concentration | 220K |

## API Endpoints

| Endpoint | Description |
|----------|-------------|
| getLigandsByUniprots | Compounds binding to UniProt targets |
| getLigandsByPDBs | Compounds binding to PDB structures |
| getTargetByCompound | Targets for a given compound |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 3.0 |
| Commercial Use | Yes |
| Attribution Required | Yes |

## Related Resources

- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Bioactivity data
- [GtoPdb](../gtopdb/_index.md) - Pharmacology database
- [TTD](../ttd/_index.md) - Therapeutic targets

## See Also

- [Schema Documentation](./schema.md)

## References

1. Liu T, et al. (2025) "BindingDB in 2024: a comprehensive database of measured binding affinities." Nucleic Acids Res. 53(D1):D1633-D1641.
