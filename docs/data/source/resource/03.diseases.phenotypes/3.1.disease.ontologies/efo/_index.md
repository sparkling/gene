---
id: efo
title: "Experimental Factor Ontology (EFO)"
type: data-source
category: diseases
subcategory: disease.ontologies
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [ontology, traits, diseases, experimental-factors, gwas, ebi]
---

# Experimental Factor Ontology (EFO)

**Category:** [Diseases & Phenotypes](../../_index.md) > [Disease Ontologies](../_index.md)

## Overview

The Experimental Factor Ontology (EFO) provides a systematic description of experimental variables available in EBI databases. Originally developed to describe experimental factors in functional genomics experiments at ArrayExpress, EFO has expanded to become a key ontology for describing diseases, phenotypes, and traits in GWAS Catalog and Open Targets.

EFO integrates terms from multiple domain ontologies including Disease Ontology, Human Phenotype Ontology, and ChEBI, providing a unified framework for annotating experimental data. The ontology is maintained by the European Bioinformatics Institute (EBI) and undergoes regular updates to incorporate new terms and improve mappings.

As a central component of the GWAS Catalog annotation pipeline, EFO enables standardized trait descriptions across thousands of genome-wide association studies. This standardization facilitates meta-analyses and cross-study comparisons while maintaining links to source ontologies through cross-references.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Terms | ~50,000+ |
| Disease Terms | ~20,000 |
| Phenotype Terms | ~10,000 |
| Release Cycle | Monthly |
| Last Update | January 2026 |

## Primary Use Cases

1. Annotation of traits and diseases in GWAS Catalog studies
2. Standardizing experimental factor descriptions in ArrayExpress
3. Disease ontology integration for Open Targets Platform
4. Cross-referencing between disease databases

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| EFO ID | `EFO_[0-9]{7}` | EFO_0000685 |
| Orphanet XREF | `Orphanet:[0-9]+` | Orphanet:558 |
| DOID XREF | `DOID:[0-9]+` | DOID:14323 |
| HP XREF | `HP:[0-9]{7}` | HP:0001250 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| OLS Browser | https://www.ebi.ac.uk/ols4/ontologies/efo | Web interface |
| GitHub | https://github.com/EBISPOT/efo | Source repository |
| OBO Foundry | https://obofoundry.org/ontology/efo.html | Ontology registry |
| PURL | http://www.ebi.ac.uk/efo/efo.owl | Direct download |

## Data Formats

| Format | File | Size |
|--------|------|------|
| OWL | efo.owl | ~150 MB |
| OBO | efo.obo | ~50 MB |
| JSON | efo.json | ~80 MB |

## License

| Aspect | Value |
|--------|-------|
| License | Apache 2.0 |
| Commercial Use | Yes |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [MONDO](../mondo/_index.md) - Related disease ontology
- [HPO](../../3.2.phenotype.databases/hpo/_index.md) - Phenotype ontology
