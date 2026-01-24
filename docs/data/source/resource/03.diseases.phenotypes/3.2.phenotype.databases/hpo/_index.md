---
id: hpo
title: "Human Phenotype Ontology (HPO)"
type: data-source
category: diseases
subcategory: phenotype.databases
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [ontology, phenotypes, rare-diseases, clinical-features, diagnosis]
---

# Human Phenotype Ontology (HPO)

**Category:** [Diseases & Phenotypes](../../_index.md) > [Phenotype Databases](../_index.md)

## Overview

The Human Phenotype Ontology (HPO) provides a standardized vocabulary of phenotypic abnormalities and clinical features encountered in human disease. It enables precise clinical phenotyping for rare disease diagnosis, genomic interpretation, and translational research.

HPO is organized hierarchically under major branches including Phenotypic abnormality, Mode of inheritance, Clinical modifier, and Clinical course. Each term has a precise definition, synonyms for clinical and lay terminology, and cross-references to SNOMED CT, UMLS, and MeSH. The ontology is maintained by a global consortium of clinical geneticists and biocurators.

Disease-phenotype annotations in HPO connect over 9,500 rare diseases to their characteristic clinical features with evidence codes and frequency information. These annotations power computational phenotype matching tools like Exomiser and LIRICAL, enabling prioritization of candidate disease genes based on patient phenotypes.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Terms | 13,000+ |
| Disease Annotations | 156,000+ |
| Annotated Diseases | ~9,500 |
| PubMed Citations | 9,573 |
| Languages | 10+ |

## Primary Use Cases

1. Clinical phenotyping for rare disease diagnosis
2. Variant prioritization in genomic analysis pipelines
3. Patient cohort matching and similarity scoring
4. Phenopacket data exchange (GA4GH standard)
5. Drug repurposing based on shared phenotypes

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| HPO ID | `HP:[0-9]{7}` | HP:0001250 (Seizure) |
| OMIM Disease | `OMIM:[0-9]{6}` | OMIM:154700 |
| Orphanet Disease | `ORPHA:[0-9]+` | ORPHA:558 |
| UMLS XREF | `UMLS:C[0-9]+` | UMLS:C0036572 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| HPO Website | https://hpo.jax.org/ | Official portal |
| GitHub | https://github.com/obophenotype/human-phenotype-ontology | Source |
| OBO Foundry | https://obofoundry.org/ontology/hp.html | Registry |
| PURL | http://purl.obolibrary.org/obo/hp.obo | Stable download |

## Data Formats

| Format | File | Size |
|--------|------|------|
| OBO | hp.obo | 10.4 MB |
| OWL | hp.owl | 75.6 MB |
| JSON | hp.json | 21.9 MB |
| Annotations | phenotype.hpoa | 35.0 MB |
| Gene Links | genes_to_phenotype.txt | 20.1 MB |

## License

| Aspect | Value |
|--------|-------|
| License | Open Access |
| Commercial Use | Yes |
| Attribution | Recommended |
| URL | https://hpo.jax.org/app/license |

## See Also

- [Schema Documentation](./schema.md)
- [OMIM](../omim/_index.md) - Mendelian disease database
- [MONDO](../../3.1.disease.ontologies/mondo/_index.md) - Disease ontology
