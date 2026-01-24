---
id: mondo
title: "MONDO Disease Ontology"
type: data-source
category: diseases
subcategory: disease.ontologies
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [ontology, diseases, rare-diseases, cross-references, monarch-initiative]
---

# MONDO Disease Ontology

**Category:** [Diseases & Phenotypes](../../_index.md) > [Disease Ontologies](../_index.md)

## Overview

MONDO (Monarch Disease Ontology) is a unified disease ontology that integrates multiple disease terminologies into a coherent, logic-based structure. Unlike loose cross-references, MONDO provides precise 1:1 equivalence axioms validated by OWL reasoning, enabling safe data propagation across OMIM, Orphanet, EFO, DOID, and NCIt.

Developed by the Monarch Initiative, MONDO addresses the long-standing challenge of disease nomenclature fragmentation by creating a single, authoritative ontology that harmonizes diverse disease classifications. Each MONDO term is carefully mapped to equivalent concepts in source databases, with explicit semantics for equivalence, narrower, and broader relationships.

MONDO is particularly valuable for rare disease research, providing comprehensive coverage of Mendelian conditions with robust links to OMIM and Orphanet. The ontology supports computational disease matching, patient cohort identification, and knowledge graph construction across the biomedical research ecosystem.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Diseases | 25,938 |
| Human Diseases | 22,977 |
| Cancer Classes | 4,728 |
| Mendelian Conditions | 11,639 |
| Database Cross-references | 129,914 |

## Primary Use Cases

1. Disease ID mapping across OMIM, Orphanet, DOID, NCIt, and EFO
2. Rare disease patient matching and cohort identification
3. Knowledge graph construction for translational research
4. Clinical data harmonization across institutions
5. Disease-gene association integration

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| MONDO ID | `MONDO:[0-9]{7}` | MONDO:0005015 |
| OMIM XREF | `OMIM:[0-9]{6}` | OMIM:154700 |
| Orphanet XREF | `Orphanet:[0-9]+` | Orphanet:558 |
| DOID XREF | `DOID:[0-9]+` | DOID:14323 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| MONDO Homepage | https://mondo.monarchinitiative.org/ | Official site |
| GitHub | https://github.com/monarch-initiative/mondo | Source repository |
| OBO Foundry | https://obofoundry.org/ontology/mondo.html | Registry entry |
| BioPortal | https://bioportal.bioontology.org/ontologies/MONDO | Browser |

## Data Formats

| Format | File | Size |
|--------|------|------|
| OWL | mondo.owl | 239.3 MB |
| OBO | mondo.obo | 50.7 MB |
| JSON | mondo.json | 102.7 MB |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required - "Mondo Disease Ontology, Monarch Initiative" |

## See Also

- [Schema Documentation](./schema.md)
- [HPO](../../3.2.phenotype.databases/hpo/_index.md) - Phenotype ontology
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/_index.md) - Rare disease database
