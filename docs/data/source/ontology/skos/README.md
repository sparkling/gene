---
id: ontology-skos
title: "SKOS Vocabularies"
type: ontology
parent: ../_index.md
last_updated: 2026-01-25
status: active
tags: [skos, vocabulary, taxonomy, thesaurus, concepts]
---

# SKOS Vocabularies

**Parent:** [Ontology](../_index.md)

SKOS (Simple Knowledge Organization System) concept schemes for taxonomies, controlled vocabularies, and thesauri.

## Purpose

- Define hierarchical concept schemes (broader/narrower)
- Manage preferred and alternative labels
- Support multilingual terminology
- Enable faceted classification

## Files

| File | Description |
|------|-------------|
| [datasource-taxonomy.ttl](./datasource-taxonomy.ttl) | Complete data source taxonomy with 6 concept schemes |

## Concept Schemes

### DataSourceTaxonomy

Main 3-level hierarchy: Category > Subcategory > Source

| Level | Count | Example |
|-------|-------|---------|
| Categories | 9 | GeneticsGenomics, CompoundsMolecules, DiseasesPhenotypes |
| Subcategories | ~45 | VariantRepositories, NaturalProducts, DiseaseOntologies |

### Auxiliary Schemes

| Scheme | Concepts | Purpose |
|--------|----------|---------|
| TierScheme | Tier1, Tier2, Tier3 | Quality/priority classification |
| StatusScheme | Draft, Active, Deprecated | Documentation lifecycle |
| AccessMethodScheme | REST API, SPARQL, FTP, S3, etc. | Access method types |
| FormatScheme | JSON, XML, CSV, VCF, FASTA, etc. | Data formats |
| LicenseTypeScheme | Public Domain, CC-BY, Academic Only, etc. | License categories |
| UpdateFrequencyScheme | Daily, Weekly, Monthly, etc. | Update schedules |

## Namespace

```turtle
@prefix tax: <https://gene.ai/taxonomy/> .
```

## Example Query

```sparql
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX tax: <https://gene.ai/taxonomy/>

SELECT ?subcategory ?label
WHERE {
    ?subcategory skos:broader tax:GeneticsGenomics ;
                 skos:prefLabel ?label .
}
```

## Alignment with Folder Structure

| SKOS Concept | Folder Path |
|--------------|-------------|
| tax:GeneticsGenomics | 01.genetics.genomics |
| tax:VariantRepositories | 01.genetics.genomics/1.1.variant.repositories |
| tax:CompoundsMolecules | 02.compounds.molecules |
