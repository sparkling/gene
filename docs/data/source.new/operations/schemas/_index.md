---
title: "Schemas"
parent: ../_index.md
last_updated: 2026-01-22
status: draft
---

# Schemas

This section contains database schema documentation and specifications for all data sources in the Gene Platform.

## Contents

The following 47 schema files document the database structures:

| Schema File | Description |
|-------------|-------------|
| [activities-schema](activities-schema.md) | Activity and compound interaction schemas |
| [activity-props-schema](activity-props-schema.md) | Activity properties and measurements |
| [activity-smid-schema](activity-smid-schema.md) | Activity source molecule identifiers |
| [activity-supp-schema](activity-supp-schema.md) | Supplementary activity data |
| [activity-supp-map-schema](activity-supp-map-schema.md) | Activity supplement mappings |
| [assay-class-schema](assay-class-schema.md) | Assay classification taxonomy |
| [assay-params-schema](assay-params-schema.md) | Assay parameter definitions |
| [assays-schema](assays-schema.md) | Assay metadata and protocols |
| [atc-class-schema](atc-class-schema.md) | Anatomical Therapeutic Chemical classification |
| [binding-sites-schema](binding-sites-schema.md) | Protein binding site annotations |
| [bio-component-seq-schema](bio-component-seq-schema.md) | Biological component sequences |
| [biotherapeutic-comps-schema](biotherapeutic-comps-schema.md) | Biotherapeutic component definitions |
| [biotherapeutics-schema](biotherapeutics-schema.md) | Biotherapeutic drug information |
| [cell-dictionary-schema](cell-dictionary-schema.md) | Cell line reference data |
| [chembl-id-lookup-schema](chembl-id-lookup-schema.md) | ChEMBL identifier cross-references |
| [component-class-schema](component-class-schema.md) | Component classification system |
| [component-domains-schema](component-domains-schema.md) | Protein domain annotations |
| [component-go-schema](component-go-schema.md) | Gene Ontology term associations |
| [component-seqs-schema](component-seqs-schema.md) | Component sequence data |
| [component-synonyms-schema](component-synonyms-schema.md) | Alternative names for components |
| [compound-props-schema](compound-props-schema.md) | Chemical compound properties |
| [compound-records-schema](compound-records-schema.md) | Compound record metadata |
| [compound-struct-alert-schema](compound-struct-alert-schema.md) | Structural alert annotations |
| [compound-structs-schema](compound-structs-schema.md) | Chemical structure data |
| [confidence-score-lookup-schema](confidence-score-lookup-schema.md) | Confidence scoring system |
| [curation-lookup-schema](curation-lookup-schema.md) | Curation status definitions |
| [data-validity-lookup-schema](data-validity-lookup-schema.md) | Data validation rules |
| [defined-daily-dose-schema](defined-daily-dose-schema.md) | WHO defined daily dose data |
| [docs-schema](docs-schema.md) | Document and literature references |
| [domains-schema](domains-schema.md) | Protein domain definitions |
| [drug-indication-schema](drug-indication-schema.md) | Drug indication associations |
| [drug-mechanism-schema](drug-mechanism-schema.md) | Drug mechanism of action data |
| [drug-warning-schema](drug-warning-schema.md) | Drug safety warnings and alerts |
| [formulations-schema](formulations-schema.md) | Drug formulation specifications |
| [frac-class-schema](frac-class-schema.md) | Fungicide Resistance Action Committee classification |
| [go-class-schema](go-class-schema.md) | Gene Ontology classification hierarchy |
| [hrac-class-schema](hrac-class-schema.md) | Herbicide Resistance Action Committee classification |
| [indication-refs-schema](indication-refs-schema.md) | Indication reference data |
| [irac-class-schema](irac-class-schema.md) | Insecticide Resistance Action Committee classification |
| [ligand-eff-schema](ligand-eff-schema.md) | Ligand efficiency metrics |
| [mechanism-refs-schema](mechanism-refs-schema.md) | Mechanism reference citations |
| [metabolism-schema](metabolism-schema.md) | Drug metabolism pathways |
| [metabolism-refs-schema](metabolism-refs-schema.md) | Metabolism reference data |
| [molecule-atc-schema](molecule-atc-schema.md) | Molecule ATC code mappings |
| [molecule-dict-schema](molecule-dict-schema.md) | Molecule dictionary and identifiers |
| [molecule-frac-schema](molecule-frac-schema.md) | Molecule FRAC classification |
| [molecule-hierarchy-schema](molecule-hierarchy-schema.md) | Chemical hierarchy relationships |

## Overview

The schema documentation provides:

- **Table Definitions**: Complete field specifications for all database tables
- **Data Types**: Field types, constraints, and validation rules
- **Relationships**: Foreign keys and table relationships
- **Indexes**: Performance optimization indexes
- **Usage Examples**: Common query patterns and examples

## Schema Categories

### Activity & Assay Schemas
Core schemas for biological activity measurements and assay protocols.

### Chemical Structure Schemas
Schemas defining chemical structures, properties, and classifications.

### Biological Component Schemas
Protein, gene, and biological entity definitions and annotations.

### Drug & Therapeutic Schemas
Drug formulations, mechanisms, indications, and safety information.

### Classification Schemas
Taxonomic and classification system hierarchies (ATC, GO, FRAC, HRAC, IRAC).

### Reference & Metadata Schemas
Literature references, documentation, and data provenance.

## Navigation

- [Parent: Operations](../_index.md)
- [Downloads](../downloads/_index.md)
- [Integration](../integration/_index.md)
- [Governance](../governance/_index.md)
