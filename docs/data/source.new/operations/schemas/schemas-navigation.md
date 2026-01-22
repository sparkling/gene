---
id: schemas-navigation
title: "Database Schemas Index"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, index, navigation, databases, pathways, diseases, ontologies]
---

**Parent:** [Schema Documentation](./_index.md)

# Database Schemas Index

**Document ID:** SCHEMAS-INDEX
**Last Updated:** January 2026
**Status:** Complete

---

## Overview

This directory contains actual schema documentation, data models, and sample data from biomedical databases including pathways, diseases, ontologies, and cross-reference systems.

---

## Contents

### Integration & Ontology Databases

| File | Database | Description |
|------|----------|-------------|
| [wikidata-schema.md](./wikidata-schema.md) | Wikidata | SPARQL endpoint schema with gene/disease/chemical properties (P351, P352, P492, etc.) |
| [uniprot-idmapping-schema.md](./uniprot-idmapping-schema.md) | UniProt ID Mapping | 22-column idmapping_selected.tab format and 286 cross-referenced databases |
| [mondo-schema.md](./mondo-schema.md) | MONDO Disease Ontology | OBO/OWL format, 26K diseases, OMIM/Orphanet equivalence mappings |
| [hpo-schema.md](./hpo-schema.md) | Human Phenotype Ontology | OBO format, phenotype.hpoa annotation format, 13K+ terms |

### Pathway Databases

| File | Database | Description |
|------|----------|-------------|
| [reactome-schema.md](./reactome-schema.md) | Reactome | Neo4j graph database schema with node types, relationships, and API endpoints |
| [wikipathways-gpml-schema.md](./wikipathways-gpml-schema.md) | WikiPathways | GPML XML schema documentation with all element types and attributes |

### Disease & Variation Databases

| File | Database | Description |
|------|----------|-------------|
| [disgenet-schema.md](./disgenet-schema.md) | DisGeNET | Gene-disease association schema with score metrics (GDA, EI, DSI, DPI) |
| [clinvar-schema.md](./clinvar-schema.md) | ClinVar | Clinical variant interpretation schema |
| [dbsnp-schema.md](./dbsnp-schema.md) | dbSNP | SNP variation schema |
| [gwas-catalog-schema.md](./gwas-catalog-schema.md) | GWAS Catalog | GWAS association schema |
| [orphanet-ordo-schema.md](./orphanet-ordo-schema.md) | Orphanet/ORDO | Rare disease ontology schema |

### Drug & Chemical Databases

| File | Database | Description |
|------|----------|-------------|
| [chembl-schema.md](./chembl-schema.md) | ChEMBL | Bioactivity and drug schema |
| [pharmgkb-schema.md](./pharmgkb-schema.md) | PharmGKB | Pharmacogenomics schema |
| [coconut-schema.md](./coconut-schema.md) | COCONUT | Natural products schema |

### Other Databases

| File | Database | Description |
|------|----------|-------------|
| [sample-data.md](./sample-data.md) | All | Actual sample data retrieved from APIs |
| [unified-schema-analysis.md](./unified-schema-analysis.md) | Cross-database | Unified schema analysis |

---

## Database Summary

### Reactome
- **URL:** https://reactome.org
- **License:** CC BY 4.0
- **Format:** Neo4j Graph Database
- **Statistics:** 2,712 human pathways, 13,872 reactions, 11,196 proteins, 1,925 small molecules
- **Key Features:**
  - Event hierarchy (TopLevelPathway > Pathway > Reaction)
  - PhysicalEntity classes (Protein, Complex, SimpleEntity)
  - GO annotations and literature references
  - 16 supported species via orthology

### WikiPathways
- **URL:** https://www.wikipathways.org
- **License:** CC0 1.0 (Public Domain)
- **Format:** GPML (XML-based)
- **Statistics:** 3,100+ pathways, 48 organisms, 955+ human pathways
- **Key Features:**
  - DataNode types (GeneProduct, Metabolite, Protein, Pathway)
  - MIM notation for interaction types
  - BridgeDb identifier mapping
  - Community curation model

### DisGeNET
- **URL:** https://www.disgenet.org
- **License:** CC BY-NC-SA 4.0 (Academic)
- **Format:** TSV, SQLite, REST API
- **Statistics:** 628,685 GDAs, 210,498 VDAs, 17,549 genes, 24,166 diseases
- **Key Features:**
  - GDA Score (confidence metric 0-1)
  - Evidence Index (publication consistency)
  - Disease Specificity Index (DSI)
  - Disease Pleiotropy Index (DPI)

---

## Cross-Reference Identifiers

### Gene/Protein Identifiers

| Database | Primary ID | Example |
|----------|------------|---------|
| Reactome | UniProt | Q9UBM7 |
| WikiPathways | Entrez Gene, Ensembl | 3156, ENSG00000113161 |
| DisGeNET | NCBI Gene | 3156 |

### Compound/Metabolite Identifiers

| Database | Primary ID | Example |
|----------|------------|---------|
| Reactome | ChEBI | CHEBI:16113 |
| WikiPathways | ChEBI, HMDB, CAS | CHEBI:16113, HMDB0000067, 57-88-5 |

### Disease Identifiers

| Database | Primary ID | Example |
|----------|------------|---------|
| DisGeNET | UMLS CUI | C0020443 |
| Cross-refs | DOID, HPO, MeSH | DOID:1168, HP:0003124 |

---

## API Access

### Reactome Content Service
```
Base URL: https://reactome.org/ContentService
Rate Limit: 100 requests/minute
Auth: None required
Format: JSON, XML
```

### WikiPathways Web Service
```
Base URL: https://webservice.wikipathways.org
Rate Limit: Reasonable use
Auth: None required
Format: JSON, XML
```

### DisGeNET API
```
Base URL: https://www.disgenet.org/api
Rate Limit: Account-based
Auth: API key required
Format: JSON, TSV, XML
```

---

## Integration Recommendations

### Primary Identifiers for Harmonization

| Entity Type | Recommended ID | Mapping Resources |
|-------------|----------------|-------------------|
| Genes | Ensembl Gene ID | BridgeDb, UniProt ID Mapping |
| Proteins | UniProt Accession | UniProt ID Mapping |
| Metabolites | ChEBI ID | ChEBI Ontology |
| Diseases | UMLS CUI | UMLS Metathesaurus |

### Data Flow

```
Gene Variant -> DisGeNET -> Disease
     |              |           |
     v              v           v
  Ensembl       UMLS CUI     DOID/HPO
     |              |           |
     v              v           v
Reactome/WikiPathways -> Pathway -> Compound (ChEBI)
```

---

## Related Documentation

- [../pathways/primary.md](../pathways/primary.md) - Primary pathway database overview
- [../pathways/disease.md](../pathways/disease.md) - Disease pathway databases
- [../integration/xrefs.md](../integration/xrefs.md) - Cross-reference mapping strategies

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `schema` | Database structure definition describing tables and fields | Reactome schema |
| `cross-reference` | Identifier mapping between different databases | UniProt to Entrez |
| `API` | Application Programming Interface for programmatic access | REST API |
| `rate_limit` | Maximum number of API requests allowed per time period | 100 requests/minute |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Reactome | Curated pathway database with Neo4j graph structure | Pathways |
| WikiPathways | Community-curated pathway resource using GPML format | Pathways |
| DisGeNET | Gene-disease association database with scoring metrics | Disease associations |
| BridgeDb | Cross-reference mapping service for identifier translation | Integration |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GDA | Gene-Disease Association | DisGeNET metric |
| GPML | Graphical Pathway Markup Language | WikiPathways format |
| MIM | Molecular Interaction Map | WikiPathways notation |
| FHIR | Fast Healthcare Interoperability Resources | Data exchange standard |
