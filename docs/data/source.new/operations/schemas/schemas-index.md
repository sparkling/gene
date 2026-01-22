---
id: schemas-schemas-index
title: "Database Schemas Index"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, index, three-worlds, genetics, traditional-medicine, nutrition, integration]
---

**Parent:** [Schema Documentation](./_index.md)

# Database Schemas Index

**Document ID:** SCHEMAS-INDEX
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 3.0

---

## Executive Summary

This directory contains **42 schema documentation files** covering databases across the THREE WORLDS framework (Modern Genetics, Traditional Medicine, Nutritional Science) plus cross-cutting resources. The research phase involved 8 parallel agents fetching actual API specifications, sample data, and statistics.

### Key Metrics

| Metric | Value |
|--------|-------|
| **Schema Documents** | 42 files |
| **Databases Documented** | 38 unique databases |
| **Total Records Accessible** | 1B+ variants, 3.6M compounds, 80K diseases |
| **CC0/Public Domain** | 12 databases (MVP-ready) |
| **CC BY/CC BY-SA** | 16 databases (commercial-friendly) |
| **Non-Commercial Only** | 5 databases |
| **Academic/Special License** | 5 databases |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "schema-001" |
| `name` | string | Entity name | "Database Name" |
| `type` | string | Record type | "schema" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## THREE WORLDS Coverage

### WORLD 1: Modern Genetics (11 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [dbsnp-schema.md](./dbsnp-schema.md) | dbSNP/ALFA R4 | 900M+ variants, 12 populations | Public Domain | Critical |
| [clinvar-schema.md](./clinvar-schema.md) | ClinVar | 2.5M+ clinical variants | Public Domain | Critical |
| [pharmgkb-schema.md](./pharmgkb-schema.md) | PharmGKB/ClinPGx | 33 CPIC guidelines | CC BY-SA 4.0 | Critical |
| [gwas-catalog-schema.md](./gwas-catalog-schema.md) | GWAS Catalog | 186K studies, 1M+ associations | CC0 | High |
| [gnomad-schema.md](./gnomad-schema.md) | gnomAD | 76K+ genomes, 730K exomes | Open Access | Critical |
| [alphamissense-schema.md](./alphamissense-schema.md) | AlphaMissense | 71M missense predictions | CC BY 4.0 | High |
| [dbnsfp-schema.md](./dbnsfp-schema.md) | dbNSFP | 84M+ variants, 46 prediction scores | Academic/Commercial | High |
| [dbvar-schema.md](./dbvar-schema.md) | dbVar | 6.5M+ structural variants | Public Domain | High |
| [spliceai-schema.md](./spliceai-schema.md) | SpliceAI | Pre-computed splice predictions | Illumina License | Medium |
| [encode-schema.md](./encode-schema.md) | ENCODE | 1M+ functional elements | Open Access | High |
| [gene-ontology-schema.md](./gene-ontology-schema.md) | Gene Ontology | 44K terms, 7M annotations | CC BY 4.0 | Critical |

### WORLD 2: Traditional Medicine (5 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [kampodb-schema.md](./kampodb-schema.md) | KampoDB (Kampo) | 481 formulas, 3K compounds, 63K targets | **CC BY-SA 4.0** | High |
| [dr-dukes-schema.md](./dr-dukes-schema.md) | Dr. Duke's (Western Herbal) | 50K entries, 1,900 activities | **CC0** | High |
| [batman-tcm-schema.md](./batman-tcm-schema.md) | BATMAN-TCM 2.0 (TCM) | 54K formulas, 39K compounds | CC BY-NC 4.0 | Medium |
| [imppat-schema.md](./imppat-schema.md) | IMPPAT 2.0 (Ayurveda) | 4K plants, 18K phytochemicals | CC BY-NC 4.0 | Medium |
| [lotus-schema.md](./lotus-schema.md) | LOTUS | 750K+ NP-organism pairs | **CC0** | High |

### WORLD 3: Nutritional Science (5 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [usda-fooddata-central.md](./usda-fooddata-central.md) | USDA FoodData Central | 20,900+ foods | **CC0** | Critical |
| [dsld-nih.md](./dsld-nih.md) | DSLD (NIH) | 200K+ supplement labels | Public Domain | High |
| [foodb.md](./foodb.md) | FooDB | 71K compounds, 778 foods | Free (citation) | High |
| [open-food-facts-schema.md](./open-food-facts-schema.md) | Open Food Facts | 3M+ food products | **ODbL** | High |
| [hmp-schema.md](./hmp-schema.md) | Human Microbiome Project | 5K+ metagenomes | CC BY 4.0 | Medium |

---

## Cross-Cutting Resources

### Natural Products & Interventions (5 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [chembl-schema.md](./chembl-schema.md) | ChEMBL 36 | 2.9M compounds, 24.3M activities | CC BY-SA 3.0 | Critical |
| [coconut-schema.md](./coconut-schema.md) | COCONUT | 716K natural products | **CC0** | Critical |
| [pubchem-schema.md](./pubchem-schema.md) | PubChem | 115M+ compounds | **Public Domain** | Critical |
| [chebi-schema.md](./chebi-schema.md) | ChEBI | 60K+ chemical entities | CC BY 4.0 | Critical |
| [dr-dukes-schema.md](./dr-dukes-schema.md) | Dr. Duke's | 50K ethnobotanical entries | **CC0** | High |

### Pathways & Networks (6 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [reactome-schema.md](./reactome-schema.md) | Reactome | 2,712 pathways, 11K proteins | CC BY 4.0 | Critical |
| [wikipathways-gpml-schema.md](./wikipathways-gpml-schema.md) | WikiPathways | 3,100+ pathways | CC BY 4.0 | High |
| [kgml-schema.md](./kgml-schema.md) | KEGG KGML | Reference format | Academic License | Reference |
| [disgenet-schema.md](./disgenet-schema.md) | DisGeNET | 628K gene-disease associations | CC BY-NC-SA 4.0 | High |
| [string-schema.md](./string-schema.md) | STRING | 68M+ protein interactions | CC BY 4.0 | Critical |
| [intact-schema.md](./intact-schema.md) | IntAct | 1.1M+ molecular interactions | CC BY 4.0 | High |

### Diseases & Phenotypes (3 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [mondo-schema.md](./mondo-schema.md) | MONDO | 26K diseases, 130K cross-refs | CC BY 4.0 | Critical |
| [hpo-schema.md](./hpo-schema.md) | HPO | 13K phenotypes, 156K annotations | Open Access | Critical |
| [orphanet-ordo-schema.md](./orphanet-ordo-schema.md) | Orphanet ORDO | 6.5K rare diseases | CC BY 4.0 | High |

### Integration & Cross-Reference (2 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [wikidata-schema.md](./wikidata-schema.md) | Wikidata | 59K genes, 200K diseases | **CC0** | Critical |
| [uniprot-idmapping-schema.md](./uniprot-idmapping-schema.md) | UniProt ID Mapping | **286 database cross-refs** | CC BY 4.0 | Critical |

---

## Hub Identifier Strategy

Based on cross-database analysis, these identifiers should serve as canonical hubs:

| Entity Type | Primary Hub ID | Rationale | Coverage |
|-------------|----------------|-----------|----------|
| **Gene** | NCBI Gene ID (Entrez) | Universal, links to dbSNP/ClinVar/DisGeNET | 100% human genes |
| **Variant** | rsID + SPDI | rsID for SNVs, SPDI for normalization | 900M+ variants |
| **Compound** | **InChIKey** | Structure-based, 27-char standard | All chemistry DBs |
| **Disease** | **MONDO ID** | Best cross-ref (OMIM, Orphanet, DOID) | 130K mappings |
| **Pathway** | Reactome ID | Highest curation quality | 2.7K pathways |
| **Phenotype** | HPO ID | GA4GH Phenopackets standard | 13K terms |
| **Protein** | UniProt ID | Universal protein identifier | 286 DB mappings |

---

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 19 |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | Markdown |
| Alternative | YAML frontmatter |
| Encoding | UTF-8 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| Schema Docs | Git | See repository |
| Documentation | HTTP | See main database |

**Access Requirements:** Open access

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| Schema Docs | CC BY 4.0 | Yes |

---

## Sample Data

### Example Record
```json
{
  "id": "schema-chembl",
  "title": "ChEMBL Schema Documentation",
  "database": "ChEMBL",
  "records": "2.9M compounds"
}
```

### Sample Query Result
| id | title | database | records |
|----|-------|----------|---------|
| schema-chembl | ChEMBL Schema Documentation | ChEMBL | 2.9M compounds |
| schema-reactome | Reactome Schema Documentation | Reactome | 2,712 pathways |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `THREE_WORLDS` | Framework organizing data into Genetics, Traditional Medicine, Nutrition | WORLD 1, WORLD 2, WORLD 3 |
| `hub_identifier` | Canonical identifier chosen for cross-database integration | InChIKey, MONDO ID |
| `schema_document` | Markdown file documenting a database's structure and fields | chembl-schema.md |
| `license` | Terms governing database usage and redistribution | CC0, CC BY-SA 4.0 |
| `priority` | Implementation order based on value and accessibility | Critical, High, Medium |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| WORLD_1 | Modern Genetics domain (variants, genes, pharmacogenomics) | dbSNP, ClinVar, gnomAD |
| WORLD_2 | Traditional Medicine domain (TCM, Kampo, Ayurveda) | BATMAN-TCM, KampoDB, IMPPAT |
| WORLD_3 | Nutritional Science domain (foods, nutrients, microbiome) | USDA FoodData, FooDB |
| CC0 | Creative Commons Zero public domain dedication | Open license |
| CC_BY | Creative Commons Attribution license | Attribution required |
| CC_BY_NC | Creative Commons Non-Commercial license | No commercial use |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GWAS | Genome-Wide Association Study | Variant-trait links |
| gnomAD | Genome Aggregation Database | Population frequencies |
| CPIC | Clinical Pharmacogenetics Implementation Consortium | Drug-gene guidelines |
| NP | Natural Product | Compound from organisms |
| TCM | Traditional Chinese Medicine | Medical system |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathways |
| KGML | KEGG Markup Language | Pathway format |
| ODbL | Open Database License | Open Food Facts |
| HMP | Human Microbiome Project | Microbiome data |
| ENCODE | Encyclopedia of DNA Elements | Functional genomics |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial index with 3 NP schemas |
| 2.0 | January 2026 | Data Engineering | Complete overhaul: 27 schemas, THREE WORLDS coverage |
| **3.0** | January 2026 | Data Engineering | **Added 14 new schemas**: Gene Ontology, ChEBI, AlphaMissense, dbNSFP, gnomAD, dbVar, STRING, IntAct, PubChem, HMP, Open Food Facts, SpliceAI, ENCODE, LOTUS. Updated coverage to 41 files. |
| **3.1** | January 2026 | Data Engineering | **Added KGML-SCHEMA.md**. Updated coverage to 42 files. |
