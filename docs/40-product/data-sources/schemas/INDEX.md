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
| [WIKIDATA-SCHEMA.md](./WIKIDATA-SCHEMA.md) | Wikidata | SPARQL endpoint schema with gene/disease/chemical properties (P351, P352, P492, etc.) |
| [UNIPROT-IDMAPPING-SCHEMA.md](./UNIPROT-IDMAPPING-SCHEMA.md) | UniProt ID Mapping | 22-column idmapping_selected.tab format and 286 cross-referenced databases |
| [MONDO-SCHEMA.md](./MONDO-SCHEMA.md) | MONDO Disease Ontology | OBO/OWL format, 26K diseases, OMIM/Orphanet equivalence mappings |
| [HPO-SCHEMA.md](./HPO-SCHEMA.md) | Human Phenotype Ontology | OBO format, phenotype.hpoa annotation format, 13K+ terms |

### Pathway Databases

| File | Database | Description |
|------|----------|-------------|
| [REACTOME-SCHEMA.md](./REACTOME-SCHEMA.md) | Reactome | Neo4j graph database schema with node types, relationships, and API endpoints |
| [WIKIPATHWAYS-GPML-SCHEMA.md](./WIKIPATHWAYS-GPML-SCHEMA.md) | WikiPathways | GPML XML schema documentation with all element types and attributes |

### Disease & Variation Databases

| File | Database | Description |
|------|----------|-------------|
| [DISGENET-SCHEMA.md](./DISGENET-SCHEMA.md) | DisGeNET | Gene-disease association schema with score metrics (GDA, EI, DSI, DPI) |
| [CLINVAR-SCHEMA.md](./CLINVAR-SCHEMA.md) | ClinVar | Clinical variant interpretation schema |
| [DBSNP-SCHEMA.md](./DBSNP-SCHEMA.md) | dbSNP | SNP variation schema |
| [GWAS-CATALOG-SCHEMA.md](./GWAS-CATALOG-SCHEMA.md) | GWAS Catalog | GWAS association schema |
| [ORPHANET-ORDO-SCHEMA.md](./ORPHANET-ORDO-SCHEMA.md) | Orphanet/ORDO | Rare disease ontology schema |

### Drug & Chemical Databases

| File | Database | Description |
|------|----------|-------------|
| [CHEMBL-SCHEMA.md](./CHEMBL-SCHEMA.md) | ChEMBL | Bioactivity and drug schema |
| [PHARMGKB-SCHEMA.md](./PHARMGKB-SCHEMA.md) | PharmGKB | Pharmacogenomics schema |
| [COCONUT-SCHEMA.md](./COCONUT-SCHEMA.md) | COCONUT | Natural products schema |

### Other Databases

| File | Database | Description |
|------|----------|-------------|
| [SAMPLE-DATA.md](./SAMPLE-DATA.md) | All | Actual sample data retrieved from APIs |
| [UNIFIED-SCHEMA-ANALYSIS.md](./UNIFIED-SCHEMA-ANALYSIS.md) | Cross-database | Unified schema analysis |

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

- [43-41-PATHWAYS-PRIMARY.md](../43-41-PATHWAYS-PRIMARY.md) - Primary pathway database overview
- [43-43-PATHWAYS-DISEASE.md](../43-43-PATHWAYS-DISEASE.md) - Disease pathway databases
- [43-86-INTEGRATION-XREFS.md](../43-86-INTEGRATION-XREFS.md) - Cross-reference mapping strategies
