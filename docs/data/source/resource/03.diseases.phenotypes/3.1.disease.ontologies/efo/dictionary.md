# EFO - Data Dictionary

## Overview

This data dictionary documents the schema for EFO (Experimental Factor Ontology).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | efo |
| **Name** | EFO |
| **Parent** | 3.1.disease.ontologies |
| **Total Fields** | 15+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### EFO Term

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | EFO identifier | EFO_0000685 |
| label | string | 1:1 | Yes | Term name | rheumatoid arthritis |
| definition | string | 1:1 | No | Textual definition | A chronic systemic... |
| synonyms | array | 1:N | No | Alternative names | [RA, rheumatoid disease] |
| xrefs | array | 1:N | No | Cross-references | [DOID:7148, Orphanet:284] |
| parents | array | 1:N | Yes | Parent term relationships | EFO_0000540 |
| obsolete | boolean | 1:1 | No | Deprecation status | false |
| replaced_by | string | 1:1 | No | Replacement term if obsolete | EFO_0001234 |

### Term Metadata

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| created_by | string | 1:1 | No | Curator ID | EFO curators |
| creation_date | datetime | 1:1 | No | When created | 2015-03-20 |
| source_ontology | string | 1:1 | No | Imported from | HP, DOID |
| subset | array | 1:N | No | Subset membership | efo_slim |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| EFO ID | EFO_[0-9]{7} | EFO_0000685 | Primary identifier |
| Orphanet | Orphanet:[0-9]+ | Orphanet:558 | Rare disease cross-ref |
| DOID | DOID:[0-9]+ | DOID:14323 | Disease ontology |
| HP | HP:[0-9]{7} | HP:0001250 | Phenotype ontology |
| MONDO | MONDO:[0-9]{7} | MONDO:0005015 | Disease ontology |
| MeSH | MESH:D[0-9]+ | MESH:D012640 | Medical subjects |
| NCIt | NCIT:C[0-9]+ | NCIT:C34807 | NCI thesaurus |

---

## Enumerations

### Term Categories

| Category | Description | Examples |
|----------|-------------|----------|
| Disease | Human disease terms | ~20,000 terms |
| Phenotype | Observable characteristics | ~10,000 terms |
| Measurement | Experimental measurements | Various |
| Material | Biological materials | Cell types, tissues |
| Process | Biological processes | Various |
| Anatomical entity | Body parts | Organs, tissues |

### Data Sources

| Source | Description |
|--------|-------------|
| GWAS Catalog | Trait annotations |
| Open Targets | Disease-target associations |
| Expression Atlas | Expression studies |
| ArrayExpress | Microarray data |

### Subset Categories

| Subset | Description |
|--------|-------------|
| efo_slim | Core slim subset |
| gwas | GWAS Catalog traits |
| disease | Disease-specific terms |

### Relationship Types

| Relationship | Description |
|--------------|-------------|
| is_a | Subclass relationship |
| part_of | Part-whole relationship |
| develops_from | Developmental relationship |
| located_in | Anatomical location |

---

## Entity Relationships

### Term to Parents
- **Cardinality:** N:M
- **Description:** Hierarchical is_a relationships
- **Key Fields:** id, parents

### Term to Cross-References
- **Cardinality:** 1:N
- **Description:** Mappings to external ontologies
- **Key Fields:** id, xrefs

### Term to Synonyms
- **Cardinality:** 1:N
- **Description:** Alternative term names
- **Key Fields:** id, synonyms

### Term to Source Ontology
- **Cardinality:** N:1
- **Description:** Terms imported from other ontologies
- **Key Fields:** id, source_ontology

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| EFO | Experimental Factor Ontology | Database name |
| EBI | European Bioinformatics Institute | Host institution |
| OLS | Ontology Lookup Service | EBI browser |
| GWAS | Genome-Wide Association Study | Key use case |
| DOID | Disease Ontology Identifier | Cross-reference |
| NCIt | NCI Thesaurus | Cancer terms |
| HP | Human Phenotype | Phenotype terms |
| OWL | Web Ontology Language | Ontology format |
| OBO | Open Biological Ontologies | Ontology format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| DOID | DOID ID | Disease ontology |
| HPO | HP ID | Phenotype terms |
| MONDO | MONDO ID | Disease ontology |
| Orphanet | ORPHA ID | Rare diseases |
| NCIt | NCIt code | Cancer terms |
| MeSH | MeSH ID | Medical vocabulary |
| SNOMED CT | SCTID | Clinical terms |
| ChEBI | ChEBI ID | Chemical entities |

---

## Data Quality Notes

1. **Term Coverage:** ~50,000+ total terms
2. **Disease Terms:** ~20,000 disease concepts
3. **Integration:** Imports from HPO, DOID, NCIt, ChEBI
4. **GWAS Catalog:** Primary trait vocabulary
5. **Open Targets:** Disease annotation standard
6. **Monthly Updates:** Regular releases
7. **OLS Browser:** https://www.ebi.ac.uk/ols4/ontologies/efo
8. **Apache 2.0:** Open source license

