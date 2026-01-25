# Gene Ontology - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gene.ontology |
| **Name** | Gene Ontology |
| **Total Fields** | 28 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| goid | String | Yes | GO term identifier | GO:0008152 |
| label | String | Yes | GO term name | metabolic process |
| namespace | Enum | Yes | GO ontology aspect | biological_process |
| aspect | Enum | No | Single-letter aspect code | P, F, C |
| definition | String | No | Full term definition | The chemical reactions... |
| is_obsolete | Boolean | No | Whether term is obsolete | false |
| synonyms | Array | No | Alternative names | ["metabolism"] |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| GO ID | GO:####### | Gene Ontology term ID | GO:0008152 |
| Alt ID | GO:####### | Alternative/merged GO IDs | GO:0044237 |
| UniProt | Alphanumeric | Protein accession | P04637 |
| Taxon | taxon:##### | NCBI Taxonomy ID | taxon:9606 |

---

## Enumerations

### Namespaces (Aspects)

| Namespace | Code | Description |
|-----------|------|-------------|
| biological_process | P | Biological objectives |
| molecular_function | F | Biochemical activities |
| cellular_component | C | Subcellular locations |

### Evidence Codes

| Code | Category | Description |
|------|----------|-------------|
| EXP | Experimental | Inferred from experiment |
| IDA | Experimental | Inferred from direct assay |
| IPI | Experimental | Inferred from physical interaction |
| IMP | Experimental | Inferred from mutant phenotype |
| IGI | Experimental | Inferred from genetic interaction |
| IEP | Experimental | Inferred from expression pattern |
| HTP | High-throughput | High-throughput experiment |
| ISS | Computational | Inferred from sequence similarity |
| ISO | Computational | Inferred from sequence orthology |
| IBA | Phylogenetic | Inferred from biological ancestry |
| TAS | Literature | Traceable author statement |
| IEA | Electronic | Inferred from electronic annotation |
| ND | No data | No biological data available |

### Relationship Types

| Type | Description |
|------|-------------|
| is_a | Subclass relationship |
| part_of | Parthood relationship |
| has_part | Inverse of part_of |
| regulates | General regulation |
| positively_regulates | Positive regulation |
| negatively_regulates | Negative regulation |

### Annotation Qualifiers

| Qualifier | Description |
|-----------|-------------|
| enables | MF annotation qualifier |
| involved_in | BP annotation qualifier |
| located_in | CC annotation qualifier |
| part_of | Component part of |
| NOT | Negation qualifier |
| contributes_to | Partial contribution |

---

## GAF 2.2 Annotation Fields

| Field | Description | Examples |
|-------|-------------|----------|
| DB | Contributing database | UniProtKB |
| DB_Object_ID | Protein/gene ID | P04637 |
| DB_Object_Symbol | Gene symbol | TP53 |
| Qualifier | Annotation qualifier | enables |
| GO_ID | GO term | GO:0003677 |
| Evidence_Code | Evidence type | IDA |
| With_From | Supporting IDs | UniProtKB:Q00987 |
| Taxon | Species | taxon:9606 |

---

## Entity Relationships

### Term Hierarchy
- **Cardinality:** N:M
- **Description:** Terms related via is_a, part_of, regulates
- **Key Fields:** goid, relationships

### Annotation
- **Cardinality:** N:M
- **Description:** Genes/proteins annotated to GO terms
- **Key Fields:** goid, annotation.DB_Object_ID

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GO | Gene Ontology | Ontology name |
| GAF | GO Annotation File | File format |
| GPAD | GO Annotation Data | File format |
| GO-CAM | GO Causal Activity Model | Activity model |
| BP | Biological Process | GO aspect |
| MF | Molecular Function | GO aspect |
| CC | Cellular Component | GO aspect |
| OBO | Open Biological Ontology | Format |
| OWL | Web Ontology Language | Format |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
