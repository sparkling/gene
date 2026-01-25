# 4.5 Gene Function & Ontology - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| Subcategory ID | 4.5 |
| Subcategory Name | Gene Function & Ontology |
| Data Sources | Gene Ontology, MSigDB |
| Schema ID | `https://gene.ai/schemas/4.5-gene-function-ontology.json` |

## Unified Fields

These fields are harmonized across all data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `term_id` | string | Required (1:1) | Primary identifier for the ontology term or gene set | Gene Ontology, MSigDB | GO: `GO:0008152`, MSigDB: `M5930` |
| `term_name` | string | Required (1:1) | Human-readable term name | Gene Ontology, MSigDB | `metabolic process`, `HALLMARK_APOPTOSIS` |
| `definition` | string | Optional (1:1) | Textual definition of the term | Gene Ontology, MSigDB | `The chemical reactions and pathways...` |
| `namespace` | string | Optional (1:1) | Ontology namespace or collection | Gene Ontology, MSigDB | `biological_process`, `H` |
| `genes` | array&lt;string&gt; | Optional (1:N) | Genes associated with the term | Gene Ontology, MSigDB | `["CASP3", "CASP8", "BAX", "BCL2"]` |

## Gene Ontology-Specific Fields

### Core Ontology Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `go_id` | string | Optional (1:1) | Gene Ontology identifier (pattern: `GO:[0-9]{7}`) | `GO:0008152` |
| `go_namespace` | enum | Optional (1:1) | GO namespace/aspect | `biological_process`, `cellular_component`, `molecular_function` |
| `is_a` | array&lt;string&gt; | Optional (1:N) | Parent terms in hierarchy (subclass relationship) | `["GO:0008150"]` |
| `part_of` | array&lt;string&gt; | Optional (1:N) | Part-of relationships | `["GO:0009987"]` |
| `regulates` | array&lt;string&gt; | Optional (1:N) | Regulatory relationships | `["GO:0006915"]` |
| `synonyms` | array&lt;string&gt; | Optional (1:N) | Alternative names (EXACT, NARROW, BROAD) | `["metabolism"]` |
| `xrefs` | array&lt;string&gt; | Optional (1:N) | Cross-references to external databases | `["Wikipedia:Metabolism"]` |
| `is_obsolete` | boolean | Optional (1:1) | Whether the term is obsolete | `false` |
| `replaced_by` | string | Optional (1:1) | Replacement term for obsolete terms | `GO:0000001` |

### GO Annotation Format (GAF 2.2) Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `DB` | string | Required (1:1) | Database contributing the annotation | `UniProtKB` |
| `DB_Object_ID` | string | Required (1:1) | Unique identifier in source database | `P04637` |
| `DB_Object_Symbol` | string | Required (1:1) | Gene/protein symbol | `TP53` |
| `Qualifier` | string | Optional (1:N) | Annotation qualifier | `NOT`, `contributes_to`, `colocalizes_with`, `enables`, `involved_in`, `located_in` |
| `GO_ID` | string | Required (1:1) | GO term identifier | `GO:0003700` |
| `DB_Reference` | string | Required (1:N) | Reference for annotation (PMID, GO_REF) | `PMID:8657117` |
| `Evidence_Code` | string | Required (1:1) | Evidence code from ECO | `IDA` |
| `With_From` | string | Optional (1:N) | Additional supporting identifiers | `UniProtKB:P04637` |
| `DB_Object_Type` | enum | Required (1:1) | Type of annotated entity | `gene`, `protein`, `transcript`, `complex` |
| `Taxon` | string | Required (1:2) | NCBI Taxonomy ID | `taxon:9606` |
| `Date` | string | Required (1:1) | Annotation date (YYYYMMDD) | `20230415` |
| `Assigned_By` | string | Required (1:1) | Group that made annotation | `UniProt` |

### GO-CAM Model Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `model_id` | string | Required (1:1) | GO-CAM model identifier | `gomodel:5fce9b7300002411` |
| `title` | string | Required (1:1) | Model title | `Insulin signaling pathway` |
| `state` | enum | Required (1:1) | Curation state | `development`, `review`, `production` |
| `activities` | array&lt;activity&gt; | Optional (1:N) | Molecular activities in the model | `[{"gene_product": "UniProtKB:P06213", "molecular_function": "GO:0004716"}]` |
| `causal_relations` | array&lt;relation&gt; | Optional (1:N) | Causal relationships between activities | `[{"predicate": "RO:0002413", "predicate_label": "directly_provides_input_for"}]` |
| `enabled_by` | string | Optional (1:1) | Gene product enabling the activity | `UniProtKB:P06213` |

## MSigDB-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `msigdb_id` | string | Optional (1:1) | MSigDB gene set identifier | `M5930` |
| `standard_name` | string | Optional (1:1) | Standardized gene set name | `HALLMARK_APOPTOSIS` |
| `systematic_name` | string | Optional (1:1) | Systematic identifier | `M5930` |
| `msigdb_collection` | enum | Optional (1:1) | MSigDB collection category | `H`, `C1`, `C2`, `C3`, `C4`, `C5`, `C6`, `C7`, `C8` |
| `msigdb_subcollection` | string | Optional (1:1) | MSigDB subcollection | `CP:KEGG`, `CP:REACTOME`, `GO:BP`, `GO:CC`, `GO:MF` |
| `brief_description` | string | Optional (1:1) | Short description of the gene set | `Genes mediating programmed cell death...` |
| `full_description` | string | Optional (1:1) | Detailed description with references | Full text with citations |
| `external_links` | array&lt;string&gt; | Optional (1:N) | URLs to external resources | `["http://..."]` |
| `pmid` | integer | Optional (1:1) | Associated PubMed ID | `26771021` |
| `gene_count` | integer | Optional (1:1) | Number of genes in the set | `161` |

## Evidence Codes

### Experimental Evidence

| Code | Name | Weight |
|------|------|--------|
| EXP | Inferred from Experiment | High |
| IDA | Inferred from Direct Assay | High |
| IPI | Inferred from Physical Interaction | High |
| IMP | Inferred from Mutant Phenotype | High |
| IGI | Inferred from Genetic Interaction | High |
| IEP | Inferred from Expression Pattern | Medium |
| HTP | High Throughput Experiment | Medium |

### Computational Evidence

| Code | Name | Weight |
|------|------|--------|
| ISS | Inferred from Sequence Similarity | Medium |
| ISO | Inferred from Sequence Orthology | Medium |
| IBA | Inferred from Biological Ancestor | Medium |
| RCA | Reviewed Computational Analysis | Medium |

### Author/Curator Evidence

| Code | Name | Weight |
|------|------|--------|
| TAS | Traceable Author Statement | Medium |
| NAS | Non-traceable Author Statement | Low |
| IC | Inferred by Curator | Medium |

### Electronic Evidence

| Code | Name | Weight |
|------|------|--------|
| IEA | Inferred from Electronic Annotation | Low |

## MSigDB Collections

| ID | Name | Description | Approximate Count |
|----|------|-------------|-------------------|
| H | Hallmark | Well-defined biological states | 50 |
| C1 | Positional | Chromosome cytogenetic bands | 299 |
| C2 | Curated | Curated pathway databases | 6,300+ |
| C3 | Motif | TF and miRNA targets | 3,700+ |
| C4 | Computational | Cancer-related signatures | 850+ |
| C5 | GO | Gene Ontology terms | 15,000+ |
| C6 | Oncogenic | Oncogenic pathway signatures | 189 |
| C7 | Immunologic | Immunologic signatures | 5,200+ |
| C8 | Cell Type | Cell type markers | 700+ |

## Source Metadata

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `_source.database` | string | Required | Source database name | `Gene Ontology`, `MSigDB` |
| `_source.version` | string | Optional | Database version | `2024-01-01`, `2023.2` |
| `_source.access_date` | date | Optional | Date data was retrieved | `2026-01-24` |
| `_source.original_id` | string | Optional | Original identifier in source | `GO:0008152` |

## Field Mappings by Source

### Gene Ontology Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `id` | `go_id` |
| `name` | `term_name` |
| `namespace` | `go_namespace` |
| `def` | `definition` |
| `is_a` | `is_a` |
| `part_of` | `part_of` |
| `regulates` | `regulates` |
| `synonym` | `synonyms` |
| `xref` | `xrefs` |
| `is_obsolete` | `is_obsolete` |
| `replaced_by` | `replaced_by` |

### MSigDB Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `GENESET_ID` | `msigdb_id` |
| `STANDARD_NAME` | `standard_name` |
| `SYSTEMATIC_NAME` | `systematic_name` |
| `COLLECTION` | `msigdb_collection` |
| `SUBCOLLECTION` | `msigdb_subcollection` |
| `BRIEF_DESCRIPTION` | `brief_description` |
| `FULL_DESCRIPTION` | `full_description` |
| `EXTERNAL_LINKS` | `external_links` |
| `PMID` | `pmid` |
| `MEMBERS` | `genes` |
| `MEMBERS_SYMBOLIZED` | `genes` |
| `SIZE` | `gene_count` |
