# Data Dictionary: 7.1 Protein Sequences & Annotations

## Overview

This data dictionary documents all fields for protein sequence and annotation data, integrating UniProt and RefSeq sources into a unified schema.

**Subcategory ID:** 7.1
**Data Sources:** UniProt, RefSeq
**Schema ID:** https://gene.taxonomy/schemas/7.1-protein-sequences-annotations

---

## Unified Fields

These fields are harmonized across all data sources with consistent semantics.

| Field | Data Type | Cardinality | Description | Sources | Example |
|-------|-----------|-------------|-------------|---------|---------|
| `accession` | string | Required (1:1) | Primary stable identifier for the protein entry. UniProt uses format like P04637, RefSeq uses prefix-based format like NP_000537.3 | UniProt, RefSeq | `P04637`, `NP_000537.3` |
| `sequence` | string | Required (1:1) | Complete amino acid sequence using standard single-letter codes (ACDEFGHIKLMNPQRSTVWY) | UniProt, RefSeq | `MEEPQSDPSVEPPLSQETFSDLWKLL...` |
| `length` | integer | Required (1:1) | Integer count of residues in the canonical protein sequence | UniProt, RefSeq | `393`, `1234` |
| `organism` | string | Required (1:1) | Scientific name of the organism from which the protein originates | UniProt, RefSeq | `Homo sapiens`, `Mus musculus` |
| `taxonomy` | array[string] | Optional (1:N) | Ordered array of taxonomic ranks representing the evolutionary classification | UniProt, RefSeq | `["Eukaryota", "Metazoa", "Chordata", "Mammalia", "Primates", "Hominidae", "Homo"]` |
| `gene` | string | Optional (1:1) | Official gene symbol (e.g., HGNC approved symbol for human genes) | UniProt, RefSeq | `TP53`, `BRCA1`, `EGFR` |
| `gene_id` | integer | Optional (1:1) | Numeric identifier from NCBI Gene database linking to gene record | UniProt, RefSeq | `7157`, `672`, `1956` |
| `features` | array[object] | Optional (1:N) | Collection of annotations mapped to specific sequence positions including domains, binding sites, PTMs, and variants | UniProt, RefSeq | `[{"type": "Domain", "location": {"start": 102, "end": 292}, "description": "p53 DNA-binding"}]` |
| `dbxrefs` | array[object] | Optional (1:N) | Links to related entries in other databases including PDB, Ensembl, HGNC, OMIM, etc. | UniProt, RefSeq | `[{"database": "PDB", "id": "1TUP"}, {"database": "HGNC", "id": "11998"}]` |
| `definition` | string | Optional (1:1) | Human-readable description of the protein including recommended name and organism | UniProt, RefSeq | `Cellular tumor antigen p53 [Homo sapiens]` |

---

## Features Object Structure

The `features` array contains objects with the following structure:

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `type` | string | Feature type (Domain, Site, Variant, etc.) | `Domain` |
| `location.start` | integer | Start position (1-based) | `102` |
| `location.end` | integer | End position (1-based) | `292` |
| `description` | string | Feature description | `p53 DNA-binding` |
| `evidence` | string | Evidence supporting the feature | `ECO:0000269` |

---

## Database Cross-References Object Structure

The `dbxrefs` array contains objects with the following structure:

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `database` | string | Database name | `PDB`, `HGNC`, `Ensembl` |
| `id` | string | Identifier in the external database | `1TUP`, `11998` |

---

## UniProt-Specific Fields

These fields are only available when the source is UniProt.

| Field | Data Type | Cardinality | Description | Example |
|-------|-----------|-------------|-------------|---------|
| `entry_name` | string | Optional (1:1) | Human-readable mnemonic combining protein and organism (id field) | `P53_HUMAN`, `BRCA1_HUMAN` |
| `reviewed` | boolean | Optional (1:1) | Indicates manual curation status - reviewed entries are in Swiss-Prot, unreviewed in TrEMBL | `true`, `false` |
| `protein_existence` | string (enum) | Optional (1:1) | PE level indicating confidence in protein existence based on evidence type (1-5) | `1: Evidence at protein level`, `3: Inferred from homology` |
| `molecular_weight` | integer | Optional (1:1) | Calculated molecular weight of the canonical sequence in Daltons | `43653`, `126000` |
| `checksum_crc64` | string | Optional (1:1) | 64-bit cyclic redundancy check for sequence verification | `AD5C149FD8106131` |
| `checksum_md5` | string | Optional (1:1) | MD5 hash for sequence integrity verification | `a1b2c3d4e5f6...` |
| `protein_description` | object | Optional (1:1) | Structured protein nomenclature including recommended, alternative, and submitted names | `{"recommended_name": "Cellular tumor antigen p53"}` |
| `genes` | array[object] | Optional (1:N) | Complete gene nomenclature including official symbol, synonyms, ordered locus names, and ORF names | `[{"gene_name": "TP53", "synonyms": ["P53"]}]` |
| `comments` | array[object] | Optional (1:N) | Curated annotations covering function, disease, subcellular location, PTM, etc. | `[{"comment_type": "FUNCTION", "text": "Acts as a tumor suppressor..."}]` |
| `keywords` | array[string] | Optional (1:N) | Standardized terms from UniProt keyword vocabulary for categorization | `["Tumor suppressor", "DNA-binding", "Nucleus"]` |
| `references` | array[object] | Optional (1:N) | Publications supporting the entry annotations | `[{"title": "...", "journal": "Nature"}]` |

### Protein Description Object Structure

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `recommended_name` | string | Official recommended protein name | `Cellular tumor antigen p53` |
| `alternative_names` | array[string] | Alternative names for the protein | `["Antigen NY-CO-13", "Phosphoprotein p53"]` |
| `submitted_names` | array[string] | Names submitted by sequence depositors | `["p53 tumor suppressor"]` |

### Genes Object Structure

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `gene_name` | string | Primary gene symbol | `TP53` |
| `synonyms` | array[string] | Gene name synonyms | `["P53", "LFS1"]` |
| `ordered_locus_names` | array[string] | Systematic locus identifiers | `["At1g01010"]` |
| `orf_names` | array[string] | Open reading frame identifiers | `["YDL131W"]` |

### Comments Object Structure

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `comment_type` | string | Type of annotation (FUNCTION, DISEASE, SUBCELLULAR LOCATION, etc.) | `FUNCTION` |
| `text` | string | Annotation text content | `Acts as a tumor suppressor in many tumor types...` |

### References Object Structure

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `title` | string | Publication title | `The p53 tumor suppressor gene...` |
| `journal` | string | Journal name | `Nature` |
| `pubmed_id` | string | PubMed identifier | `12345678` |

---

## RefSeq-Specific Fields

These fields are only available when the source is RefSeq.

| Field | Data Type | Cardinality | Description | Example |
|-------|-----------|-------------|-------------|---------|
| `gi` | integer | Optional (1:1) | Deprecated numeric identifier from NCBI GenInfo system (legacy) | `4557757` |
| `locus` | string | Optional (1:1) | Locus identifier from the source genome annotation | `TP53`, `NP_000537` |
| `coded_by` | string | Optional (1:1) | Reference to the nucleotide sequence encoding this protein with CDS coordinates | `NM_000546.6:203..1384` |
| `chromosome` | string | Optional (1:1) | Chromosome on which the encoding gene is located | `17`, `X`, `MT` |
| `map_location` | string | Optional (1:1) | Cytogenetic band location on the chromosome | `17p13.1`, `Xq28` |
| `annotation_status` | string (enum) | Optional (1:1) | Level of manual review and validation applied to the record | `VALIDATED`, `REVIEWED`, `PREDICTED` |
| `version` | integer | Optional (1:1) | Integer version incremented when sequence changes | `1`, `3`, `5` |

### Annotation Status Values

| Value | Description |
|-------|-------------|
| `VALIDATED` | Highest level of validation with experimental support |
| `REVIEWED` | Manually reviewed by NCBI staff |
| `PREDICTED` | Computationally predicted, not manually reviewed |

---

## Source Metadata

The `_source` object provides metadata about data provenance.

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `primary_source` | string | Name of the primary data source | `UniProt`, `RefSeq` |
| `source_id` | string | Original identifier in the source | `P04637` |
| `extraction_date` | string (date) | Date data was extracted | `2026-01-24` |
| `source_version` | string | Version of the source database | `2026_01` |

---

## Field Mappings

### UniProt to Unified Schema

| UniProt Field | Unified Field |
|---------------|---------------|
| `accession` | `accession` |
| `sequence.value` | `sequence` |
| `sequence.length` | `length` |
| `organism.scientificName` | `organism` |
| `organism.lineage` | `taxonomy` |
| `genes[0].geneName.value` | `gene` |
| `dbReferences[database=GeneID].id` | `gene_id` |
| `features` | `features` |
| `dbReferences` | `dbxrefs` |
| `proteinDescription.recommendedName.fullName` | `definition` |
| `id` | `entry_name` |
| `reviewed` | `reviewed` |
| `proteinExistence` | `protein_existence` |
| `sequence.molWeight` | `molecular_weight` |
| `sequence.crc64` | `checksum_crc64` |
| `sequence.md5` | `checksum_md5` |
| `proteinDescription` | `protein_description` |
| `genes` | `genes` |
| `comments` | `comments` |
| `keywords` | `keywords` |
| `references` | `references` |

### RefSeq to Unified Schema

| RefSeq Field | Unified Field |
|--------------|---------------|
| `accession` | `accession` |
| `sequence` | `sequence` |
| `length` | `length` |
| `organism` | `organism` |
| `taxonomy` | `taxonomy` |
| `gene` | `gene` |
| `db_xref[GeneID]` | `gene_id` |
| `features` | `features` |
| `db_xref` | `dbxrefs` |
| `definition` | `definition` |
| `gi` | `gi` |
| `locus` | `locus` |
| `coded_by` | `coded_by` |
| `chromosome` | `chromosome` |
| `map` | `map_location` |
| `status` | `annotation_status` |
| `version` | `version` |

---

## Required Fields

The following fields are required for a valid record:

- `accession`
- `sequence`
- `length`
- `organism`
