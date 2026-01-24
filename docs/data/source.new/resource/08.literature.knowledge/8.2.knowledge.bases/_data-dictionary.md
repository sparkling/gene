# 8.2 Knowledge Bases - Data Dictionary

**Subcategory ID:** 8.2
**Subcategory Name:** Knowledge Bases
**Data Sources:** Wikidata, Wikipedia

## Overview

This subcategory integrates structured knowledge from Wikidata (a free knowledge graph) and Wikipedia (the encyclopedia), providing access to entity information, cross-references to biomedical databases, and encyclopedic content.

---

## Unified Fields

### Core Entity Information

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `entity_id` | string | Yes | Unique identifier for the knowledge base entity | `Q123456` (Wikidata), `1234567` (Wikipedia) |
| `name` | string | Yes | Primary name or title of the entity | `TP53` |
| `description` | string | No | Brief description of the entity | `Human gene and protein` |
| `wikidata_qid` | string | No | Wikidata item identifier (Q-number) | `Q21173022` |
| `aliases` | array[string] | No | Alternative names for the entity | `["p53", "Tumor protein p53"]` |

**Source Mappings - Core Entity:**

| Field | Wikidata | Wikipedia |
|-------|----------|-----------|
| `entity_id` | id | pageid |
| `name` | label | title |
| `description` | description | extract |
| `wikidata_qid` | id | wikibase_item |
| `aliases` | aliases | - |

---

### Categories and Classification

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `categories` | array[string] | No | Categories or classes the entity belongs to | `["Tumor suppressor genes", "Genes on human chromosome 17"]` |
| `external_urls` | array[string] | No | External links to other resources | `["https://www.ncbi.nlm.nih.gov/gene/7157"]` |
| `last_modified` | datetime | No | Timestamp of last modification | `2024-01-15T10:30:00Z` |

**Source Mappings - Categories:**

| Field | Wikidata | Wikipedia |
|-------|----------|-----------|
| `categories` | instance_of (P31) | categories |
| `external_urls` | url properties | extlinks |
| `last_modified` | modified | touched |

---

### Biomedical Identifiers

| Field | Data Type | Wikidata Property | Description | Example |
|-------|-----------|-------------------|-------------|---------|
| `biomedical_identifiers.entrez_gene_id` | string | P351 | NCBI Entrez Gene identifier | `7157` |
| `biomedical_identifiers.uniprot_id` | string | P352 | UniProt protein accession | `P04637` |
| `biomedical_identifiers.hgnc_symbol` | string | P353 | HGNC gene symbol | `TP53` |
| `biomedical_identifiers.hgnc_id` | string | P354 | HGNC identifier with prefix | `HGNC:11998` |
| `biomedical_identifiers.ensembl_gene_id` | string | P594 | Ensembl gene identifier | `ENSG00000141510` |
| `biomedical_identifiers.refseq_protein_id` | string | P639 | NCBI RefSeq protein accession | `NP_000537.3` |
| `biomedical_identifiers.refseq_rna_id` | string | P637 | NCBI RefSeq RNA accession | `NM_000546.6` |

**Biomedical Identifiers Object Structure:**

```json
{
  "entrez_gene_id": "7157",
  "uniprot_id": "P04637",
  "hgnc_symbol": "TP53",
  "hgnc_id": "HGNC:11998",
  "ensembl_gene_id": "ENSG00000141510",
  "refseq_protein_id": "NP_000537.3",
  "refseq_rna_id": "NM_000546.6"
}
```

---

### Disease Identifiers

| Field | Data Type | Wikidata Property | Description | Example |
|-------|-----------|-------------------|-------------|---------|
| `disease_identifiers.disease_ontology_id` | string | P699 | Disease Ontology identifier | `DOID:9952` |
| `disease_identifiers.omim_id` | string | P492 | OMIM identifier | `191170` |
| `disease_identifiers.orphanet_id` | string | P1395 | Orphanet rare disease identifier | `1234` |
| `disease_identifiers.icd9` | string | P493 | ICD-9-CM code | `250.00` |
| `disease_identifiers.icd10` | string | P494 | ICD-10 code | `E11` |
| `disease_identifiers.mondo_id` | string | P5270 | MONDO disease ontology identifier | `MONDO:0005148` |
| `disease_identifiers.mesh_id` | string | P486 | MeSH descriptor identifier | `D003920` |
| `disease_identifiers.umls_cui` | string | P2892 | UMLS Concept Unique Identifier | `C0011849` |

---

### Chemical Identifiers

| Field | Data Type | Wikidata Property | Description | Example |
|-------|-----------|-------------------|-------------|---------|
| `chemical_identifiers.cas_number` | string | P231 | Chemical Abstracts Service registry number | `50-78-2` |
| `chemical_identifiers.pubchem_cid` | string | P662 | PubChem Compound identifier | `2244` |
| `chemical_identifiers.chembl_id` | string | P592 | ChEMBL compound identifier | `CHEMBL25` |
| `chemical_identifiers.chebi_id` | string | P652 | ChEBI ontology identifier | `CHEBI:15365` |
| `chemical_identifiers.drugbank_id` | string | P715 | DrugBank drug identifier | `DB00945` |
| `chemical_identifiers.inchi` | string | P234 | IUPAC InChI chemical identifier | `InChI=1S/C9H8O4/...` |
| `chemical_identifiers.inchikey` | string | P235 | InChI hash key | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` |
| `chemical_identifiers.smiles` | string | P233 | SMILES chemical notation | `CC(=O)OC1=CC=CC=C1C(=O)O` |

---

### Biological Relationships

| Field | Data Type | Wikidata Property | Description | Example |
|-------|-----------|-------------------|-------------|---------|
| `biological_relationships.found_in_taxon` | array[string] | P703 | Taxa where entity is found | `["Homo sapiens", "Mus musculus"]` |
| `biological_relationships.genetic_association` | array[string] | P2293 | Gene-disease associations | `["breast cancer", "Li-Fraumeni syndrome"]` |

---

## Wikipedia-Specific Fields

### Metadata

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `wikipedia_metadata.page_id` | integer | Wikipedia page identifier | `1234567` |
| `wikipedia_metadata.revision_id` | integer | Current revision identifier | `987654321` |
| `wikipedia_metadata.namespace` | integer | Wikipedia namespace code (0=article, 14=category) | `0` |
| `wikipedia_metadata.content_model` | string | Content type | `wikitext` |
| `wikipedia_metadata.page_length` | integer | Article length in bytes | `45678` |

### Content Elements

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `infobox_data` | object | Parsed infobox template data | `{Symbol: "TP53", Entrez: "7157"}` |
| `page_links` | array[string] | Links to other Wikipedia articles | `["Cancer", "Tumor suppressor gene"]` |
| `section_titles` | array[string] | Article section headings | `["Function", "Clinical significance", "References"]` |
| `images` | array[object] | Images included in the article | `[{title, url}]` |
| `dbpedia_uri` | string | DBpedia resource URI | `http://dbpedia.org/resource/TP53` |

**Image Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `title` | string | Image file title |
| `url` | string | Image URL |

---

## Wikidata Property Reference

### Gene/Protein Properties

| Property | Name | Description | Format |
|----------|------|-------------|--------|
| P351 | Entrez Gene ID | NCBI Gene identifier | `[0-9]+` |
| P352 | UniProt protein ID | UniProtKB accession | `[A-Z][0-9][A-Z0-9]{3}[0-9]` |
| P353 | HGNC gene symbol | Official gene symbol | `[A-Z0-9]+` |
| P354 | HGNC ID | HGNC identifier | `HGNC:[0-9]+` |
| P594 | Ensembl gene ID | Ensembl gene identifier | `ENSG[0-9]{11}` |
| P637 | RefSeq RNA ID | RefSeq RNA accession | `[NX]M_[0-9]+` |
| P639 | RefSeq protein ID | RefSeq protein accession | `[NX]P_[0-9]+` |

### Disease Properties

| Property | Name | Description | Format |
|----------|------|-------------|--------|
| P486 | MeSH descriptor ID | Medical Subject Headings | `[A-Z][0-9]+` |
| P492 | OMIM ID | Online Mendelian Inheritance | `[0-9]{6}` |
| P493 | ICD-9-CM | International Classification of Diseases 9 | varies |
| P494 | ICD-10 | International Classification of Diseases 10 | varies |
| P699 | Disease Ontology ID | Disease Ontology | `DOID:[0-9]+` |
| P1395 | Orphanet ID | Rare disease identifier | `[0-9]+` |
| P2892 | UMLS CUI | Unified Medical Language System | `C[0-9]+` |
| P5270 | MONDO ID | Mondo Disease Ontology | `MONDO:[0-9]+` |

### Chemical Properties

| Property | Name | Description | Format |
|----------|------|-------------|--------|
| P231 | CAS Registry Number | Chemical Abstracts Service | `[0-9]+-[0-9]+-[0-9]+` |
| P233 | canonical SMILES | Simplified molecular-input line-entry | varies |
| P234 | InChI | IUPAC International Chemical Identifier | `InChI=...` |
| P235 | InChIKey | InChI hash key | `[A-Z]{14}-[A-Z]{10}-[A-Z]` |
| P592 | ChEMBL ID | ChEMBL compound identifier | `CHEMBL[0-9]+` |
| P652 | ChEBI ID | Chemical Entities of Biological Interest | `[0-9]+` |
| P662 | PubChem CID | PubChem Compound ID | `[0-9]+` |
| P715 | DrugBank ID | DrugBank drug identifier | `DB[0-9]+` |

### Biological Relationship Properties

| Property | Name | Description |
|----------|------|-------------|
| P31 | instance of | Classification of entity type |
| P279 | subclass of | Hierarchical classification |
| P703 | found in taxon | Organism where entity exists |
| P2293 | genetic association | Gene-disease relationship |

---

## Source Mappings Summary

### Wikidata Field Mappings

| Unified Field | Wikidata Source |
|---------------|-----------------|
| `entity_id` | id |
| `name` | label |
| `description` | description |
| `aliases` | aliases |
| `last_modified` | modified |
| `categories` | P31 (instance_of) |
| `external_urls` | url properties |
| `biomedical_identifiers.*` | P351, P352, P353, P354, P594, P637, P639 |
| `disease_identifiers.*` | P486, P492, P493, P494, P699, P1395, P2892, P5270 |
| `chemical_identifiers.*` | P231, P233, P234, P235, P592, P652, P662, P715 |
| `biological_relationships.*` | P703, P2293 |

### Wikipedia Field Mappings

| Unified Field | Wikipedia Source |
|---------------|------------------|
| `entity_id` | pageid |
| `name` | title |
| `description` | extract |
| `wikidata_qid` | wikibase_item |
| `last_modified` | touched |
| `categories` | categories |
| `external_urls` | extlinks |
| `wikipedia_metadata.revision_id` | revid |
| `wikipedia_metadata.namespace` | ns |
| `wikipedia_metadata.content_model` | contentmodel |
| `wikipedia_metadata.page_length` | length |
| `infobox_data` | infobox |
| `page_links` | links |
| `section_titles` | sections |
| `images` | images |
| `dbpedia_uri` | dbpedia_uri |

---

## Data Source Metadata

| Field | Data Type | Description |
|-------|-----------|-------------|
| `_source.primary_source` | string | Name of the primary data source (Wikidata or Wikipedia) |
| `_source.source_id` | string | Original ID in the source system |
| `_source.extraction_date` | date | Date when data was extracted |
| `_source.source_version` | string | Version of the source data/API |

---

## Query Endpoints

### Wikidata

- **SPARQL Endpoint:** `https://query.wikidata.org/sparql`
- **REST API:** `https://www.wikidata.org/wiki/Special:EntityData/{QID}.json`

### Wikipedia

- **MediaWiki API:** `https://en.wikipedia.org/w/api.php`
- **REST API:** `https://en.wikipedia.org/api/rest_v1/page/summary/{title}`
