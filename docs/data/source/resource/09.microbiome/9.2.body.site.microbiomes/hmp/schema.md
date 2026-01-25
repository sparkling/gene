---
id: schema-hmp
title: "HMP (Human Microbiome Project) Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, database, microbiome, body-sites, reference, nih, multi-omics]
---

# HMP Schema Documentation

**Document ID:** SCHEMA-HMP-BODY-SITES
**Version:** 2024.01
**Source Version:** HMP1 + iHMP Phase 2

---

## TL;DR

The Human Microbiome Project (HMP) provides the definitive reference dataset for human-associated microbial communities across 18 body sites. Includes 48+ TB of multi-omics data (16S, WGS, metabolomics, proteomics) from 300+ healthy subjects. Features OSDF (Open Science Data Framework) schemas for standardized data organization, enabling cross-site microbiome comparisons and reference-based analyses.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Body sites | 18 | HMP Protocol |
| Subjects (HMP1) | 300 | HMP Statistics |
| Total samples | 11,000+ | HMP Portal |
| Reference genomes | 3,000+ | HMP Genomes |
| Sequence data | 5.2 Tbp | HMP Statistics |
| Total data volume | 48+ TB | AWS Registry |

---

## Entity Relationship Overview

```
HMP Data Model
  ├── Project
  │     └── Study (HMP1, IBDMDB, T2D, MOMS-PI)
  │           └── Subject
  │                 └── Visit (longitudinal)
  │                       └── Sample
  │                             ├── Body Site
  │                             ├── 16S Data
  │                             ├── WGS Data
  │                             └── Multi-omics
  ├── Body Site Hierarchy
  │     ├── Supersite (5 major regions)
  │     └── Specific site (18 locations)
  └── Reference Data
        ├── Reference genomes
        ├── Taxonomic profiles
        └── Functional profiles
```

---

## Body Site Classification

### Five Major Supersites

| Supersite | Specific Body Sites | Sample Count |
|-----------|---------------------|--------------|
| **Airways** | Anterior nares, Throat | ~2,000 |
| **Gastrointestinal** | Stool | ~3,500 |
| **Oral** | Buccal mucosa, Hard palate, Keratinized gingiva, Palatine tonsils, Saliva, Subgingival plaque, Supragingival plaque, Throat, Tongue dorsum | ~4,000 |
| **Skin** | Left retroauricular crease, Right retroauricular crease, Left antecubital fossa, Right antecubital fossa | ~1,000 |
| **Urogenital** | Mid vagina, Posterior fornix, Vaginal introitus | ~500 |

### Body Site Ontology Mapping

| Body Site | FMA Term | UBERON ID |
|-----------|----------|-----------|
| Stool | Feces | UBERON:0001988 |
| Anterior nares | Nasal cavity | UBERON:0001707 |
| Buccal mucosa | Buccal mucosa | UBERON:0006956 |
| Tongue dorsum | Dorsum of tongue | UBERON:0009471 |
| Supragingival plaque | Supragingival dental plaque | UBERON:0016485 |
| Subgingival plaque | Subgingival dental plaque | UBERON:0016486 |
| Posterior fornix | Posterior vaginal fornix | UBERON:0012247 |

---

## Core Tables/Entities

### Subject Schema

**Description:** HIPAA-compliant subject metadata.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| rand_subject_id | string | Yes | Randomized subject ID (1-32 chars) |
| gender | enum | Yes | male, female, unknown |
| race | enum | No | Ethnicity classification |
| subtype | string | Yes | Subject classification |
| tags | array[string] | Yes | Descriptive tags |

### Sample Schema

**Description:** Physical sample with body site annotation.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| sample_id | string | Yes | Unique sample identifier |
| fma_body_site | string | Yes | FMA ontology term for body site |
| supersite | string | Yes | Major body region |
| body_site | string | Yes | Specific anatomical location |
| int_sample_id | string | No | Center-specific sample ID |
| mixs | object | Yes | MIxS metadata fields |
| collection_date | string | No | ISO format date |
| subtype | string | No | Sample classification |
| tags | array[string] | Yes | Descriptive tags |

### 16S Sequence Set Schema

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| checksums | object | Yes | md5, sha256 hashes |
| format | enum | Yes | fasta, fastq, sff |
| size | integer | Yes | File size in bytes |
| exp_length | integer | Yes | Expected base count |
| seq_model | string | Yes | Sequencing platform |
| study | reference | Yes | Study reference |
| urls | array[string] | Yes | Download URLs |
| variable_region | string | No | V1-V9 region |
| subtype | string | Yes | Sequence classification |

### WGS Sequence Set Schema

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| checksums | object | Yes | File integrity hashes |
| format | enum | Yes | fasta, fastq |
| size | integer | Yes | File size in bytes |
| exp_length | integer | Yes | Expected read count |
| seq_model | string | Yes | Sequencing platform |
| lib_layout | string | Yes | paired, fragment |
| insert_size | integer | No | Library insert size |
| study | reference | Yes | Study reference |
| urls | array[string] | Yes | Download URLs |

---

## Multi-omics Data Types

### Data Type Coverage by Study

| Data Type | HMP1 | IBDMDB | T2D | MOMS-PI |
|-----------|------|--------|-----|---------|
| 16S rRNA | Yes | Yes | Yes | Yes |
| WGS/Metagenomics | Yes | Yes | Yes | Yes |
| Metatranscriptomics | No | Yes | Yes | Yes |
| Host Transcriptomics | No | Yes | Yes | Yes |
| Proteomics | No | Yes | No | Yes |
| Metabolomics | No | Yes | Yes | Yes |
| Lipidomics | No | Yes | No | No |
| Host Cytokines | No | Yes | No | No |

### Abundance Matrix Schema

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| matrix_type | enum | Yes | 16s_community, wgs_community, wgs_functional |
| format | enum | Yes | biom, tbl, csv |
| size | integer | Yes | File size |
| checksums | object | Yes | Integrity hashes |
| parameters | string | No | Generation parameters |
| sop | string | No | SOP documentation URL |

**Matrix Types:**
- `16s_community` - 16S taxonomic profiles
- `wgs_community` - WGS taxonomic profiles
- `wgs_functional` - WGS functional annotations (HUMAnN)
- `microb_metatranscriptome` - Expression profiles
- `host_transcriptome` - Host gene expression
- `microb_metabolome` - Metabolite profiles

---

## Sample Records

### Sample JSON Example

```json
{
  "sample_id": "SRS011084",
  "rand_subject_id": "HMP_SUBJ_001",
  "fma_body_site": "UBERON:0009471",
  "body_site": "Tongue dorsum",
  "supersite": "Oral",
  "collection_date": "2010-03-15",
  "mixs": {
    "biome": "ENVO:00002030",
    "body_product": "UBERON:0000165",
    "env_package": "human-oral",
    "geo_loc_name": "USA: Missouri",
    "lat_lon": "38.6270 -90.1994",
    "material": "oral mucosa"
  },
  "subtype": "tongue_swab",
  "tags": ["HMP1", "oral", "healthy"]
}
```

### Cross-Site Comparison

| Body Site | Avg Shannon | Dominant Phylum | Core Genera |
|-----------|-------------|-----------------|-------------|
| Stool | 3.2-4.0 | Firmicutes | Bacteroides, Faecalibacterium |
| Tongue dorsum | 2.5-3.5 | Firmicutes | Streptococcus, Veillonella |
| Buccal mucosa | 2.0-3.0 | Firmicutes | Streptococcus, Haemophilus |
| Supragingival plaque | 3.0-4.0 | Actinobacteria | Actinomyces, Corynebacterium |
| Anterior nares | 1.5-2.5 | Actinobacteria | Cutibacterium, Corynebacterium |
| Posterior fornix | 0.5-2.0 | Firmicutes | Lactobacillus |

---

## Reference Genome Collection

### Genome Statistics

| Metric | Value |
|--------|-------|
| Total genomes | 3,000+ |
| Complete genomes | 800+ |
| Draft assemblies | 2,200+ |
| Unique species | 1,500+ |
| Average genome size | 3.2 Mb |

### Reference Genome Schema

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| genome_id | string | Yes | HMP genome identifier |
| ncbi_assembly | string | No | NCBI Assembly accession |
| organism_name | string | Yes | Species name |
| strain | string | Yes | Strain designation |
| isolation_source | string | No | Body site of isolation |
| genome_size | integer | Yes | Size in bp |
| gc_content | float | Yes | GC percentage |
| completeness | string | Yes | Complete, Draft |
| gene_count | integer | No | Predicted genes |

---

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| FASTQ | Raw sequences | Sequence reads |
| FASTA | Processed sequences | Assembled data |
| BIOM | Biological Observation Matrix | OTU/ASV tables |
| TSV | Tab-separated values | Metadata, abundance |
| JSON | OSDF schema format | Structured metadata |
| gzip | Compression | File transfer |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/query` | GET | Query OSDF nodes |
| `/api/search/f` | GET | Faceted search |
| `/api/node/{id}` | GET | Get specific node |
| `/api/files/{id}` | GET | File metadata |

### Query Examples

```bash
# Find all oral samples
GET /api/query?supersite=Oral&node_type=sample

# Find 16S data for specific body site
GET /api/query?body_site=Tongue%20dorsum&node_type=16s_raw_seq_set

# Get sample metadata
GET /api/node/SRS011084
```

---

## Cross-References

| Database | ID Type | Usage |
|----------|---------|-------|
| NCBI BioProject | PRJNA | Project identifiers |
| NCBI SRA | SRS/SRR | Sample/run accessions |
| NCBI Taxonomy | TaxID | Organism identification |
| UBERON | Ontology ID | Body site annotation |
| ENVO | Ontology ID | Environment classification |
| FMA | Ontology ID | Anatomical terms |

---

## Glossary

| Term | Definition |
|------|------------|
| Supersite | Major body region (5 categories) |
| Body site | Specific anatomical sampling location |
| MIxS | Minimum Information about any (x) Sequence |
| OSDF | Open Science Data Framework |
| FMA | Foundational Model of Anatomy |
| Alpha diversity | Within-sample diversity |
| Beta diversity | Between-sample diversity |
| Core microbiome | Taxa present across individuals |
| Dysbiosis | Altered microbiome composition |

---

## References

1. Human Microbiome Project Consortium. (2012). Structure, function and diversity of the healthy human microbiome. Nature, 486(7402), 207-214.
2. Human Microbiome Project Consortium. (2012). A framework for human microbiome research. Nature, 486(7402), 215-221.
3. Integrative HMP Research Network Consortium. (2019). The Integrative Human Microbiome Project. Nature, 569(7758), 641-648.
4. OSDF Schemas: https://github.com/ihmpdcc/osdf-schemas
5. HMP DACC Portal: https://hmpdacc.org/
