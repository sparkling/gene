# 9.1 Gut Microbiome - Data Dictionary

## Overview

This data dictionary documents the unified schema for gut microbiome data, integrating GMrepo, gutMGene, HMP, and MetaHIT data sources.

**Subcategory ID:** 9.1
**Subcategory Name:** Gut Microbiome
**Data Sources:** GMrepo, gutMGene, HMP, MetaHIT

## Required Fields

| Field | Type | Description | Sources |
|-------|------|-------------|---------|
| sample_id | string | BioSample accession for individual biological samples | GMrepo, HMP, MetaHIT |
| taxon_id | integer | NCBI Taxonomy identifier for microbial species or genus | GMrepo, gutMGene, HMP, MetaHIT |
| taxon_name | string | Scientific name of the microbial taxon | GMrepo, gutMGene, HMP, MetaHIT |

## Unified Fields

### Sample Identifiers

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| project_id | string | Optional | BioProject accession identifying a collection of related samples | `PRJEB6070` | GMrepo, HMP, MetaHIT |
| sample_id | string | Required | BioSample accession for individual biological samples | `SAMN1234567` | GMrepo, HMP, MetaHIT |
| run_id | string | Optional | Sequencing run accession identifier from SRA/ENA archives | `SRR1234567` | GMrepo, HMP |
| subject_id | string | Optional | Randomized identifier for study participant (HIPAA-compliant) | `HMP2_J12345` | HMP |

### Taxonomic Information

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| taxon_id | integer | Required | NCBI Taxonomy identifier for microbial species or genus | `816` | GMrepo, gutMGene, HMP, MetaHIT |
| taxon_name | string | Required | Scientific name of the microbial taxon | `Bacteroides` | GMrepo, gutMGene, HMP, MetaHIT |
| rank | enum | Optional | Taxonomic rank (kingdom to species). Values: `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` | `genus` | GMrepo, MetaHIT |
| enterotype | enum | Optional | Gut microbiome community type classification. Values: `Bacteroides`, `Prevotella`, `Ruminococcus` | `Bacteroides` | MetaHIT |

### Abundance Data

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| relative_abundance | float | Optional | Proportion of a taxon in a sample (0-1 scale) | `0.2534` | GMrepo, HMP, MetaHIT |
| read_count | integer | Optional | Number of sequencing reads mapped to a taxon | `15234` | GMrepo, HMP, MetaHIT |

### Host Metadata

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| host_metadata.age | integer | Optional | Age of the human subject in years | `45` | GMrepo, HMP, MetaHIT |
| host_metadata.sex | enum | Optional | Biological sex of the human subject. Values: `male`, `female`, `unknown` | `male` | GMrepo, HMP, MetaHIT |
| host_metadata.bmi | float | Optional | Body mass index | `28.5` | GMrepo, MetaHIT |
| host_metadata.country | string | Optional | Country of sample origin | `USA` | GMrepo, HMP, MetaHIT |
| host_metadata.race | string | Optional | Ethnicity classification | `caucasian` | HMP |

### Phenotype Information

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| phenotype | string | Optional | Disease or health status annotation for a sample | `Type 2 diabetes` | GMrepo, HMP |
| phenotype_mesh_id | string | Optional | Medical Subject Headings identifier for standardized phenotype | `D003924` | GMrepo |
| ibd_status | enum | Optional | Inflammatory bowel disease status. Values: `CD`, `UC`, `healthy` | `healthy` | MetaHIT |

### Phenotype Associations (Array)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| phenotype_associations[].phenotype | string | Optional | Associated phenotype name | `Type 2 diabetes` | GMrepo |
| phenotype_associations[].direction | enum | Optional | Direction of change. Values: `increased`, `decreased`, `altered` | `increased` | GMrepo |
| phenotype_associations[].effect_size | float | Optional | Magnitude of association | `0.45` | GMrepo |
| phenotype_associations[].p_value | float | Optional | Statistical significance | `0.001` | GMrepo |
| phenotype_associations[].study_count | integer | Optional | Number of supporting studies | `5` | GMrepo |
| phenotype_associations[].marker_count | integer | Optional | Number of marker taxa for phenotype | `12` | GMrepo |

### Gene Associations (Array)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| gene_associations[].gene_id | string | Optional | Gene catalog identifier for microbial gene | `MH0001_GL0000001` | gutMGene |
| gene_associations[].gene_symbol | string | Optional | Host gene symbol (official symbol) | `IL10` | gutMGene |
| gene_associations[].direction | enum | Optional | Direction of change. Values: `up`, `down`, `increased`, `decreased`, `altered` | `up` | gutMGene |
| gene_associations[].fold_change | float | Optional | Expression fold change between conditions | `2.5` | gutMGene |
| gene_associations[].p_value | float | Optional | Statistical significance | `0.001` | gutMGene |
| gene_associations[].tissue | string | Optional | Tissue or organ where gene expression was measured | `colon` | gutMGene |
| gene_associations[].experiment_type | string | Optional | Type of experimental approach used | `mono-colonization` | gutMGene |
| gene_associations[].model | string | Optional | Model system used in the experiment | `germ-free mice` | gutMGene |
| gene_associations[].pmid | integer | Optional | PubMed identifier | `28123456` | gutMGene |

### Sequencing Information

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| sequencing_type | enum | Optional | Method used for microbiome profiling. Values: `16S`, `WGS` | `16S` | GMrepo, HMP, MetaHIT |
| matrix_type | enum | Optional | Classification of processed abundance data type. Values: `16s_community`, `wgs_community`, `wgs_functional` | `16s_community` | HMP |
| visit_number | integer | Optional | Sequential number indicating order of clinical visits | `1` | HMP |
| fma_body_site | string | Optional | Foundational Model of Anatomy term for sample location | `UBERON:0001988` | HMP |

### Sequencing Metadata (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| sequencing_metadata.seq_model | string | Optional | Sequencing instrument model | `Illumina MiSeq` | HMP |
| sequencing_metadata.exp_length | integer | Optional | Expected number of bases or reads in sequence file | `10000000` | HMP |
| sequencing_metadata.coverage | float | Optional | Read coverage for gene | `25.5` | MetaHIT |
| sequencing_metadata.sequencing_depth | integer | Optional | Total reads in sample | `1000000` | MetaHIT |

### MIxS Metadata (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| mixs_metadata.biome | string | Optional | ENVO biome term | `ENVO:00000446` | HMP |
| mixs_metadata.collection_date | string | Optional | Date of sample collection | `2020-01-15` | HMP |
| mixs_metadata.geo_loc_name | string | Optional | Geographic location name | `USA: Missouri` | HMP |
| mixs_metadata.env_package | string | Optional | Environmental package type | `human-gut` | HMP |

### Functional Annotations (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| functional_annotations.kegg_ko | array[string] | Optional | KEGG Ortholog identifiers | `["K00001", "K00002"]` | MetaHIT |
| functional_annotations.eggnog | array[string] | Optional | eggNOG ortholog group annotations | `["COG0001"]` | MetaHIT |
| functional_annotations.cog | array[string] | Optional | COG category assignments | `["COG0001"]` | MetaHIT |

### Microbe Characteristics (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| microbe_characteristics.gram_stain | enum | Optional | Gram staining classification. Values: `positive`, `negative`, `variable` | `negative` | gutMGene |
| microbe_characteristics.oxygen_req | enum | Optional | Oxygen requirements of microbe. Values: `aerobic`, `anaerobic`, `facultative` | `anaerobic` | gutMGene |
| microbe_characteristics.host_species | string | Optional | Host organism species | `Homo sapiens` | gutMGene |

### File Integrity

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| checksums.md5 | string | Optional | MD5 hash for file integrity | `d41d8cd98f00b204e9800998ecf8427e` | HMP |
| checksums.sha1 | string | Optional | SHA1 hash for file integrity | - | HMP |
| checksums.sha256 | string | Optional | SHA256 hash for file integrity | - | HMP |

### Other Fields

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| antibiotics_usage | enum | Optional | Recent antibiotic use status. Values: `yes`, `no`, `unknown` | `no` | GMrepo |
| pmid | array[integer] | Optional | PubMed identifiers for literature evidence | `[28123456]` | gutMGene, HMP |

### Source Metadata

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| _source.primary_source | string | Optional | Name of the data source | `GMrepo` | All |
| _source.source_id | string | Optional | Original identifier in source | `GMrepo_001` | All |
| _source.extraction_date | string (date) | Optional | Date of data extraction | `2026-01-24` | All |
| _source.source_version | string | Optional | Version of the source database | `1.0` | All |

## Source-Specific Fields

### GMrepo-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| antibiotics_usage | string | Recent antibiotic use status | `yes`, `no`, `unknown` |
| marker_count | integer | Number of marker taxa for a phenotype | `12` |
| study_count | integer | Number of supporting studies for marker association | `5` |

### gutMGene-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| association_id | integer | Unique microbe-gene association identifier | `1234` |
| host_species | string | Host organism species | `Homo sapiens`, `Mus musculus` |
| gram_stain | string | Gram staining classification of microbe | `positive`, `negative` |
| oxygen_req | string | Oxygen requirements of microbe | `aerobic`, `anaerobic`, `facultative` |
| condition | string | Experimental condition description | - |

### HMP-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| race | enum | Ethnicity classification | `caucasian`, `african_american`, `asian` |
| interval | integer | Days since last visit | `30` |
| clinic_id | string | Clinical site identifier | `CLINIC01` |
| seq_model | string | Sequencing instrument model | `Illumina MiSeq`, `Illumina HiSeq` |
| exp_length | integer | Expected number of bases or reads in sequence file | `10000000` |
| subtype | string | Classification distinguishing host vs microbiome data | - |
| pride_id | string | PRIDE archive identifier for proteomics data | - |
| analyzer | string | Mass analyzer type for proteomics | - |

### MetaHIT-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| ibd_status | string | Inflammatory bowel disease status | `CD`, `UC`, `healthy` |
| cog | string | COG category assignment | `COG0001` |
| sample_origin | string | Source sample ID for gene catalog entry | - |
| coverage | float | Read coverage for gene | `25.5` |
| sequencing_depth | integer | Total reads in sample | `1000000` |

## Cross-References

| Database | ID Field | Usage |
|----------|----------|-------|
| NCBI Taxonomy | taxon_id | Microbial species identification |
| MeSH | mesh_id | Disease/phenotype standardization |
| NCBI SRA | run_id | Raw sequence data access |
| NCBI BioSample | sample_id | Sample metadata |
| NCBI Gene | gene_id | Host gene identification |
| PubMed | pmid | Literature evidence |
| KEGG | kegg_ko | Functional annotation |

## Field Mappings by Source

### GMrepo Mappings

| Source Field | Unified Field |
|--------------|---------------|
| project_id | project_id |
| run_id | run_id |
| taxon_id | taxon_id |
| taxon_name | taxon_name |
| relative_abundance | relative_abundance |
| read_count | read_count |
| phenotype | phenotype |
| mesh_id | phenotype_mesh_id |
| sequencing_type | sequencing_type |
| host_age | host_metadata.age |
| host_sex | host_metadata.sex |
| host_bmi | host_metadata.bmi |
| country | host_metadata.country |
| direction | phenotype_associations[].direction |
| effect_size | phenotype_associations[].effect_size |
| p_value | phenotype_associations[].p_value |
| marker_count | phenotype_associations[].marker_count |
| study_count | phenotype_associations[].study_count |
| rank | rank |
| antibiotics_usage | antibiotics_usage |

### gutMGene Mappings

| Source Field | Unified Field |
|--------------|---------------|
| taxon_id | taxon_id |
| taxon_name | taxon_name |
| gene_id | gene_associations[].gene_id |
| gene_symbol | gene_associations[].gene_symbol |
| direction | gene_associations[].direction |
| fold_change | gene_associations[].fold_change |
| p_value | gene_associations[].p_value |
| tissue | gene_associations[].tissue |
| experiment_type | gene_associations[].experiment_type |
| model | gene_associations[].model |
| pmid | pmid |
| gram_stain | microbe_characteristics.gram_stain |
| oxygen_req | microbe_characteristics.oxygen_req |
| host_species | microbe_characteristics.host_species |

### HMP Mappings

| Source Field | Unified Field |
|--------------|---------------|
| project_id | project_id |
| sample_id | sample_id |
| run_id | run_id |
| taxon_id | taxon_id |
| taxon_name | taxon_name |
| relative_abundance | relative_abundance |
| read_count | read_count |
| phenotype | phenotype |
| sequencing_type | sequencing_type |
| host_age | host_metadata.age |
| host_sex | host_metadata.sex |
| country | host_metadata.country |
| race | host_metadata.race |
| rand_subject_id | subject_id |
| visit_number | visit_number |
| fma_body_site | fma_body_site |
| mixs | mixs_metadata |
| checksums | checksums |
| matrix_type | matrix_type |
| seq_model | sequencing_metadata.seq_model |
| exp_length | sequencing_metadata.exp_length |
| pmid | pmid |

### MetaHIT Mappings

| Source Field | Unified Field |
|--------------|---------------|
| project_id | project_id |
| sample_id | sample_id |
| taxon_id | taxon_id |
| taxon_name | taxon_name |
| relative_abundance | relative_abundance |
| read_count | read_count |
| sequencing_type | sequencing_type |
| host_age | host_metadata.age |
| host_sex | host_metadata.sex |
| host_bmi | host_metadata.bmi |
| country | host_metadata.country |
| rank | rank |
| enterotype | enterotype |
| kegg_ko | functional_annotations.kegg_ko |
| eggnog | functional_annotations.eggnog |
| cog | functional_annotations.cog |
| ibd_status | ibd_status |
| coverage | sequencing_metadata.coverage |
| sequencing_depth | sequencing_metadata.sequencing_depth |
