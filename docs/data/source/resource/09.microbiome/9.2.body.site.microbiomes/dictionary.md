# 9.2 Body Site Microbiomes - Data Dictionary

## Overview

This data dictionary documents the unified schema for body site microbiome data, integrating HMP, HOMD (Human Oral Microbiome Database), and mBodyMap data sources.

**Subcategory ID:** 9.2
**Subcategory Name:** Body Site Microbiomes
**Data Sources:** HMP, HOMD, mBodyMap

## Required Fields

| Field | Type | Description | Sources |
|-------|------|-------------|---------|
| sample_id | string | Unique sample identifier | HMP, mBodyMap |
| body_site | string | Specific anatomical location for microbiome sampling | HMP, HOMD, mBodyMap |
| taxon_id | integer | NCBI Taxonomy identifier for microbial species | HMP, HOMD, mBodyMap |
| taxon_name | string | Scientific name of the microbial taxon | HMP, HOMD, mBodyMap |

## Unified Fields

### Sample and Subject Identifiers

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| sample_id | string | Required | Unique sample identifier | `SRS011084` | HMP, mBodyMap |
| subject_id | string | Optional | Subject/participant identifier | `HMP_SUBJ_001` | HMP, mBodyMap |

### Body Site Information

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| body_site | string | Required | Specific anatomical location for microbiome sampling | `Tongue dorsum` | HMP, HOMD, mBodyMap |
| supersite | enum | Optional | Major body region grouping multiple body sites. Values: `Oral`, `Gastrointestinal`, `Skin`, `Airways`, `Urogenital` | `Oral` | HMP, mBodyMap |
| site_code | string | Optional | Hierarchical site code for body location | `ORAL_TONGUE_DORSUM` | mBodyMap |
| fma_body_site | string | Optional | Foundational Model of Anatomy ontology term for body site | `UBERON:0001988` | HMP |

### Taxonomic Information

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| taxon_id | integer | Required | NCBI Taxonomy identifier for microbial species | `1309` | HMP, HOMD, mBodyMap |
| taxon_name | string | Required | Scientific name of the microbial taxon | `Streptococcus mutans` | HMP, HOMD, mBodyMap |
| homt_id | string | Optional | Human Oral Microbiome Taxon identifier | `HMT-096` | HOMD |
| cultivation_status | enum | Optional | Cultivation status of oral taxon. Values: `Named`, `Unnamed`, `Provisional` | `Named` | HOMD |

### Taxonomy Classification (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| taxonomy.phylum | string | Optional | Phylum classification | `Firmicutes` | HOMD, mBodyMap |
| taxonomy.class | string | Optional | Class classification | `Bacilli` | HOMD, mBodyMap |
| taxonomy.order | string | Optional | Order classification | `Lactobacillales` | HOMD, mBodyMap |
| taxonomy.family | string | Optional | Family classification | `Streptococcaceae` | HOMD, mBodyMap |
| taxonomy.genus | string | Optional | Genus classification | `Streptococcus` | HOMD, mBodyMap |
| taxonomy.species | string | Optional | Species classification | `mutans` | HOMD, mBodyMap |

### Abundance Data

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| relative_abundance | float | Optional | Proportion of taxon in sample (0-1 scale) | `0.15` | HMP, mBodyMap |
| absolute_count | integer | Optional | Absolute read count for taxon | `5000` | mBodyMap |

### Diversity Metrics (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| diversity_metrics.shannon | float | Optional | Shannon diversity index | `3.8` | mBodyMap |
| diversity_metrics.simpson | float | Optional | Simpson diversity index | `0.92` | mBodyMap |
| diversity_metrics.chao1 | float | Optional | Chao1 richness estimate | `250.5` | mBodyMap |
| diversity_metrics.observed_species | integer | Optional | Number of species detected | `185` | mBodyMap |
| diversity_metrics.evenness | float | Optional | Pielou evenness index | `0.85` | mBodyMap |

### Health and Sequencing

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| health_status | string | Optional | Health or disease status of subject | `healthy` | HMP, mBodyMap |
| sequencing_method | string | Optional | Sequencing technology used | `16S_V4` | HMP, mBodyMap |

### Microbial Phenotype (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| microbial_phenotype.gram_stain | enum | Optional | Gram staining classification. Values: `positive`, `negative`, `variable` | `positive` | HOMD |
| microbial_phenotype.oxygen_req | enum | Optional | Oxygen requirements for growth. Values: `aerobic`, `anaerobic`, `facultative` | `facultative` | HOMD |
| microbial_phenotype.cell_shape | string | Optional | Morphological cell shape | `coccus` | HOMD |
| microbial_phenotype.motility | boolean | Optional | Whether organism is motile | `false` | HOMD |
| microbial_phenotype.spore_forming | boolean | Optional | Ability to form spores | `false` | HOMD |
| microbial_phenotype.growth_temp | string | Optional | Optimal growth temperature | `37C` | HOMD |

### Genome Information (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| genome_info.genome_id | string | Optional | Reference genome identifier | `HOMD_GEN_096_001` | HMP, HOMD |
| genome_info.genome_size | integer | Optional | Genome size in base pairs | `2032807` | HMP, HOMD |
| genome_info.gc_content | float | Optional | GC percentage of genome | `36.8` | HMP, HOMD |
| genome_info.assembly_level | string | Optional | Completeness level of genome assembly | `Complete` | HMP, HOMD |
| genome_info.gene_count | integer | Optional | Number of predicted genes in genome | `1963` | HMP, HOMD |
| genome_info.ncbi_assembly | string | Optional | NCBI Assembly accession | `GCF_000007465.2` | HOMD |

### Strain and Disease Information

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| type_strain | string | Optional | Type strain designation for species definition | `ATCC 25175` | HOMD |
| disease_associations | array[string] | Optional | Associated diseases or conditions | `["dental caries", "endocarditis"]` | HOMD |
| oral_sites | array[string] | Optional | Oral cavity locations where taxon is found | `["teeth", "plaque", "saliva"]` | HOMD |

### Sequence Information (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| sequence_info.sequence_id | string | Optional | 16S rRNA sequence identifier | `SEQ001` | HOMD |
| sequence_info.clone_id | string | Optional | Clone identifier for sequence | `CLONE001` | HOMD |
| sequence_info.genbank_acc | string | Optional | GenBank accession for sequence | `AB123456` | HOMD |

### HMP-Specific Metadata (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| hmp_specific.rand_subject_id | string | Optional | Randomized HIPAA-compliant subject ID | `HMP_001` | HMP |
| hmp_specific.variable_region | string | Optional | 16S rRNA variable region sequenced | `V1-V3` | HMP |
| hmp_specific.lib_layout | string | Optional | Library layout (paired or fragment) | `paired` | HMP |
| hmp_specific.insert_size | integer | Optional | Library insert size | `300` | HMP |
| hmp_specific.matrix_type | string | Optional | Abundance matrix type classification | `16s_community` | HMP |

### mBodyMap-Specific Metadata (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| mbodymap_specific.site_category | string | Optional | Major body region category | `Oral` | mBodyMap |
| mbodymap_specific.anatomy_term | string | Optional | Anatomical ontology term (UBERON) | `UBERON:0001088` | mBodyMap |
| mbodymap_specific.sample_count | integer | Optional | Number of samples for body site | `500` | mBodyMap |
| mbodymap_specific.study_count | integer | Optional | Number of studies for body site | `15` | mBodyMap |
| mbodymap_specific.study_id | string | Optional | Source study identifier | `STUDY001` | mBodyMap |
| mbodymap_specific.dominant_phyla | array[string] | Optional | Most abundant phyla at body site | `["Firmicutes", "Bacteroidetes"]` | mBodyMap |
| mbodymap_specific.typical_shannon | float | Optional | Typical Shannon diversity for site | `3.5` | mBodyMap |
| mbodymap_specific.typical_richness | integer | Optional | Typical species richness for site | `150` | mBodyMap |

### Source Metadata

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| _source.primary_source | string | Optional | Name of the data source | `HMP` | All |
| _source.source_id | string | Optional | Original identifier in source | `HMP_001` | All |
| _source.extraction_date | string (date) | Optional | Date of data extraction | `2026-01-24` | All |
| _source.source_version | string | Optional | Version of the source database | `2.0` | All |

## Source-Specific Fields

### HMP-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| rand_subject_id | string | Randomized HIPAA-compliant subject ID | `HMP_001` |
| variable_region | string | 16S rRNA variable region sequenced | `V1-V3`, `V3-V5`, `V4` |
| lib_layout | string | Library layout (paired or fragment) | `paired` |
| insert_size | integer | Library insert size | `300` |
| matrix_type | enum | Abundance matrix type classification | `16s_community` |

### HOMD-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| sequence_id | string | 16S rRNA sequence identifier | `SEQ001` |
| clone_id | string | Clone identifier for sequence | `CLONE001` |
| genbank_acc | string | GenBank accession for sequence | `AB123456` |
| ncbi_assembly | string | NCBI Assembly accession | `GCF_000007465.2` |
| motility | boolean | Whether organism is motile | `true`, `false` |
| spore_forming | boolean | Ability to form spores | `true`, `false` |
| growth_temp | string | Optimal growth temperature | `37C` |

### mBodyMap-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| site_category | string | Major body region category | `Oral` |
| anatomy_term | string | Anatomical ontology term (UBERON) | `UBERON:0001088` |
| sample_count | integer | Number of samples for body site | `500` |
| study_count | integer | Number of studies for body site | `15` |
| study_id | string | Source study identifier | `STUDY001` |
| dominant_phyla | array | Most abundant phyla at body site | `["Firmicutes", "Bacteroidetes"]` |
| typical_shannon | float | Typical Shannon diversity for site | `3.5` |
| typical_richness | integer | Typical species richness for site | `150` |
| chao1 | float | Chao1 richness estimate | `250.5` |
| evenness | float | Pielou evenness index | `0.85` |
| absolute_count | integer | Absolute read count for taxon | `5000` |

## Cross-References

| Database | ID Field | Usage |
|----------|----------|-------|
| NCBI Taxonomy | taxon_id | Species identification |
| UBERON | fma_body_site | Body site ontology |
| FMA | fma_body_site | Anatomical terms |
| NCBI Assembly | ncbi_assembly | Genome links |
| GenBank | genbank_acc | Sequence links |
| ATCC/DSMZ | type_strain | Type strain references |
| SILVA | sequence_id | 16S database |

## Field Mappings by Source

### HMP Mappings

| Source Field | Unified Field |
|--------------|---------------|
| sample_id | sample_id |
| rand_subject_id | subject_id |
| body_site | body_site |
| body_supersite | supersite |
| fma_body_site | fma_body_site |
| taxon_id | taxon_id |
| taxon_name | taxon_name |
| relative_abundance | relative_abundance |
| health_status | health_status |
| sequencing_type | sequencing_method |
| variable_region | hmp_specific.variable_region |
| lib_layout | hmp_specific.lib_layout |
| insert_size | hmp_specific.insert_size |
| matrix_type | hmp_specific.matrix_type |
| genome_id | genome_info.genome_id |
| genome_size | genome_info.genome_size |
| gc_content | genome_info.gc_content |
| assembly_level | genome_info.assembly_level |
| gene_count | genome_info.gene_count |

### HOMD Mappings

| Source Field | Unified Field |
|--------------|---------------|
| homt_id | homt_id |
| taxon_id | taxon_id |
| taxon_name | taxon_name |
| cultivation_status | cultivation_status |
| phylum | taxonomy.phylum |
| gram_stain | microbial_phenotype.gram_stain |
| oxygen_req | microbial_phenotype.oxygen_req |
| cell_shape | microbial_phenotype.cell_shape |
| motility | microbial_phenotype.motility |
| spore_forming | microbial_phenotype.spore_forming |
| growth_temp | microbial_phenotype.growth_temp |
| type_strain | type_strain |
| disease_association | disease_associations |
| oral_sites | oral_sites |
| sequence_id | sequence_info.sequence_id |
| clone_id | sequence_info.clone_id |
| genbank_acc | sequence_info.genbank_acc |
| ncbi_assembly | genome_info.ncbi_assembly |
| genome_size | genome_info.genome_size |
| gc_content | genome_info.gc_content |
| gene_count | genome_info.gene_count |

### mBodyMap Mappings

| Source Field | Unified Field |
|--------------|---------------|
| sample_id | sample_id |
| subject_id | subject_id |
| body_site | body_site |
| supersite | supersite |
| site_code | site_code |
| taxon_id | taxon_id |
| taxon_name | taxon_name |
| relative_abundance | relative_abundance |
| absolute_count | absolute_count |
| shannon | diversity_metrics.shannon |
| simpson | diversity_metrics.simpson |
| chao1 | diversity_metrics.chao1 |
| observed_species | diversity_metrics.observed_species |
| evenness | diversity_metrics.evenness |
| health_status | health_status |
| sequencing_method | sequencing_method |
| phylum | taxonomy.phylum |
| site_category | mbodymap_specific.site_category |
| anatomy_term | mbodymap_specific.anatomy_term |
| sample_count | mbodymap_specific.sample_count |
| study_count | mbodymap_specific.study_count |
| study_id | mbodymap_specific.study_id |
| dominant_phyla | mbodymap_specific.dominant_phyla |
| typical_shannon | mbodymap_specific.typical_shannon |
| typical_richness | mbodymap_specific.typical_richness |

## Body Site Reference

### HMP Body Supersites

| Supersite | Example Body Sites |
|-----------|-------------------|
| Oral | Tongue dorsum, Buccal mucosa, Hard palate, Saliva, Subgingival plaque |
| Gastrointestinal | Stool |
| Skin | Retroauricular crease, Antecubital fossa, Volar forearm |
| Airways | Anterior nares |
| Urogenital | Posterior fornix, Mid vagina, Vaginal introitus |

### HOMD Oral Sites

| Site Category | Examples |
|---------------|----------|
| Teeth | Supragingival plaque, Subgingival plaque |
| Soft tissues | Tongue, Buccal mucosa, Hard palate |
| Fluids | Saliva |
