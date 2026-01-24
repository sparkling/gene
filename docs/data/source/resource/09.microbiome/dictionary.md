# Category 09: Microbiome - Data Dictionary

## Overview

This data dictionary documents the unified data schema for Category 09: Microbiome, which encompasses 10 data sources across 3 subcategories focused on human microbiome research, body site-specific communities, and microbe-host interactions.

**Category ID:** 09
**Category Name:** Microbiome
**Total Sources:** 10
**Total Unique Fields:** 147
**Total Unified Fields:** 89
**Total Source-Specific Fields:** 58

## Subcategories

| ID | Name | Data Sources | Description |
|----|------|--------------|-------------|
| 9.1 | Gut Microbiome | GMrepo, gutMGene, HMP, MetaHIT | Gut microbiome composition, taxa-phenotype associations, and gene expression |
| 9.2 | Body Site Microbiomes | HMP, HOMD, mBodyMap | Multi-site microbiome profiling across oral, skin, and other body sites |
| 9.3 | Microbe-Host Interactions | gutMDisorder, MASI, VMH | Disease associations, metabolite signaling, and metabolic modeling |

## Cross-Subcategory Shared Fields

The following fields are present across all three subcategories:

| Field Pattern | Description | Present In |
|---------------|-------------|------------|
| taxon_id / microbe_taxid | NCBI Taxonomy identifier - universal microbial identifier | 9.1, 9.2, 9.3 |
| taxon_name / microbe_name | Scientific name of microbial organism | 9.1, 9.2, 9.3 |
| phylum | Major taxonomic division | 9.1, 9.2, 9.3 |
| gram_stain | Cell wall classification | 9.1, 9.2, 9.3 |
| pmid | PubMed literature reference | 9.1, 9.2, 9.3 |

## Semantic Categories

### Microbial Taxonomy
Fields for identifying and classifying microbial organisms.

| Field | Type | Description |
|-------|------|-------------|
| taxon_id | integer | NCBI Taxonomy identifier |
| taxon_name | string | Binomial/genus name |
| rank | string | Taxonomic rank level |
| phylum | string | Major taxonomic division |
| homt_id | string | HOMD-specific identifier |

### Sample Metadata
Fields describing sample and subject characteristics.

| Field | Type | Description |
|-------|------|-------------|
| sample_id | string | Unique sample identifier |
| subject_id | string | De-identified subject ID |
| host_age | integer | Subject age in years |
| host_sex | string | Biological sex |
| host_bmi | float | Body mass index |
| country | string | Geographic origin |
| health_status | string | Clinical status |
| sequencing_method | string | Profiling technology |

### Abundance Data
Quantitative abundance and expression measurements.

| Field | Type | Description |
|-------|------|-------------|
| relative_abundance | float | Proportion (0-1 scale) |
| read_count | integer | Raw sequencing counts |
| absolute_count | integer | Absolute read counts |
| fold_change | float | Expression ratio |
| effect_size | float | Statistical effect magnitude |

### Host-Microbe Interactions
Fields describing interactions between microbes and host.

| Field | Type | Description |
|-------|------|-------------|
| direction | enum | Change direction (increased/decreased/altered) |
| effect | enum | Molecular effect type |
| disease_name | string | Associated disease |
| target_name | string | Host molecular target |
| pathway | string | Biological pathway |
| metabolite_name | string | Metabolite involved |

### Body Site Classification
Anatomical location and body site classification.

| Field | Type | Description |
|-------|------|-------------|
| body_site | string | Specific sampling location |
| supersite | string | Major body region |
| site_code | string | Hierarchical site code |
| fma_body_site | string | FMA/UBERON ontology term |
| oral_sites | array | Oral cavity locations |
| tissue | string | Tissue type |

### Metabolic Modeling
Fields for genome-scale metabolic reconstruction.

| Field | Type | Description |
|-------|------|-------------|
| rxn_id | string | Reaction identifier |
| metabolite_id | string | Metabolite identifier |
| gpr | string | Gene-Protein-Reaction rule |
| subsystem | string | Metabolic subsystem |
| lb | float | Lower flux bound |
| ub | float | Upper flux bound |
| model_id | string | AGORA model ID |

### Diversity Metrics
Alpha diversity and richness measurements.

| Field | Type | Description |
|-------|------|-------------|
| shannon | float | Shannon diversity index |
| simpson | float | Simpson diversity index |
| chao1 | float | Chao1 richness estimate |
| observed_species | integer | Species count |
| evenness | float | Pielou evenness index |

## Common Ontologies

| Ontology | Usage | ID Format |
|----------|-------|-----------|
| NCBI Taxonomy | Microbial species identification | Integer TaxID |
| MeSH | Disease/phenotype standardization | D###### format |
| UBERON | Body site anatomy | UBERON:####### format |
| FMA | Anatomical terms | FMA ID |
| ICD-10 | Disease classification | Letter+digits |
| Disease Ontology | Disease standardization | DOID:###### format |
| ChEBI | Chemical entities | Integer ChEBI ID |
| KEGG | Pathways and compounds | C#####, K##### format |
| ENVO | Environment classification | ENVO:######## format |

## Common Cross-References

| Database | ID Field | Usage |
|----------|----------|-------|
| NCBI Taxonomy | taxon_id, microbe_taxid | Microbial species identification |
| NCBI SRA | run_id | Raw sequence data access |
| NCBI BioSample | sample_id | Sample metadata |
| NCBI Gene | gene_id | Host gene identification |
| PubMed | pmid | Literature evidence |
| KEGG | kegg_ko, kegg_id | Functional/compound annotation |
| PubChem | pubchem_cid | Compound identification |
| HMDB | hmdb_id | Metabolite data |
| ChEBI | chebi_id | Chemical ontology |
| UniProt | uniprot_id | Protein targets |
| BiGG Models | bigg_id | Model compatibility |

## Data Formats

| Category | Supported Formats |
|----------|-------------------|
| Primary Formats | JSON, TSV, FASTA, SBML, BIOM |
| Encoding | UTF-8 |
| Compression | gzip |
| API Response | JSON |

## Related Dictionaries

- [9.1 Gut Microbiome Data Dictionary](./9.1.gut.microbiome/_data-dictionary.md)
- [9.2 Body Site Microbiomes Data Dictionary](./9.2.body.site.microbiomes/_data-dictionary.md)
- [9.3 Microbe-Host Interactions Data Dictionary](./9.3.microbe.host.interactions/_data-dictionary.md)
