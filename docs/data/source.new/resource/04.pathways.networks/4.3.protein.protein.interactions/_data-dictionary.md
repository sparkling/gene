# 4.3 Protein-Protein Interactions - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| Subcategory ID | 4.3 |
| Subcategory Name | Protein-Protein Interactions |
| Data Sources | BioGRID, IntAct, STRING |
| Schema ID | `https://gene.example.org/schemas/4.3-protein-protein-interactions.json` |

## Unified Fields

These fields are harmonized across all data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `interactor_a_id` | string | Required (1:1) | Identifier for first interacting protein | BioGRID, IntAct, STRING | BioGRID: `7157`, IntAct: `uniprotkb:P04637`, STRING: `9606.ENSP00000269305` |
| `interactor_b_id` | string | Required (1:1) | Identifier for second interacting protein | BioGRID, IntAct, STRING | BioGRID: `3845`, IntAct: `uniprotkb:Q00987`, STRING: `9606.ENSP00000357274` |
| `gene_symbol_a` | string | Optional (1:1) | Gene symbol for interactor A | BioGRID, IntAct, STRING | `TP53` |
| `gene_symbol_b` | string | Optional (1:1) | Gene symbol for interactor B | BioGRID, IntAct, STRING | `KRAS`, `MDM2` |
| `detection_method` | string | Optional (1:1) | Experimental method used to detect interaction | BioGRID, IntAct, STRING | `Two-hybrid`, `psi-mi:"MI:0018"(two hybrid)`, `experimental (escore)` |
| `interaction_type` | string | Optional (1:1) | Type of molecular interaction | BioGRID, IntAct, STRING | `physical`, `psi-mi:"MI:0915"(physical association)`, `functional` |
| `publication_ids` | array&lt;string&gt; | Optional (1:N) | PubMed IDs of supporting publications | BioGRID, IntAct, STRING | `["1535557"]`, `["pubmed:1535557"]` |
| `organism_a` | integer | Optional (1:1) | NCBI taxonomy ID for interactor A | BioGRID, IntAct, STRING | `9606` |
| `organism_b` | integer | Optional (1:1) | NCBI taxonomy ID for interactor B | BioGRID, IntAct, STRING | `9606` |
| `confidence_score` | float | Optional (1:1) | Confidence score for the interaction | BioGRID, IntAct, STRING | IntAct: `0.73`, STRING: `0.999` |

## BioGRID-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `biogrid_interaction_id` | integer | Optional (1:1) | Unique BioGRID interaction identifier | `103` |
| `biogrid_id_a` | integer | Optional (1:1) | BioGRID internal ID for interactor A | `112315` |
| `biogrid_id_b` | integer | Optional (1:1) | BioGRID internal ID for interactor B | `107485` |
| `systematic_name_a` | string | Optional (1:1) | Systematic gene name for A | `YDR420W` |
| `synonyms_a` | array&lt;string&gt; | Optional (1:N) | Alias names (pipe-separated in source) | `["p53", "LFS1"]` |
| `experimental_system` | string | Optional (1:1) | Specific experimental method | `Affinity Capture-MS` |
| `experimental_system_type` | enum | Optional (1:1) | Broad interaction type | `physical`, `genetic` |
| `throughput` | enum | Optional (1:1) | Experiment scale | `high`, `low` |
| `modification` | string | Optional (1:1) | Post-translational modification info | `Phosphorylation` |
| `ontology_term_ids` | array&lt;string&gt; | Optional (1:N) | GO terms (pipe-separated in source, TAB 3.0) | `["GO:0005515", "GO:0005634"]` |
| `swissprot_accession_a` | string | Optional (1:1) | UniProt Swiss-Prot accession for A | `P04637` |

## IntAct-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `intact_interaction_ac` | string | Optional (1:1) | IntAct interaction accession | `EBI-77734` |
| `intact_miscore` | float | Optional (1:1) | IntAct molecular interaction confidence score (0-1) | `0.73` |
| `expansion_method` | string | Optional (1:1) | Complex expansion method | `psi-mi:"MI:1060"(spoke)`, `matrix` |
| `biological_role_a` | string | Optional (1:1) | Biological role of participant A | `psi-mi:"MI:0499"(unspecified role)` |
| `experimental_role_a` | string | Optional (1:1) | Experimental role of participant A | `psi-mi:"MI:0496"(bait)` |
| `experimental_role_b` | string | Optional (1:1) | Experimental role of participant B | `psi-mi:"MI:0498"(prey)` |
| `interactor_type_a` | string | Optional (1:1) | Molecule type of interactor A | `psi-mi:"MI:0326"(protein)`, `protein`, `dna`, `rna`, `small molecule`, `complex`, `gene` |
| `feature_a` | array&lt;string&gt; | Optional (1:N) | Binding regions or modifications | `["binding region:23-56"]` |
| `stoichiometry_a` | integer | Optional (1:1) | Number of molecules of participant A | `2` |
| `negative` | boolean | Optional (1:1) | Flag for negative (non-)interaction | `false` |
| `host_organism` | string | Optional (1:1) | Experimental host organism | `taxid:9606(human)` |
| `parameters_interaction` | array&lt;string&gt; | Optional (1:N) | Kinetic parameters (Kd, Ki, etc.) | `["kd:2.0e-9"]` |
| `checksum_interaction` | string | Optional (1:1) | RIGID interaction checksum | `rigid:4hL...` |

## STRING-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `string_id_a` | string | Optional (1:1) | STRING protein identifier (taxid.ENSP) | `9606.ENSP00000269305` |
| `string_id_b` | string | Optional (1:1) | STRING protein identifier (taxid.ENSP) | `9606.ENSP00000357274` |
| `string_combined_score` | float | Optional (1:1) | Combined confidence score (0-1) | `0.999` |
| `string_nscore` | float | Optional (1:1) | Neighborhood score (gene proximity) | `0.0` |
| `string_fscore` | float | Optional (1:1) | Fusion score (gene fusion events) | `0.0` |
| `string_pscore` | float | Optional (1:1) | Phylogenetic profile score (co-occurrence) | `0.0` |
| `string_ascore` | float | Optional (1:1) | Co-expression score | `0.076` |
| `string_escore` | float | Optional (1:1) | Experimental score (lab evidence) | `0.994` |
| `string_dscore` | float | Optional (1:1) | Database score (curated sources) | `0.900` |
| `string_tscore` | float | Optional (1:1) | Text-mining score (literature co-mentions) | `0.969` |
| `network_type` | enum | Optional (1:1) | Type of network | `functional`, `physical` |

## PSI-MI Detection Methods

Common detection methods using PSI-MI controlled vocabulary.

| MI ID | Term | Description |
|-------|------|-------------|
| MI:0006 | anti bait coimmunoprecipitation | IP with antibody against bait |
| MI:0018 | two hybrid | Yeast two-hybrid |
| MI:0019 | coimmunoprecipitation | Co-IP |
| MI:0096 | pull down | Pull-down assay |
| MI:0114 | x-ray crystallography | Crystal structure |
| MI:0676 | tandem affinity purification | TAP tag purification |
| MI:0004 | affinity capture-MS | Mass spec pull-down |
| MI:1313 | proximity label-MS | BioID, APEX |

## PSI-MI Interaction Types

| MI ID | Term | Description |
|-------|------|-------------|
| MI:0407 | direct interaction | Direct physical contact |
| MI:0403 | colocalization | Same cellular location |
| MI:0915 | physical association | Physical contact |

## Source Metadata

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `_source.database` | string | Required | Source database name | `BioGRID`, `IntAct`, `STRING` |
| `_source.version` | string | Optional | Database version | `4.4.220` |
| `_source.access_date` | date | Optional | Date data was retrieved | `2026-01-24` |
| `_source.original_id` | string | Optional | Original identifier in source | `103` |

## Field Mappings by Source

### BioGRID Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `BioGRID_Interaction_ID` | `biogrid_interaction_id` |
| `BioGRID_ID_A` | `biogrid_id_a` |
| `BioGRID_ID_B` | `biogrid_id_b` |
| `Entrez_Gene_A` | `interactor_a_id` |
| `Entrez_Gene_B` | `interactor_b_id` |
| `Official_Symbol_A` | `gene_symbol_a` |
| `Official_Symbol_B` | `gene_symbol_b` |
| `Systematic_Name_A` | `systematic_name_a` |
| `Synonyms_A` | `synonyms_a` |
| `Experimental_System` | `experimental_system` |
| `Experimental_System_Type` | `experimental_system_type` |
| `Throughput` | `throughput` |
| `Modification` | `modification` |
| `Ontology_Term_IDs` | `ontology_term_ids` |
| `SWISS-PROT_Accession_A` | `swissprot_accession_a` |
| `Organism_ID_A` | `organism_a` |
| `Organism_ID_B` | `organism_b` |
| `Pubmed_ID` | `publication_ids` |

### IntAct Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `interactionAc` | `intact_interaction_ac` |
| `idA` | `interactor_a_id` |
| `idB` | `interactor_b_id` |
| `altA` | `gene_symbol_a` |
| `altB` | `gene_symbol_b` |
| `detMethod` | `detection_method` |
| `type` | `interaction_type` |
| `pubid` | `publication_ids` |
| `taxidA` | `organism_a` |
| `taxidB` | `organism_b` |
| `intact-miscore` | `intact_miscore` |
| `expansion` | `expansion_method` |
| `bioRoleA` | `biological_role_a` |
| `expRoleA` | `experimental_role_a` |
| `expRoleB` | `experimental_role_b` |
| `typeA` | `interactor_type_a` |
| `xrefA` | `feature_a` |
| `stoichiometryA` | `stoichiometry_a` |
| `negative` | `negative` |
| `hostOrganism` | `host_organism` |
| `parameters` | `parameters_interaction` |
| `checksum` | `checksum_interaction` |

### STRING Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `protein1` | `string_id_a` |
| `protein2` | `string_id_b` |
| `combined_score` | `string_combined_score` |
| `nscore` | `string_nscore` |
| `fscore` | `string_fscore` |
| `pscore` | `string_pscore` |
| `ascore` | `string_ascore` |
| `escore` | `string_escore` |
| `dscore` | `string_dscore` |
| `tscore` | `string_tscore` |
| `network_type` | `network_type` |
