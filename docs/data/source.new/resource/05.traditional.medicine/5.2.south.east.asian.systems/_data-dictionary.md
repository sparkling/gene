# 5.2 South/East Asian Systems - Data Dictionary

## Overview

This data dictionary documents the unified schema for South/East Asian traditional medicine systems, integrating data from four databases covering Ayurveda, Kampo, and multi-cultural Asian medicinal plant data.

**Subcategory ID:** 5.2
**Subcategory Name:** South/East Asian Systems
**Data Sources:** NPACT, TM-MC, IMPPAT, KampoDB

---

## Data Sources Summary

| Database | Focus Area | Key Contributions |
|----------|------------|-------------------|
| NPACT | Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target | Cancer cell line activity data |
| TM-MC | Traditional Medicine Multi-Culture | Cross-system validation, Ayurvedic properties |
| IMPPAT | Indian Medicinal Plants, Phytochemistry And Therapeutics | Indian plants, comprehensive ADMET, vernacular names |
| KampoDB | Kampo Database | Japanese Kampo formulas, docking-based target predictions |

---

## Unified Fields

### Plant Identification

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| plant_id | string | 1:1 | Yes | Unique identifier for medicinal plant | NPACT, TM-MC, IMPPAT | TMMC_P_001234, IMPPAT_PLANT_03985, NAP_ORG_045678 |
| botanical_name | string | 1:1 | No | Scientific Latin binomial name | All 4 sources | - |
| family | string | 1:1 | No | Botanical family classification | NPACT, TM-MC, IMPPAT | Zingiberaceae, Solanaceae |
| traditional_system | array[enum] | 1:N | No | Traditional medicine system of origin | TM-MC, IMPPAT | Ayurveda, TCM, Unani, Siddha, Kampo, Jamu, Sowa-Rigpa, Homeopathy, Western Herbal |
| plant_part | array[string] | 1:N | No | Plant parts used therapeutically | NPACT, TM-MC, IMPPAT | root, leaf, bark, seed, flower, fruit, stem, rhizome, tuber, whole plant |

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| compound_id | string | 1:1 | No | Unique compound identifier | All 4 sources | - |
| compound_name | string | 1:1 | No | Chemical compound name | All 4 sources | - |
| pubchem_cid | integer | 1:1 | No | PubChem Compound Identifier | All 4 sources | - |
| smiles | string | 1:1 | No | SMILES molecular notation | NPACT, TM-MC, IMPPAT | - |
| molecular_formula | string | 1:1 | No | Chemical formula | All 4 sources | - |
| molecular_weight | float | 1:1 | No | Molecular weight in Daltons | NPACT, TM-MC, IMPPAT | - |

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| target_id | string | 1:1 | No | Target protein identifier | All 4 sources | - |
| gene_symbol | string | 1:1 | No | HGNC gene symbol | All 4 sources | - |
| uniprot_id | string | 1:1 | No | UniProt protein accession | NPACT, TM-MC, IMPPAT | - |

---

## Source-Specific Fields

### NPACT Database (Cancer Activity)

| Field Name | Data Type | Cardinality | Description | Examples/Values |
|------------|-----------|-------------|-------------|-----------------|
| npact_id | string | 1:1 | NPACT internal compound identifier | NPACT001234 |
| npact_activity_id | string | 1:1 | Unique activity measurement identifier | - |
| npact_activity_type | enum | 1:1 | Type of bioactivity measurement | IC50, EC50, GI50, LC50, CC50 |
| npact_activity_value | float | 1:1 | Quantitative activity value | - |
| npact_activity_unit | string | 1:1 | Unit for activity measurement | nM, uM, mM, ug/mL |
| npact_cell_line_name | string | 1:1 | Cancer cell line tested | MCF-7, A549, HeLa, HT-29 |
| npact_cancer_type | string | 1:1 | Cancer type for cell line | Breast, Lung, Colon, Leukemia, Liver |
| npact_assay_type | string | 1:1 | Type of assay performed | - |
| npact_pubmed_id | integer | 1:1 | PubMed reference ID | - |
| npact_bioavailability_score | float | 1:1 | Predicted oral bioavailability score | - |

### TM-MC Database (Multi-Culture)

| Field Name | Data Type | Cardinality | Description | Examples/Values |
|------------|-----------|-------------|-------------|-----------------|
| tmmc_system_code | enum | 1:1 | Traditional system abbreviation | AYU, TCM, UNA, SID, KAM, JAM |
| tmmc_local_name | string | 1:1 | Plant name in local language | - |
| tmmc_convergent_uses | array[object] | 1:N | Therapeutic uses shared across systems | {"indication": "...", "systems_agreeing": [...], "confidence": "high/medium/low"} |
| tmmc_cross_system_validation | object | 1:1 | Cross-cultural validation metadata | {"systems_present": [...], "shared_compounds": int, "shared_uses": [...]} |
| tmmc_ayurveda_properties | object | 1:1 | Ayurvedic property classification | See structure below |
| tmmc_indication_id | string | 1:1 | Therapeutic indication identifier | - |
| tmmc_traditional_terms | array[object] | 1:N | Indication terminology by system | {"system": "...", "term": "..."} |

#### Ayurveda Properties Structure

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| rasa | array[string] | Taste classification |
| guna | array[string] | Physical qualities |
| vipaka | string | Post-digestive effect |
| virya | string | Potency (heating/cooling) |
| dosha_effect | array[string] | Effect on doshas (Vata/Pitta/Kapha) |

### IMPPAT Database (Indian Medicinal Plants)

#### Identifiers

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| imppat_plant_id | string | 1:1 | IMPPAT plant identifier | IMPPAT_PLANT_XXXXX |
| imppat_id | string | 1:1 | IMPPAT phytochemical identifier | IMPPAT_CHEM_XXXXX, IMPHY###### |
| imppat_iupac_name | string | 1:1 | Systematic IUPAC chemical name | - |
| imppat_inchi | string | 1:1 | Full InChI string | - |
| imppat_inchi_key | string | 1:1 | InChIKey hash (27 characters) | - |
| imppat_deepsmiles | string | 1:1 | DeepSMILES notation for ML | - |
| imppat_hgnc_id | string | 1:1 | HUGO Gene Nomenclature Committee ID | - |

#### Vernacular Names (10 Indian Languages)

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| imppat_vernacular_names.hindi | string | Hindi name |
| imppat_vernacular_names.tamil | string | Tamil name |
| imppat_vernacular_names.telugu | string | Telugu name |
| imppat_vernacular_names.kannada | string | Kannada name |
| imppat_vernacular_names.malayalam | string | Malayalam name |
| imppat_vernacular_names.bengali | string | Bengali name |
| imppat_vernacular_names.gujarati | string | Gujarati name |
| imppat_vernacular_names.marathi | string | Marathi name |
| imppat_vernacular_names.punjabi | string | Punjabi name |
| imppat_vernacular_names.oriya | string | Oriya name |

#### Plant Information

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| imppat_source_plant_parts | array[string] | 1:N | Plant parts where compound was isolated |
| imppat_iucn_conservation_status | string | 1:1 | IUCN Red List conservation status |
| imppat_therapeutic_use | array[string] | 1:N | Traditional therapeutic application |

#### Drug-Likeness Filters

| Field Name | Data Type | Description | Structure |
|------------|-----------|-------------|-----------|
| imppat_lipinski_rule_of_five | object | Lipinski RO5 assessment | {passes, violations, mw_ok, logp_ok, hbd_ok, hba_ok} |
| imppat_ghose_filter | object | Ghose filter compliance | {passes, violations} |
| imppat_veber_filter | object | Veber oral bioavailability | {passes, rotatable_bonds_ok, tpsa_ok} |
| imppat_qed_score | float | Quantitative Estimate of Drug-likeness (0-1) | - |
| imppat_np_likeness_score | float | Natural Product likeness score | - |

#### ADMET Properties

| Field Name | Data Type | Description | Structure |
|------------|-----------|-------------|-----------|
| imppat_admet_absorption | object | Absorption properties | {gi_absorption, pgp_substrate, water_solubility, caco2_permeability} |
| imppat_admet_distribution | object | Distribution properties | {bbb_permeant, vdss} |
| imppat_admet_metabolism | object | CYP450 inhibition predictions | {cyp1a2_inhibitor, cyp2c19_inhibitor, cyp2c9_inhibitor, cyp2d6_inhibitor, cyp3a4_inhibitor} |
| imppat_admet_toxicity | object | Toxicity predictions | {ames_mutagenicity, hepatotoxicity, herg_inhibition, skin_sensitization} |

#### Molecular Descriptors and Classification

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| imppat_descriptor_count | integer | 1:1 | Number of computed descriptors (up to 1,875) |
| imppat_combined_score | integer | 1:1 | STITCH combined confidence (0-1000, >=700 high confidence) |
| imppat_murcko_scaffold | object | 1:1 | Murcko scaffold decomposition |
| imppat_npclassifier | object | 1:1 | Natural product classification |

#### Murcko Scaffold Structure

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| generic_node_bond | string | Generic node-bond scaffold |
| generic_node | string | Generic node scaffold |
| graph | string | Graph representation |

#### NPClassifier Structure

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| biosynthetic_pathway | string | Biosynthetic pathway |
| superclass | string | NP superclass |
| class | string | NP class |
| subclass | string | NP subclass |

### KampoDB Database (Japanese Kampo)

#### Formula Information

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| kampo_formula_id | string | 1:1 | Kampo formula code identifier | KT, YKS, TSS |
| kampo_formula_name | string | 1:1 | Romanized Kampo formula name | Kakkonto, Yokukansankan |
| kampo_formula_name_jp | string | 1:1 | Japanese kanji formula name | 葛根湯, 抑肝散 |
| kampo_phonetic | string | 1:1 | Japanese kana pronunciation | かっこんとう |

#### Crude Drug Information

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| kampo_crude_id | integer | 1:1 | Crude drug identifier (0-179) | - |
| kampo_crude_name | string | 1:1 | English crude drug name | Cinnamon Bark, Ephedra Herb, Ginger |
| kampo_crude_name_jp | string | 1:1 | Japanese crude drug name | 桂皮, 麻黄, 生姜 |
| kampo_origin | string | 1:1 | Botanical/biological source | Cinnamomum cassia, Ephedra sinica |

#### Target and Docking Information

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| kampo_protein_id | integer | 1:1 | NCBI Gene ID for target protein | - |
| kampo_protein_name | string | 1:1 | Gene symbol/name | MTOR, CYP3A4, BCL2 |
| kampo_protein_aliases | string | 1:1 | Comma-separated alias names | - |
| kampo_domain_id | integer | 1:1 | Protein structural domain ID | - |
| kampo_affinity_kcal_mol | float | 1:1 | Predicted binding affinity (AutoDock Vina) | -7.8, -6.5 (more negative = stronger) |
| kampo_ligand_atoms | integer | 1:1 | Number of atoms in ligand | - |
| kampo_protein_atoms | integer | 1:1 | Number of atoms in protein domain | - |

#### Pathway and Annotation

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| kampo_kegg_pathway_id | string | 1:1 | KEGG pathway identifier | hsa04150, hsa04151 |
| kampo_kegg_disease_id | string | 1:1 | KEGG disease identifier | - |
| kampo_go_process | array[string] | 1:N | GO Biological Process terms | - |
| kampo_go_function | array[string] | 1:N | GO Molecular Function terms | - |

---

## Source Field Mappings

### NPACT Mappings
| Original Field | Unified Field |
|----------------|---------------|
| compound_id | npact_id |
| plant_id | plant_id |
| scientific_name | botanical_name |
| family | family |
| activity_id | npact_activity_id |
| activity_type | npact_activity_type |
| activity_value | npact_activity_value |
| activity_unit | npact_activity_unit |
| cell_line_name | npact_cell_line_name |
| cancer_type | npact_cancer_type |

### TM-MC Mappings
| Original Field | Unified Field |
|----------------|---------------|
| system_code | tmmc_system_code |
| local_name | tmmc_local_name |
| convergent_uses | tmmc_convergent_uses |
| cross_system_validation | tmmc_cross_system_validation |
| ayurveda_properties | tmmc_ayurveda_properties |
| indication_id | tmmc_indication_id |
| traditional_terms | tmmc_traditional_terms |

### IMPPAT Mappings
| Original Field | Unified Field |
|----------------|---------------|
| plant_id | imppat_plant_id |
| compound_id | imppat_id |
| vernacular_names | imppat_vernacular_names |
| iupac_name | imppat_iupac_name |
| inchi | imppat_inchi |
| inchi_key | imppat_inchi_key |
| deepsmiles | imppat_deepsmiles |
| lipinski_rule_of_five | imppat_lipinski_rule_of_five |
| qed_score | imppat_qed_score |
| admet_absorption | imppat_admet_absorption |
| admet_toxicity | imppat_admet_toxicity |

### KampoDB Mappings
| Original Field | Unified Field |
|----------------|---------------|
| formula_id | kampo_formula_id |
| formula_name | kampo_formula_name |
| formula_name_jp | kampo_formula_name_jp |
| phonetic | kampo_phonetic |
| crude_id | kampo_crude_id |
| crude_name | kampo_crude_name |
| affinity_kcal_mol | kampo_affinity_kcal_mol |
| kegg_pathway_id | kampo_kegg_pathway_id |

---

## Metadata Field

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| _source.database | string | Source database name |
| _source.version | string | Database version |
| _source.access_date | date | Date of data access |
| _source.original_id | string | Original identifier in source database |
