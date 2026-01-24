# 5.1 Traditional Chinese Medicine - Data Dictionary

## Overview

This data dictionary documents the unified schema for Traditional Chinese Medicine (TCM) data, integrating information from six major TCM databases.

**Subcategory ID:** 5.1
**Subcategory Name:** Traditional Chinese Medicine
**Data Sources:** HERB, ETCM, TCMBank, TCMSID, SymMap, BATMAN-TCM

---

## Data Sources Summary

| Database | Focus Area | Key Contributions |
|----------|------------|-------------------|
| HERB | High-throughput Experiment and Reference Database | Expression experiments, connectivity scores |
| ETCM | Encyclopedia of Traditional Chinese Medicine | Formulas, therapeutic actions, indications |
| TCMBank | TCM Bank Database | Drug-likeness, ADMET properties |
| TCMSID | TCM Systems Pharmacology Database | OB/DL screening, KEGG pathways |
| SymMap | Symptom Mapping Database | Symptoms, HPO mappings, patterns |
| BATMAN-TCM | Bioinformatics Analysis Tool for TCM | Target predictions, interaction evidence |

---

## Unified Fields

### Herb Identification

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| herb_id | string | 1:1 | Yes | Unique identifier for a TCM herb | All 6 sources | HERB001234, H001, SMHB00001 |
| herb_name_chinese | string | 1:1 | No | Chinese character name of the herb | All 6 sources | - |
| herb_name_pinyin | string | 1:1 | No | Romanized pronunciation using Pinyin | All 6 sources | Huang Qi, Jiang Huang, Gan Cao |
| herb_name_english | string | 1:1 | No | English common name | All 6 sources | - |
| herb_name_latin | string | 1:1 | No | Botanical Latin binomial name | All 6 sources | Astragalus membranaceus, Panax ginseng |

### TCM Properties

| Field Name | Data Type | Cardinality | Required | Description | Sources | Allowed Values |
|------------|-----------|-------------|----------|-------------|---------|----------------|
| tcm_nature | enum | 1:1 | No | TCM thermal nature classification (Xing/性) | HERB, ETCM, SymMap, BATMAN-TCM | cold, cool, neutral, warm, hot |
| tcm_flavor | array[string] | 1:N | No | TCM flavor classification (Wei/味) | HERB, ETCM, SymMap, BATMAN-TCM | sweet, bitter, sour, salty, pungent, astringent |
| meridian_tropism | array[string] | 1:N | No | TCM organ/meridian affinity (Gui Jing/归经) | HERB, ETCM, SymMap, BATMAN-TCM | lung, heart, liver, spleen, kidney, stomach, gallbladder, bladder, large_intestine, small_intestine, triple_burner, pericardium |

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| compound_id | string | 1:1 | No | Unique identifier for chemical compound | All 6 sources | C001, TCM001234, MOL000001, SMCP00001 |
| compound_name | string | 1:1 | No | Chemical compound name | All 6 sources | - |
| pubchem_cid | integer | 1:1 | No | PubChem Compound Identifier | All 6 sources | 5280343, 969516 |
| smiles | string | 1:1 | No | SMILES molecular notation | All 6 sources | - |
| inchi_key | string | 1:1 | No | InChIKey hash (27 characters) | HERB, TCMBank, BATMAN-TCM | - |
| molecular_formula | string | 1:1 | No | Chemical formula | All 6 sources | C21H20O6, C28H38O6 |
| molecular_weight | float | 1:1 | No | Molecular weight in Daltons | All 6 sources | - |
| source_herbs | array[string] | 1:N | No | Herb IDs containing this compound | All 6 sources | - |

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| target_id | string | 1:1 | No | Identifier for target protein/gene | All 6 sources | - |
| gene_symbol | string | 1:1 | No | HGNC-approved gene symbol | All 6 sources | TP53, EGFR, MTOR |
| uniprot_id | string | 1:1 | No | UniProt protein accession | All 6 sources | P04637, P00533 |
| entrez_gene_id | integer | 1:1 | No | NCBI Entrez Gene identifier | HERB, SymMap, BATMAN-TCM | - |
| predicted_targets | array[string] | 1:N | No | Computationally predicted target IDs | ETCM, TCMBank, TCMSID, BATMAN-TCM | - |

### Disease Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| disease_id | string | 1:1 | No | Identifier for disease entity | All 6 sources | - |
| disease_name | string | 1:1 | No | Disease or condition name | All 6 sources | - |
| mesh_id | string | 1:1 | No | Medical Subject Headings identifier | HERB, ETCM, TCMBank, SymMap, BATMAN-TCM | - |

### Formula Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| formula_id | string | 1:1 | No | Unique identifier for TCM formula | ETCM, TCMBank, BATMAN-TCM | - |
| formula_name_chinese | string | 1:1 | No | Chinese character formula name | ETCM, TCMBank, BATMAN-TCM | - |
| formula_name_pinyin | string | 1:1 | No | Romanized formula name | ETCM, TCMBank, BATMAN-TCM | Bu Zhong Yi Qi Tang, Liu Wei Di Huang Wan |
| classical_source | string | 1:1 | No | Historical text reference | ETCM, TCMBank, BATMAN-TCM | - |
| herb_role | enum | 1:1 | No | Role in formula (Jun-Chen-Zuo-Shi/君臣佐使) | ETCM, BATMAN-TCM | jun (sovereign), chen (minister), zuo (assistant), shi (envoy) |

### Pathway Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| kegg_pathways | array[object] | 1:N | No | KEGG pathway associations | TCMSID, BATMAN-TCM | {"pathway_id": "hsa04151", "pathway_name": "PI3K-Akt signaling pathway"} |

---

## Source-Specific Fields

### HERB Database

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| herb_geo_accession | string | 1:1 | Gene Expression Omnibus accession number | GSE12345 |
| herb_expression_experiments | array[string] | 1:N | GEO accessions for herb expression studies | - |
| herb_connectivity_score | float | 1:1 | Signature connectivity score (-1 to +1) | - |
| herb_differential_genes | object | 1:1 | Up/down-regulated genes from treatment | {"upregulated": [...], "downregulated": [...]} |

### ETCM Database

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| etcm_actions | array[string] | 1:N | Therapeutic actions (TCM theory) | tonifies qi, clears heat |
| etcm_indications | array[string] | 1:N | Disease patterns/conditions treated | - |
| etcm_therapeutic_category | string | 1:1 | Formula classification by function | tonifying formulas, heat-clearing formulas |
| etcm_herb_dosage | string | 1:1 | Recommended dosage in formula | - |

### TCMBank Database

| Field Name | Data Type | Cardinality | Description | Examples/Values |
|------------|-----------|-------------|-------------|-----------------|
| tcmbank_chembl_id | string | 1:1 | ChEMBL compound identifier | - |
| tcmbank_cas_number | string | 1:1 | CAS registry number | - |
| tcmbank_logp | float | 1:1 | Partition coefficient (lipophilicity) | - |
| tcmbank_tpsa | float | 1:1 | Topological Polar Surface Area (A^2) | - |
| tcmbank_hydrogen_bond_donors | integer | 1:1 | H-bond donors (Lipinski <=5) | - |
| tcmbank_hydrogen_bond_acceptors | integer | 1:1 | H-bond acceptors (Lipinski <=10) | - |
| tcmbank_rotatable_bonds | integer | 1:1 | Rotatable bonds (Veber <=10) | - |
| tcmbank_lipinski_violations | integer | 1:1 | Lipinski Rule of Five violations (0-4) | - |
| tcmbank_pains_alerts | integer | 1:1 | PAINS structural alert count | - |
| tcmbank_gi_absorption | enum | 1:1 | Predicted GI absorption | High, Low |
| tcmbank_bbb_permeant | enum | 1:1 | Predicted BBB permeability | Yes, No |
| tcmbank_cyp_inhibitors | array[string] | 1:N | CYP450 enzymes inhibited | CYP1A2, CYP2C9, CYP2D6, CYP3A4 |
| tcmbank_ames_mutagenicity | enum | 1:1 | Predicted Ames mutagenicity | Yes, No |
| tcmbank_hepatotoxicity_risk | enum | 1:1 | Predicted hepatotoxicity risk | High, Medium, Low |

### TCMSID Database

| Field Name | Data Type | Cardinality | Description | Threshold |
|------------|-----------|-------------|-------------|-----------|
| tcmsid_oral_bioavailability | float | 1:1 | Oral Bioavailability percentage | >=30% favorable |
| tcmsid_drug_likeness | float | 1:1 | Drug-Likeness score (0-1) | >=0.18 drug-like |
| tcmsid_caco2 | float | 1:1 | Caco-2 cell permeability | >=-0.4 good |
| tcmsid_bbb | float | 1:1 | Blood-Brain Barrier score | >=-0.3 penetration |
| tcmsid_half_life | enum | 1:1 | Pharmacokinetic half-life | Long, Short |

### SymMap Database

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| symmap_symptom_id | string | 1:1 | Unique TCM symptom identifier | SMSY00001 |
| symmap_symptom_name_chinese | string | 1:1 | Chinese symptom name | - |
| symmap_symptom_name_english | string | 1:1 | English symptom translation | - |
| symmap_tcm_pattern | string | 1:1 | TCM syndrome/pattern (Zheng/证) | - |
| symmap_hpo_mappings | array[object] | 1:N | Human Phenotype Ontology mappings | {"hpo_id": "HP:0002315", "hpo_term": "Headache", "mapping_confidence": 0.95} |
| symmap_symptom_herb_associations | array[string] | 1:N | Herbs treating this symptom | - |
| symmap_symptom_gene_associations | array[string] | 1:N | Genes associated via HPO | - |
| symmap_omim_id | integer | 1:1 | OMIM disease/phenotype identifier | - |

### BATMAN-TCM Database

| Field Name | Data Type | Cardinality | Description | Values/Range |
|------------|-----------|-------------|-------------|--------------|
| batman_tti_id | string | 1:1 | Target-TCM Interaction identifier | - |
| batman_interaction_type | enum | 1:1 | Interaction classification | known, predicted |
| batman_confidence_score | float | 1:1 | Prediction confidence (0-1) | - |
| batman_evidence_source | string | 1:1 | Source of known interaction | ChEMBL, CTD, STITCH, DGIdb, literature |
| batman_binding_data | object | 1:1 | Quantitative binding affinity | {"activity_type": "IC50", "activity_value": 50, "activity_unit": "nM"} |
| batman_pubmed_ids | array[string] | 1:N | PubMed literature references | - |
| batman_stitch_score | integer | 1:1 | STITCH combined score (0-1000) | - |
| batman_enrichment_ratio | float | 1:1 | KEGG pathway enrichment ratio | - |
| batman_pharmacopoeia_status | enum | 1:1 | Chinese Pharmacopoeia status | CP2020, not_listed |
| batman_herb_processing | string | 1:1 | Processing method for preparation | - |
| batman_drugbank_id | string | 1:1 | DrugBank database identifier | - |

---

## Source Field Mappings

### HERB Mappings
| Original Field | Unified Field |
|----------------|---------------|
| herb_id | herb_id |
| Herb_cn_name | herb_name_chinese |
| Herb_pinyin_name | herb_name_pinyin |
| Herb_en_name | herb_name_english |
| Herb_latin_name | herb_name_latin |
| nature | tcm_nature |
| flavor | tcm_flavor |
| meridians | meridian_tropism |
| Ingredient_id | compound_id |
| Molecule_name | compound_name |
| PubChem_id | pubchem_cid |
| SMILES | smiles |
| InChIKey | inchi_key |
| Formula | molecular_formula |
| MW | molecular_weight |
| Target_id | target_id |
| Gene_name | gene_symbol |
| Uniprot_id | uniprot_id |
| Entrez_id | entrez_gene_id |
| Disease_id | disease_id |
| Disease_name | disease_name |
| MeSH_id | mesh_id |

### ETCM Mappings
| Original Field | Unified Field |
|----------------|---------------|
| Herb_id | herb_id |
| Herb_cn | herb_name_chinese |
| Herb_pinyin | herb_name_pinyin |
| Herb_en | herb_name_english |
| Herb_latin | herb_name_latin |
| Nature | tcm_nature |
| Flavor | tcm_flavor |
| Meridian | meridian_tropism |
| Formula_id | formula_id |
| Formula_cn | formula_name_chinese |
| Formula_pinyin | formula_name_pinyin |
| Classical_Source | classical_source |
| Herb_Role | herb_role |
| Actions | etcm_actions |
| Indications | etcm_indications |
| Therapeutic_Category | etcm_therapeutic_category |
| Herb_Dosage | etcm_herb_dosage |

---

## Metadata Field

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| _source.database | string | Source database name |
| _source.version | string | Database version |
| _source.access_date | date | Date of data access |
| _source.original_id | string | Original identifier in source database |
