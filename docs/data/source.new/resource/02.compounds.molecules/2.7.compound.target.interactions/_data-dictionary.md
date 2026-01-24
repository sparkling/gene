# 2.7 Compound-Target Interactions - Data Dictionary

## Overview

This data dictionary documents the unified schema for compound-target interaction data from four major databases: BindingDB, DGIdb, Guide to Pharmacology (GtoPdb), and Therapeutic Target Database (TTD).

**Subcategory ID:** 2.7
**Subcategory Name:** Compound-Target Interactions
**Data Sources:** BindingDB, DGIdb, GtoPdb, TTD

---

## Unified Fields

### Compound Identification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| compound_id | string/integer | 1:1, Optional | Primary identifier for compound/ligand | `50001234`, `CHEMBL25` | BindingDB: monomer_id, GtoPdb: gtopdb_ligand_id, TTD: ttd_drug_id |
| compound_name | string | 1:1, Required | Drug/compound name | `Aspirin`, `Ibuprofen` | BindingDB, DGIdb, GtoPdb, TTD |

### Compound Structure Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| smiles | string | 1:1, Optional | SMILES structure | `CC(=O)Oc1ccccc1C(=O)O` | BindingDB, GtoPdb |
| inchi_key | string | 1:1, Optional | InChI Key | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` | BindingDB, GtoPdb |

### Target Identification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| target_id | string | 1:1, Optional | Target identifier | `P23219`, `TTDT00001` | BindingDB, DGIdb, GtoPdb, TTD |
| target_name | string | 1:1, Required | Target protein/gene name | `Cyclooxygenase-1`, `EGFR` | BindingDB, DGIdb, GtoPdb, TTD |
| uniprot_id | string | 1:1, Optional | UniProt accession | `P23219`, `P35354` | BindingDB, GtoPdb, TTD |

### Interaction Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| interaction_types | string[] | 1:N, Optional | Types of interaction | `["inhibitor", "agonist", "antagonist", "blocker", "modulator"]` | BindingDB, DGIdb, GtoPdb, TTD |
| affinity_value | decimal | 1:1, Optional | Binding affinity measurement | `50`, `0.5`, `1000` | BindingDB, GtoPdb |
| affinity_type | string | 1:1, Optional | Type of affinity measurement | `IC50`, `Ki`, `Kd`, `EC50`, `pKi` | BindingDB, GtoPdb |
| affinity_unit | string | 1:1, Optional | Unit of measurement | `nM`, `uM`, `mM` | BindingDB |

---

## Source-Specific Fields

### BindingDB

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| monomer_id | integer | 1:1, Required | BindingDB ligand ID | `50001234` |
| bindingdb_id | integer | 1:1, Optional | BindingDB interaction record ID | `12345678` |
| affinity_relation | string | 1:1, Optional | Affinity relation qualifier | `=`, `<`, `>`, `~` |
| assay_ph | decimal | 1:1, Optional | Assay pH | `7.4` |
| assay_temperature | decimal | 1:1, Optional | Temperature (Celsius) | `25`, `37` |
| pdb_ids | string[] | 1:N, Optional | PDB structure IDs | `["1PRH", "3N86"]` |
| pmid | integer | 1:1, Optional | PubMed ID for source publication | `12345678` |

### DGIdb

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| concept_id | string | 1:1, Optional | Normalized gene/drug ID | `ncbigene:1956`, `chembl:CHEMBL25` |
| gene_categories | string[] | 1:N, Optional | Druggability categories | `["KINASE", "DRUGGABLE GENOME", "TUMOR SUPPRESSOR"]` |
| interaction_score | decimal | 1:1, Optional | Confidence score for interaction | `5.0`, `10.0` |
| dgidb_sources | string[] | 1:N, Optional | Source databases (40+ aggregated sources) | `["ChEMBL", "DrugBank", "PharmGKB"]` |
| approved | boolean | 1:1, Optional | FDA approval status | `true`, `false` |
| chembl_id | string | 1:1, Optional | ChEMBL identifier | `CHEMBL25` |
| entrez_id | integer | 1:1, Optional | NCBI Entrez Gene ID | `1956` |

### Guide to Pharmacology (GtoPdb)

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| gtopdb_target_id | integer | 1:1, Required | GtoPdb target ID | `1234` |
| gtopdb_ligand_id | integer | 1:1, Required | GtoPdb ligand ID | `5678` |
| systematic_name | string | 1:1, Optional | IUPHAR systematic name | `2-acetoxybenzoic acid` |
| target_family | string | 1:1, Optional | Target family classification | `GPCR`, `ion_channel`, `kinase`, `nuclear_receptor`, `enzyme`, `transporter` |
| selectivity | string | 1:1, Optional | Selectivity classification | `Selective`, `Non-selective` |
| is_endogenous | boolean | 1:1, Optional | Endogenous ligand flag | `true`, `false` |
| inn | string | 1:1, Optional | International nonproprietary name | `aspirin` |
| hgnc_symbol | string | 1:1, Optional | HGNC gene symbol | `PTGS1`, `EGFR` |

### Therapeutic Target Database (TTD)

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| ttd_target_id | string | 1:1, Required | TTD target ID (TTDT prefix) | `TTDT00001` |
| ttd_drug_id | string | 1:1, Required | TTD drug ID (6 characters) | `D00ABC` |
| target_validation | string | 1:1, Optional | Target validation status | `Successful target`, `Clinical trial target`, `Literature-reported target` |
| development_status | string | 1:1, Optional | Drug development status | `Approved`, `Phase III`, `Phase II`, `Phase I`, `Preclinical`, `Discontinued` |
| drug_class | string | 1:1, Optional | Drug modality class | `Small molecule`, `Antibody`, `Protein`, `Peptide`, `Oligonucleotide` |
| indication | string | 1:1, Optional | Therapeutic indication | `Inflammation`, `Cancer`, `Hypertension` |
| icd_11_code | string | 1:1, Optional | ICD-11 disease classification | `BA00`, `2A00` |
| pathway_associations | object[] | 1:N, Optional | KEGG, Reactome pathway links | `[{"pathway_id": "hsa04010", "source": "KEGG"}]` |

---

## Metadata Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| _source.database | string | 1:1, Required | Source database name | `BindingDB`, `DGIdb`, `GtoPdb`, `TTD` |
| _source.version | string | 1:1, Optional | Database version | `2024.1`, `4.2.0` |
| _source.access_date | date | 1:1, Optional | Date data was accessed | `2026-01-24` |
| _source.original_id | string | 1:1, Optional | Original identifier from source | - |

---

## Field Mappings by Source

### BindingDB to Unified Schema

| BindingDB Field | Unified Field |
|-----------------|---------------|
| monomer_id | monomer_id |
| ligand_name | compound_name |
| smiles | smiles |
| inchi_key | inchi_key |
| target_name | target_name |
| uniprot_id | uniprot_id |
| activity_type | affinity_type |
| activity_value | affinity_value |
| activity_unit | affinity_unit |
| activity_relation | affinity_relation |
| ph | assay_ph |
| temperature | assay_temperature |
| pdb_ids | pdb_ids |
| pmid | pmid |
| bindingdb_id | bindingdb_id |

### DGIdb to Unified Schema

| DGIdb Field | Unified Field |
|-------------|---------------|
| drug_name | compound_name |
| gene_name | target_name |
| concept_id | concept_id |
| gene_categories | gene_categories |
| interaction_types | interaction_types |
| interaction_score | interaction_score |
| sources | dgidb_sources |
| approved | approved |
| chembl_id | chembl_id |
| entrez_id | entrez_id |

### GtoPdb to Unified Schema

| GtoPdb Field | Unified Field |
|--------------|---------------|
| ligand_id | gtopdb_ligand_id |
| ligand_name | compound_name |
| smiles | smiles |
| inchi_key | inchi_key |
| target_id | gtopdb_target_id |
| target_name | target_name |
| uniprot_id | uniprot_id |
| systematic_name | systematic_name |
| target_family | target_family |
| action | interaction_types |
| selectivity | selectivity |
| endogenous | is_endogenous |
| inn | inn |
| hgnc_symbol | hgnc_symbol |
| affinity_type | affinity_type |
| affinity_value | affinity_value |

### TTD to Unified Schema

| TTD Field | Unified Field |
|-----------|---------------|
| drug_id | ttd_drug_id |
| drug_name | compound_name |
| target_id | ttd_target_id |
| target_name | target_name |
| uniprot_id | uniprot_id |
| target_validation | target_validation |
| development_status | development_status |
| drug_class | drug_class |
| indication | indication |
| icd_11_code | icd_11_code |
| pathway_associations | pathway_associations |

---

## Data Quality Notes

- **Required Fields:** `compound_name` and `target_name` are required across all sources
- **Affinity Data:** Quantitative binding data primarily from BindingDB and GtoPdb
- **Aggregation:** DGIdb aggregates from 40+ source databases
- **Clinical Data:** TTD provides drug development status and therapeutic indications
- **Structures:** PDB structures available from BindingDB for co-crystal complexes

### Interaction Types

| Type | Description |
|------|-------------|
| inhibitor | Reduces target activity |
| agonist | Activates target (mimics endogenous ligand) |
| antagonist | Blocks target activation (competitive with agonist) |
| partial agonist | Partially activates target |
| inverse agonist | Reduces constitutive target activity |
| blocker | Physically blocks target (e.g., ion channel) |
| modulator | Allosteric modification of target activity |
| activator | Increases target activity (non-competitive) |
| inducer | Increases target expression |
| substrate | Metabolized by target enzyme |

### Affinity Measurement Types

| Type | Description | Typical Range |
|------|-------------|---------------|
| IC50 | Half-maximal inhibitory concentration | nM - uM |
| Ki | Inhibition constant | nM - uM |
| Kd | Dissociation constant | nM - uM |
| EC50 | Half-maximal effective concentration | nM - uM |
| pKi | -log10(Ki) | 5-10 |
| pIC50 | -log10(IC50) | 5-10 |
| Kon | Association rate constant | M-1s-1 |
| Koff | Dissociation rate constant | s-1 |

### Target Families (GtoPdb)

| Family | Description | Example Targets |
|--------|-------------|-----------------|
| GPCR | G protein-coupled receptors | Adrenergic receptors, Dopamine receptors |
| ion_channel | Ion channels | Sodium channels, Potassium channels |
| kinase | Protein kinases | EGFR, BCR-ABL |
| nuclear_receptor | Nuclear hormone receptors | Estrogen receptor, Glucocorticoid receptor |
| enzyme | Enzymes (non-kinase) | COX-1, COX-2, ACE |
| transporter | Membrane transporters | P-glycoprotein, SERT |
| other | Other target types | - |

### DGIdb Gene Categories

| Category | Description |
|----------|-------------|
| KINASE | Protein kinase |
| DRUGGABLE GENOME | Genes with drug-like binding properties |
| TUMOR SUPPRESSOR | Tumor suppressor gene |
| ONCOGENE | Cancer-promoting gene |
| TRANSCRIPTION FACTOR | DNA-binding transcription factor |
| DRUG RESISTANCE | Associated with drug resistance |
| CLINICALLY ACTIONABLE | Target with clinical relevance |
| CELL SURFACE | Cell surface protein |
| G PROTEIN COUPLED RECEPTOR | GPCR family member |
| ION CHANNEL | Ion channel protein |

### TTD Target Validation Levels

| Level | Description |
|-------|-------------|
| Successful target | Target with approved drug |
| Clinical trial target | Target in clinical trials |
| Literature-reported target | Target with preclinical evidence |
| Patented target | Target mentioned in patents |

### TTD Development Status

| Status | Description |
|--------|-------------|
| Approved | FDA/EMA approved drug |
| Phase III | Phase 3 clinical trials |
| Phase II | Phase 2 clinical trials |
| Phase I | Phase 1 clinical trials |
| Preclinical | Preclinical development |
| Discovery | Early discovery stage |
| Discontinued | Development discontinued |
| Withdrawn | Withdrawn from market |
