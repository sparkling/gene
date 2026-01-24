# 9.3 Microbe-Host Interactions - Data Dictionary

## Overview

This data dictionary documents the unified schema for microbe-host interaction data, integrating gutMDisorder, MASI (Microbiome-Associated Signaling Interactions), and VMH (Virtual Metabolic Human) data sources.

**Subcategory ID:** 9.3
**Subcategory Name:** Microbe-Host Interactions
**Data Sources:** gutMDisorder, MASI, VMH

## Required Fields

| Field | Type | Description | Sources |
|-------|------|-------------|---------|
| microbe_name | string | Scientific name of microbe | gutMDisorder, MASI, VMH |

## Unified Fields

### Interaction Identifiers

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| interaction_id | string | Optional | Unique identifier for the microbe-host interaction | `GMD_00001` | gutMDisorder, MASI |
| microbe_taxid | integer | Optional | NCBI Taxonomy identifier for microbe | `816` | gutMDisorder, VMH |
| microbe_name | string | Required | Scientific name of microbe | `Bacteroides vulgatus` | gutMDisorder, MASI, VMH |
| microbe_rank | enum | Optional | Taxonomic rank of microbe. Values: `species`, `genus`, `family`, `order`, `phylum` | `species` | gutMDisorder |

### Disease Association (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| disease_association.disease_name | string | Optional | Name of disease or disorder | `Type 2 Diabetes` | gutMDisorder, VMH |
| disease_association.disease_category | string | Optional | Broad disease category classification | `Metabolic disorders` | gutMDisorder |
| disease_association.icd10_code | string | Optional | ICD-10 disease classification code | `E11` | gutMDisorder |
| disease_association.mesh_id | string | Optional | Medical Subject Headings identifier | `D003924` | gutMDisorder |
| disease_association.doid | string | Optional | Disease Ontology identifier | `DOID:9352` | gutMDisorder |
| disease_association.direction | enum | Optional | Direction of change in disease context. Values: `increased`, `decreased`, `altered` | `increased` | gutMDisorder |

### Statistical Evidence (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| statistical_evidence.effect_size | float | Optional | Magnitude of association | `0.45` | gutMDisorder |
| statistical_evidence.p_value | float | Optional | Statistical significance | `0.001` | gutMDisorder |
| statistical_evidence.sample_size | integer | Optional | Total number of subjects in study | `200` | gutMDisorder |
| statistical_evidence.cases | integer | Optional | Number of disease cases in study | `100` | gutMDisorder |
| statistical_evidence.controls | integer | Optional | Number of healthy controls in study | `100` | gutMDisorder |

### Study Metadata (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| study_metadata.study_design | string | Optional | Type of study design | `case-control` | gutMDisorder |
| study_metadata.method | string | Optional | Detection or sequencing method | `16S rRNA sequencing` | gutMDisorder |
| study_metadata.pmid | integer | Optional | PubMed identifier | `28123456` | gutMDisorder |

### Metabolite Signaling (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| metabolite_signaling.metabolite_id | string | Optional | Metabolite identifier | `but[e]` | MASI |
| metabolite_signaling.metabolite_name | string | Optional | Name of metabolite | `Butyrate` | MASI |
| metabolite_signaling.metabolite_class | string | Optional | Chemical class of metabolite | `Short-chain fatty acid` | MASI |
| metabolite_signaling.microbe_source | string | Optional | Microbe(s) producing the metabolite | `Faecalibacterium prausnitzii` | MASI |
| metabolite_signaling.pubchem_cid | integer | Optional | PubChem compound identifier | `264` | MASI |
| metabolite_signaling.hmdb_id | string | Optional | Human Metabolome Database identifier | `HMDB0000039` | MASI |
| metabolite_signaling.chebi_id | string | Optional | ChEBI identifier | `17634` | MASI |
| metabolite_signaling.smiles | string | Optional | SMILES structure notation | `CCCC(=O)[O-]` | MASI |

### Host Target (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| host_target.target_name | string | Optional | Host target protein/receptor name | `GPR43` | MASI |
| host_target.target_type | string | Optional | Type of host target | `GPCR` | MASI |
| host_target.symbol | string | Optional | Gene/protein symbol for target | `FFAR2` | MASI |
| host_target.uniprot_id | string | Optional | UniProt protein accession | `O15552` | MASI |
| host_target.effect | enum | Optional | Effect type on host target. Values: `activate`, `inhibit`, `modulate`, `substrate` | `activate` | MASI |
| host_target.tissue_expression | array[string] | Optional | Tissues where target is expressed | `["colon", "adipose tissue"]` | MASI |

### Pathway Effects (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| pathway_effects.pathway | string | Optional | Signaling or metabolic pathway | `Prolactin signaling` | MASI, VMH |
| pathway_effects.pathway_id | string | Optional | KEGG or Reactome pathway ID | `hsa04917` | MASI |
| pathway_effects.pathway_category | string | Optional | Pathway classification category | `Signal transduction` | MASI |
| pathway_effects.downstream_effects | array[string] | Optional | Biological outcomes of interaction | `["anti-inflammatory", "gut barrier"]` | MASI |

### Binding Affinity (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| binding_affinity.kd | string | Optional | Dissociation constant for binding | `1.2 uM` | MASI |
| binding_affinity.mechanism | string | Optional | Molecular mechanism description | `Competitive agonist` | MASI |

### Metabolic Model (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| metabolic_model.model_id | string | Optional | Metabolic model identifier | `Bacteroides_thetaiotaomicron_VPI_5482` | VMH |
| metabolic_model.organism_name | string | Optional | Full organism name for model | `Bacteroides thetaiotaomicron` | VMH |
| metabolic_model.strain | string | Optional | Strain designation | `VPI-5482` | VMH |
| metabolic_model.reactions_count | integer | Optional | Number of reactions in model | `1498` | VMH |
| metabolic_model.metabolites_count | integer | Optional | Number of metabolites in model | `1371` | VMH |
| metabolic_model.genes_count | integer | Optional | Number of genes in model | `872` | VMH |

### Metabolite Information (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| metabolite_info.metabolite_id | string | Optional | Metabolite identifier | `but[e]` | VMH |
| metabolite_info.metabolite_name | string | Optional | Name of metabolite | `Butyrate` | VMH |
| metabolite_info.formula | string | Optional | Chemical formula | `C4H7O2` | VMH |
| metabolite_info.charge | integer | Optional | Molecular charge of metabolite | `-1` | VMH |
| metabolite_info.compartment | string | Optional | Cellular compartment for metabolite | `e` | VMH |
| metabolite_info.kegg_id | string | Optional | KEGG compound identifier | `C00246` | VMH |
| metabolite_info.hmdb_id | string | Optional | Human Metabolome Database identifier | `HMDB0000039` | VMH |
| metabolite_info.chebi_id | string | Optional | ChEBI identifier | `17634` | VMH |
| metabolite_info.pubchem_cid | integer | Optional | PubChem compound identifier | `264` | VMH |
| metabolite_info.inchi | string | Optional | InChI structure representation | `InChI=1S/C4H8O2/...` | VMH |

### Reaction Information (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| reaction_info.rxn_id | string | Optional | Metabolic reaction identifier | `HEX1` | VMH |
| reaction_info.equation | string | Optional | Human-readable reaction equation | `atp[c] + glc_D[c] -> adp[c] + g6p[c] + h[c]` | VMH |
| reaction_info.subsystem | array[string] | Optional | Metabolic subsystem classification | `["Glycolysis/gluconeogenesis"]` | VMH |
| reaction_info.gpr | string | Optional | Gene-Protein-Reaction association rule | `(HKDC1) or (GCK)` | VMH |
| reaction_info.lb | float | Optional | Lower bound for reaction flux | `0` | VMH |
| reaction_info.ub | float | Optional | Upper bound for reaction flux | `1000` | VMH |
| reaction_info.reversibility | boolean | Optional | Whether reaction is reversible | `true` | VMH |
| reaction_info.ec_number | string | Optional | Enzyme Commission classification | `2.7.1.1` | VMH |

### Microbe Phenotype (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| microbe_phenotype.gram_stain | string | Optional | Gram staining classification | `negative` | VMH |
| microbe_phenotype.oxygen_status | string | Optional | Oxygen requirement status | `anaerobic` | VMH |

### Cross-References (Object)

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| cross_references.uniprot_id | string | Optional | UniProt protein accession | `O15552` | MASI, VMH |
| cross_references.entrez_id | integer | Optional | NCBI Entrez Gene ID | `12345` | VMH |
| cross_references.bigg_id | string | Optional | BiGG Models database identifier | `but` | VMH |

### Literature Evidence

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| pmid | array[integer] | Optional | PubMed identifiers for literature evidence | `[28123456]` | gutMDisorder, MASI |

### Source Metadata

| Field | Type | Required | Description | Example | Sources |
|-------|------|----------|-------------|---------|---------|
| _source.primary_source | string | Optional | Name of the data source | `gutMDisorder` | All |
| _source.source_id | string | Optional | Original identifier in source | `GMD_001` | All |
| _source.extraction_date | string (date) | Optional | Date of data extraction | `2026-01-24` | All |
| _source.source_version | string | Optional | Version of the source database | `1.0` | All |

## Source-Specific Fields

### gutMDisorder-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| association_id | integer | Unique microbe-disease association identifier | `12345` |
| microbe_rank | string | Taxonomic rank of microbe | `species`, `genus`, `family` |
| disease_id | string | Internal disease identifier | `DIS001` |

### MASI-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| interaction_id | integer | Unique signaling interaction identifier | `54321` |
| microbe_source | string | Microbe(s) producing the metabolite | `Faecalibacterium prausnitzii` |
| metabolite_class | string | Chemical class of metabolite | `Short-chain fatty acid`, `Indole`, `Bile acid derivative` |
| smiles | string | SMILES structure notation | `CCCC(=O)[O-]` |
| target_id | string | Internal target identifier | `TGT001` |
| symbol | string | Gene/protein symbol for target | `FFAR2` |
| tissue_expression | array | Tissues where target is expressed | `["colon", "adipose"]` |
| pathway_id | string | KEGG or Reactome pathway ID | `hsa04917` |
| pathway_category | string | Pathway classification category | `Signal transduction` |
| downstream_effects | array | Biological outcomes of interaction | `["anti-inflammatory"]` |
| kd | string | Dissociation constant for binding | `1.2 uM` |
| mechanism | string | Molecular mechanism description | `Competitive agonist` |

### VMH-Specific

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| organism_name | string | Full organism name for model | `Bacteroides thetaiotaomicron` |
| strain | string | Strain designation | `VPI-5482` |
| gram_stain | string | Gram staining classification | `negative` |
| oxygen_status | string | Oxygen requirement status | `anaerobic` |
| charge | integer | Molecular charge of metabolite | `-1` |
| inchi | string | InChI structure representation | `InChI=1S/...` |
| ec_number | string | Enzyme Commission classification | `2.7.1.1` |
| entrez_id | integer | NCBI Entrez Gene ID | `12345` |
| bigg_id | string | BiGG Models database identifier | `but` |

## Cross-References

| Database | ID Field | Usage |
|----------|----------|-------|
| NCBI Taxonomy | microbe_taxid | Microbe identification |
| ICD-10 | icd10_code | Disease classification |
| MeSH | mesh_id | Medical subject headings |
| Disease Ontology | doid | Disease ontology |
| PubMed | pmid | Literature evidence |
| PubChem | pubchem_cid | Compound identification |
| HMDB | hmdb_id | Metabolite data |
| ChEBI | chebi_id | Chemical ontology |
| UniProt | uniprot_id | Protein targets |
| KEGG | kegg_id | Pathway/compound mapping |
| Reactome | pathway_id | Pathway details |
| BiGG Models | bigg_id | Model compatibility |

## Field Mappings by Source

### gutMDisorder Mappings

| Source Field | Unified Field |
|--------------|---------------|
| association_id | interaction_id |
| microbe_taxid | microbe_taxid |
| microbe_name | microbe_name |
| microbe_rank | microbe_rank |
| disease_name | disease_association.disease_name |
| disease_category | disease_association.disease_category |
| icd10_code | disease_association.icd10_code |
| mesh_id | disease_association.mesh_id |
| doid | disease_association.doid |
| direction | disease_association.direction |
| effect_size | statistical_evidence.effect_size |
| p_value | statistical_evidence.p_value |
| sample_size | statistical_evidence.sample_size |
| cases | statistical_evidence.cases |
| controls | statistical_evidence.controls |
| study_design | study_metadata.study_design |
| method | study_metadata.method |
| pmid | pmid |

### MASI Mappings

| Source Field | Unified Field |
|--------------|---------------|
| interaction_id | interaction_id |
| microbe_name | microbe_name |
| microbe_source | metabolite_signaling.microbe_source |
| metabolite_id | metabolite_signaling.metabolite_id |
| metabolite_name | metabolite_signaling.metabolite_name |
| metabolite_class | metabolite_signaling.metabolite_class |
| pubchem_cid | metabolite_signaling.pubchem_cid |
| hmdb_id | metabolite_signaling.hmdb_id |
| chebi_id | metabolite_signaling.chebi_id |
| smiles | metabolite_signaling.smiles |
| target_name | host_target.target_name |
| target_type | host_target.target_type |
| symbol | host_target.symbol |
| uniprot_id | host_target.uniprot_id |
| effect | host_target.effect |
| tissue_expression | host_target.tissue_expression |
| pathway | pathway_effects.pathway |
| pathway_id | pathway_effects.pathway_id |
| pathway_category | pathway_effects.pathway_category |
| downstream_effects | pathway_effects.downstream_effects |
| kd | binding_affinity.kd |
| mechanism | binding_affinity.mechanism |
| pmid | pmid |

### VMH Mappings

| Source Field | Unified Field |
|--------------|---------------|
| microbe_taxid | microbe_taxid |
| microbe_name | microbe_name |
| model_id | metabolic_model.model_id |
| organism_name | metabolic_model.organism_name |
| strain | metabolic_model.strain |
| reactions_count | metabolic_model.reactions_count |
| metabolites_count | metabolic_model.metabolites_count |
| genes_count | metabolic_model.genes_count |
| metabolite_id | metabolite_info.metabolite_id |
| metabolite_name | metabolite_info.metabolite_name |
| formula | metabolite_info.formula |
| charge | metabolite_info.charge |
| compartment | metabolite_info.compartment |
| kegg_id | metabolite_info.kegg_id |
| hmdb_id | metabolite_info.hmdb_id |
| chebi_id | metabolite_info.chebi_id |
| pubchem_cid | metabolite_info.pubchem_cid |
| inchi | metabolite_info.inchi |
| rxn_id | reaction_info.rxn_id |
| equation | reaction_info.equation |
| subsystem | reaction_info.subsystem |
| gpr | reaction_info.gpr |
| lb | reaction_info.lb |
| ub | reaction_info.ub |
| reversibility | reaction_info.reversibility |
| ec_number | reaction_info.ec_number |
| gram_stain | microbe_phenotype.gram_stain |
| oxygen_status | microbe_phenotype.oxygen_status |
| uniprot_id | cross_references.uniprot_id |
| entrez_id | cross_references.entrez_id |
| bigg_id | cross_references.bigg_id |
| pathway | pathway_effects.pathway |

## VMH Compartment Codes

| Code | Compartment |
|------|-------------|
| c | Cytoplasm |
| e | Extracellular |
| m | Mitochondria |
| n | Nucleus |
| r | Endoplasmic reticulum |
| g | Golgi apparatus |
| l | Lysosome |
| x | Peroxisome |

## Disease Categories (gutMDisorder)

| Category | Examples |
|----------|----------|
| Metabolic disorders | Type 2 Diabetes, Obesity, NAFLD |
| Gastrointestinal | Crohn's Disease, Ulcerative Colitis, IBS |
| Neurological | Parkinson's Disease, Autism, Depression |
| Cardiovascular | Atherosclerosis, Hypertension |
| Autoimmune | Rheumatoid Arthritis, Multiple Sclerosis |
| Cancer | Colorectal Cancer, Liver Cancer |

## Metabolite Classes (MASI)

| Class | Examples |
|-------|----------|
| Short-chain fatty acids | Butyrate, Propionate, Acetate |
| Indoles | Indole, Indole-3-aldehyde, Indole-3-propionic acid |
| Bile acid derivatives | Deoxycholic acid, Lithocholic acid |
| Amino acid metabolites | p-Cresol, Phenylacetate |
| Vitamins | Folate, Vitamin K2, Biotin |
| Polyamines | Putrescine, Spermidine |
