# 1.6 Cancer Genomics - Data Dictionary

## Overview

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 1.6 |
| Subcategory Name | Cancer Genomics |
| Data Sources | BRCA Exchange, Cancer Gene Census, cBioPortal, CIViC, COSMIC, OncoKB |
| Schema ID | https://gene.example.org/schemas/1.6-cancer-genomics.json |

---

## Unified Fields

Core fields that are mapped and normalized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| gene | string | Required | HGNC gene symbol | `BRAF`, `TP53`, `BRCA1` | BRCA Exchange, Cancer Gene Census, cBioPortal, CIViC, COSMIC, OncoKB |
| variant | string | Optional | Protein-level variant notation | `V600E`, `p.Glu23ValfsTer17`, `c.1799T>A` | BRCA Exchange, cBioPortal, CIViC, COSMIC, OncoKB |
| clinical_significance | string | Optional | Pathogenicity or oncogenicity classification | `Pathogenic`, `Oncogenic`, `Likely Oncogenic` | BRCA Exchange, CIViC, OncoKB |
| cancer_types | array[string] | Optional | Associated cancer type(s) | `Melanoma`, `Colorectal`, `Breast` | Cancer Gene Census, cBioPortal, CIViC, COSMIC, OncoKB |
| evidence_level | string | Optional | Tier classification for clinical evidence | `A`, `B`, `1`, `2`, `3A` | CIViC, OncoKB |

---

## BRCA Exchange-Specific Fields

### Variant Identification

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| brca_genomic_coordinate | string | Optional | chr:pos:ref:alt format coordinate (hg38) | Genomic_Coordinate_hg38 |
| brca_hgvs_cdna | string | Optional | Coding DNA change | HGVS_cDNA |
| brca_hgvs_protein | string | Optional | Protein change | HGVS_Protein |
| brca_reference_sequence | string | Optional | Transcript reference | Reference_Sequence |

### Pathogenicity Classifications

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| brca_pathogenicity_expert | string | Optional | ENIGMA expert classification | Pathogenicity_expert |
| brca_pathogenicity_aggregated | string | Optional | Aggregated classification from all sources | Pathogenicity_all |
| brca_clinvar_significance | string | Optional | ClinVar classification | Clinical_significance_ClinVar |

### Metadata

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| brca_date_evaluated | date | Optional | ENIGMA review date | Date_last_evaluated_ENIGMA |
| brca_source_databases | string | Optional | Contributing databases | Source |
| brca_clinvar_scv | string | Optional | ClinVar accession | SCV_ClinVar |

---

## Cancer Gene Census-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| cgc_gene_name | string | Optional | Full gene name | Name |
| entrez_gene_id | integer | Optional | NCBI Entrez gene ID | Entrez GeneId |
| cgc_genome_location | string | Optional | Chromosomal location | Genome Location |
| cgc_tier | integer | Optional | Tier: 1 (strong evidence) or 2 (emerging) | Tier |
| cgc_hallmark | string | Optional | Cancer hallmarks | Hallmark |
| role_in_cancer | string | Optional | Role: `oncogene`, `TSG`, `fusion` | Role in Cancer |
| mutation_types | string | Optional | Mutation type codes (Mis, N, F, S, D, T, A) | Mutation Types |
| translocation_partners | string | Optional | Fusion partner genes | Translocation Partner |
| cgc_synonyms | string | Optional | Alternative gene names | Synonyms |

---

## cBioPortal-Specific Fields

### Study Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| cbio_study_id | string | Optional | Unique study identifier | studyId |
| cbio_cancer_type_id | string | Optional | Cancer type classification | cancerTypeId |
| cbio_molecular_profile_id | string | Optional | Molecular data profile ID | molecularProfileId |
| cbio_alteration_type | string | Optional | Type: `MUTATION_EXTENDED`, `COPY_NUMBER_ALTERATION`, etc. | molecularAlterationType |
| cbio_datatype | string | Optional | Data type: `MAF`, `DISCRETE`, `LOG2-VALUE` | datatype |
| cbio_sample_count | integer | Optional | Total samples in study | allSampleCount |
| reference_genome | string | Optional | Reference genome: `hg19` or `hg38` | referenceGenome |

### Mutation Details

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| chromosome | string | Optional | Chromosome number | Chromosome |
| start_position | integer | Optional | Start genomic coordinate | Start_Position |
| end_position | integer | Optional | End genomic coordinate | End_Position |
| variant_classification | string | Optional | Mutation functional class | Variant_Classification |
| variant_type | string | Optional | Variant type: `SNP`, `INS`, `DEL` | Variant_Type |
| reference_allele | string | Optional | Reference nucleotide(s) | Reference_Allele |
| tumor_allele | string | Optional | Tumor allele | Tumor_Seq_Allele2 |
| sample_barcode | string | Optional | Sample identifier | Tumor_Sample_Barcode |

---

## CIViC-Specific Fields

### Entity Identification

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| civic_id | integer | Optional | CIViC entity ID | id |
| civic_variant_name | string | Optional | Variant name in CIViC | name |
| civic_variant_types | array[string] | Optional | Sequence Ontology terms | variant_types |
| civic_coordinates | object | Optional | Genomic coordinates | coordinates |
| civic_hgvs_expressions | array[string] | Optional | HGVS nomenclature | hgvs_expressions |

### Evidence Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| civic_evidence_type | string | Optional | Type: `Predictive`, `Diagnostic`, `Prognostic`, `Predisposing` | evidence_type |
| civic_evidence_direction | string | Optional | Direction: `Supports`, `Does Not Support` | evidence_direction |
| civic_clinical_significance | string | Optional | Significance: `Sensitivity`, `Resistance`, etc. | clinical_significance |

### Clinical Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| civic_disease | object | Optional | Cancer type with DOID | disease |
| associated_drugs | array[string] | Optional | Associated therapies | drugs |
| civic_source | object | Optional | Publication reference | source |
| amp_level | string | Optional | AMP/ASCO/CAP tier classification | amp_level |

---

## COSMIC-Specific Fields

### Mutation Identification

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| cosmic_cosv_id | string | Optional | Genomic mutation ID (COSV + integer) | COSV_ID |
| cosmic_cosm_id | string | Optional | Legacy mutation ID (COSM + integer) | COSM_ID |
| cosmic_transcript | string | Optional | Ensembl transcript ID | Transcript |
| cosmic_cds_mutation | string | Optional | CDS-level HGVS | CDS_Mutation |
| cosmic_aa_mutation | string | Optional | Protein-level HGVS | AA_Mutation |
| cosmic_mutation_type | string | Optional | Type: `Substitution`, `Insertion`, `Deletion` | Mutation_Type |
| cosmic_mutation_description | string | Optional | Functional consequence | Mutation_Description |

### Sample Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| cosmic_sample_id | integer | Optional | Internal sample ID | Sample_ID |
| cosmic_sample_name | string | Optional | Sample identifier | Sample_Name |
| primary_site | string | Optional | Anatomical site | Primary_Site |
| primary_histology | string | Optional | Histological type | Primary_Histology |
| sample_type | string | Optional | Type: `cell line`, `tumour`, etc. | Sample_Type |
| patient_age | integer | Optional | Patient age | Age |

### Mutational Signatures

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| mutational_signature | string | Optional | Mutational signature ID (SBS/DBS/ID) | Signature |
| signature_type | string | Optional | Type: `SBS`, `DBS`, `ID` | Signature_Type |
| signature_contribution | float | Optional | Signature proportion (0-1) | Contribution |
| proposed_etiology | string | Optional | Mutagenic process | Proposed_Etiology |

### Gene Fusions

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| fusion_id | integer | Optional | Gene fusion record ID | Fusion_ID |
| fusion_gene_5prime | string | Optional | 5' fusion partner | Gene_5prime |
| fusion_gene_3prime | string | Optional | 3' fusion partner | Gene_3prime |

---

## OncoKB-Specific Fields

### Gene Classification

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| oncokb_is_oncogene | boolean | Optional | Is oncogene | oncogene |
| oncokb_is_tsg | boolean | Optional | Is tumor suppressor gene | tsg |
| oncokb_gene_aliases | array[string] | Optional | Alternative gene names | geneAliases |
| oncokb_background | string | Optional | Gene description | background |

### Alteration Details

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| oncokb_alteration | string | Optional | Alteration name (e.g., V600E) | alteration |
| oncokb_alteration_type | string | Optional | Type: `Mutation`, `CNA`, `Fusion` | alterationType |
| oncokb_consequence | string | Optional | Protein consequence | consequence |
| oncokb_protein_start | integer | Optional | Protein position start | proteinStart |
| oncokb_protein_end | integer | Optional | Protein position end | proteinEnd |

### Oncogenicity and Effect

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| oncokb_oncogenic | string | Optional | Classification: `Oncogenic`, `Likely Oncogenic`, etc. | oncogenic |
| oncokb_mutation_effect | object | Optional | Effect: `Gain of function`, `Loss of function` | mutationEffect |
| oncokb_level | string | Optional | Evidence level: `1`, `2`, `3A`, `3B`, `4`, `R1`, `R2` | level |

### Supporting Evidence

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| oncokb_pmids | array[integer] | Optional | Supporting PubMed IDs | pmids |
| oncokb_abstracts | array[string] | Optional | Conference abstracts | abstracts |
| oncokb_description | string | Optional | Clinical context description | description |

---

## Field Mapping Reference

### BRCA Exchange Mappings

| Original Field | Unified Field |
|----------------|---------------|
| Gene_Symbol | gene |
| Genomic_Coordinate_hg38 | brca_genomic_coordinate |
| HGVS_cDNA | brca_hgvs_cdna |
| HGVS_Protein | brca_hgvs_protein |
| Reference_Sequence | brca_reference_sequence |
| Pathogenicity_expert | brca_pathogenicity_expert |
| Pathogenicity_all | brca_pathogenicity_aggregated |
| Clinical_significance_ENIGMA | clinical_significance |
| Clinical_significance_ClinVar | brca_clinvar_significance |
| Date_last_evaluated_ENIGMA | brca_date_evaluated |
| Source | brca_source_databases |
| SCV_ClinVar | brca_clinvar_scv |

### Cancer Gene Census Mappings

| Original Field | Unified Field |
|----------------|---------------|
| Gene Symbol | gene |
| Name | cgc_gene_name |
| Entrez GeneId | entrez_gene_id |
| Genome Location | cgc_genome_location |
| Tier | cgc_tier |
| Hallmark | cgc_hallmark |
| Role in Cancer | role_in_cancer |
| Mutation Types | mutation_types |
| Translocation Partner | translocation_partners |
| Cancer Types | cancer_types |
| Synonyms | cgc_synonyms |

### cBioPortal Mappings

| Original Field | Unified Field |
|----------------|---------------|
| studyId | cbio_study_id |
| cancerTypeId | cbio_cancer_type_id |
| molecularProfileId | cbio_molecular_profile_id |
| molecularAlterationType | cbio_alteration_type |
| datatype | cbio_datatype |
| allSampleCount | cbio_sample_count |
| referenceGenome | reference_genome |
| Hugo_Symbol | gene |
| Entrez_Gene_Id | entrez_gene_id |
| Chromosome | chromosome |
| Start_Position | start_position |
| End_Position | end_position |
| Variant_Classification | variant_classification |
| Variant_Type | variant_type |
| Reference_Allele | reference_allele |
| Tumor_Seq_Allele2 | tumor_allele |
| Tumor_Sample_Barcode | sample_barcode |

### CIViC Mappings

| Original Field | Unified Field |
|----------------|---------------|
| id | civic_id |
| name | civic_variant_name |
| entrez_id | entrez_gene_id |
| variant_types | civic_variant_types |
| coordinates | civic_coordinates |
| hgvs_expressions | civic_hgvs_expressions |
| evidence_type | civic_evidence_type |
| evidence_level | evidence_level |
| evidence_direction | civic_evidence_direction |
| clinical_significance | civic_clinical_significance |
| disease | civic_disease |
| drugs | associated_drugs |
| source | civic_source |
| amp_level | amp_level |

### COSMIC Mappings

| Original Field | Unified Field |
|----------------|---------------|
| COSV_ID | cosmic_cosv_id |
| COSM_ID | cosmic_cosm_id |
| Transcript | cosmic_transcript |
| CDS_Mutation | cosmic_cds_mutation |
| AA_Mutation | cosmic_aa_mutation |
| Mutation_Type | cosmic_mutation_type |
| Mutation_Description | cosmic_mutation_description |
| GRCh38_Position | chromosome |
| Sample_ID | cosmic_sample_id |
| Sample_Name | cosmic_sample_name |
| Primary_Site | primary_site |
| Primary_Histology | primary_histology |
| Sample_Type | sample_type |
| Age | patient_age |
| census_tier | cgc_tier |
| Role_in_Cancer | role_in_cancer |
| Signature | mutational_signature |
| Signature_Type | signature_type |
| Contribution | signature_contribution |
| Proposed_Etiology | proposed_etiology |
| Fusion_ID | fusion_id |
| Gene_5prime | fusion_gene_5prime |
| Gene_3prime | fusion_gene_3prime |

### OncoKB Mappings

| Original Field | Unified Field |
|----------------|---------------|
| hugoSymbol | gene |
| oncogene | oncokb_is_oncogene |
| tsg | oncokb_is_tsg |
| geneAliases | oncokb_gene_aliases |
| background | oncokb_background |
| alteration | oncokb_alteration |
| alterationType | oncokb_alteration_type |
| consequence | oncokb_consequence |
| proteinStart | oncokb_protein_start |
| proteinEnd | oncokb_protein_end |
| oncogenic | oncokb_oncogenic |
| mutationEffect | oncokb_mutation_effect |
| level | oncokb_level |
| cancerTypes | cancer_types |
| drugs | associated_drugs |
| pmids | oncokb_pmids |
| abstracts | oncokb_abstracts |
| description | oncokb_description |

---

## Interpretation Guide

### OncoKB Evidence Levels

| Level | Description |
|-------|-------------|
| 1 | FDA-recognized biomarker predictive of response to FDA-approved drug |
| 2 | Standard care biomarker recommended by NCCN/other expert panels |
| 3A | Compelling clinical evidence supports use in this indication |
| 3B | Standard care or investigational biomarker in another indication |
| 4 | Compelling biological evidence supports the biomarker |
| R1 | Standard care biomarker predictive of resistance to FDA-approved drug |
| R2 | Compelling clinical evidence supports the biomarker as resistance |

### CIViC Evidence Levels

| Level | Description |
|-------|-------------|
| A | Validated association in clinical practice |
| B | Clinical evidence from trials |
| C | Case study |
| D | Preclinical evidence |
| E | Inferential association |

### Cancer Gene Census Tiers

| Tier | Description |
|------|-------------|
| 1 | Strong evidence of cancer causation |
| 2 | Emerging evidence, requires further investigation |

### Mutation Type Codes

| Code | Description |
|------|-------------|
| Mis | Missense |
| N | Nonsense |
| F | Frameshift |
| S | Splice site |
| D | Deletion |
| T | Translocation |
| A | Amplification |

### Common Mutational Signatures

| Signature | Proposed Etiology |
|-----------|------------------|
| SBS1 | Age-related (spontaneous deamination) |
| SBS2/13 | APOBEC activity |
| SBS4 | Tobacco smoking |
| SBS6 | MMR deficiency |
| SBS7 | UV light exposure |
| SBS22 | Aristolochic acid exposure |

---

*Generated: 2026-01-24*
