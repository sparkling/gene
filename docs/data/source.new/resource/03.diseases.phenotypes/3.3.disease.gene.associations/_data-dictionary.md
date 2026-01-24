# 3.3 Disease-Gene Associations & Pharmacogenomics - Data Dictionary

## Overview

This subcategory contains data linking genes to diseases with evidence scores, as well as pharmacogenomic information about drug-gene relationships and clinical recommendations.

**Data Sources:** DisGeNET, Open Targets, PharmGKB

---

## Unified Fields

These fields are harmonized across multiple data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `gene_id` | string | Optional (1:1) | Gene identifier (NCBI Gene ID or Ensembl ID) | DisGeNET, Open Targets, PharmGKB | `3157`, `ENSG00000146648`, `PA126` |
| `gene_symbol` | string | Required (1:1) | HGNC approved gene symbol | DisGeNET, Open Targets, PharmGKB | `EGFR`, `BRCA1`, `CYP2D6` |
| `disease_id` | string | Optional (1:1) | Disease identifier (UMLS CUI, EFO ID, etc.) | DisGeNET, Open Targets | `C0020443`, `EFO_0000311` |
| `disease_name` | string | Required (1:1) | Disease or phenotype name | DisGeNET, Open Targets, PharmGKB | `Hypercholesterolemia`, `Breast cancer` |
| `score` | float | Optional (1:1) | Association confidence score (0-1) | DisGeNET, Open Targets | `0.85`, `0.42`, `0.97` |
| `evidence_level` | string | Optional (1:1) | Evidence strength classification | DisGeNET, PharmGKB | `1A`, `2B`, `high confidence` |
| `pmids` | array[integer] | Optional (1:N) | PubMed IDs supporting the association | DisGeNET, Open Targets, PharmGKB | `[12345678, 23456789]` |

---

## Source-Specific Fields

### DisGeNET

| Field Name | Data Type | Cardinality | Description | Range/Pattern | Example Values |
|------------|-----------|-------------|-------------|---------------|----------------|
| `ei` | float | Optional | Evidence Index measuring consistency across publications | 0-1 | `0.95`, `0.72` |
| `dsi` | float | Optional | Disease Specificity Index (inverse of associated diseases) | 0-1 | `0.88`, `0.45` |
| `dpi` | float | Optional | Disease Pleiotropy Index (spread across disease classes) | 0-1 | `0.23`, `0.67` |
| `association_type` | enum | Optional | Type of gene-disease relationship | `Therapeutic`, `Biomarker`, `GeneticVariation`, `Causal`, `SomaticMutation`, `GermlineMutation`, `Susceptibility` | `GeneticVariation` |
| `snp_id` | string | Optional | dbSNP variant identifier for VDAs | `rs[0-9]+` | `rs12345678` |
| `chromosome` | string | Optional | Chromosome location | - | `chr17`, `chrX` |
| `position` | integer | Optional | GRCh38 genomic position | - | `43044295` |
| `consequence` | string | Optional | VEP variant consequence type | - | `missense_variant`, `stop_gained` |
| `number_of_pmids` | integer | Optional | Number of supporting publications | - | `15`, `142` |

---

### Open Targets

| Field Name | Data Type | Cardinality | Description | Pattern/Range | Example Values |
|------------|-----------|-------------|-------------|---------------|----------------|
| `target_id` | string | Optional | Ensembl gene ID for target | `ENSG[0-9]{11}` | `ENSG00000146648` |
| `ot_disease_id` | string | Optional | EFO disease identifier | `EFO_[0-9]+` | `EFO_0000311` |
| `datatype_scores` | object | Optional | Scores per data type (genetic, drug, pathway, etc.) | - | `{"genetic": 0.85, "drug": 0.72}` |
| `evidence_count` | integer | Optional | Total number of evidence strings | - | `245`, `38` |
| `tractability` | object | Optional | Druggability assessment by modality | - | `{"smallMolecule": true, "antibody": false}` |
| `safety` | object | Optional | Safety liability information | - | `{"hasSafetyLiability": true}` |
| `mechanism_of_action` | string | Optional | Drug mechanism of action | - | `EGFR inhibitor`, `HMG-CoA reductase inhibitor` |
| `clinical_phase` | integer | Optional | Clinical trial phase for known drugs | 0-4 | `3`, `4`, `1` |

---

### PharmGKB

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `pgkb_id` | string | Optional | PharmGKB accession for entities | `PA[0-9]+` | `PA126`, `PA447999` |
| `haplotype` | string | Optional | Star allele haplotype | - | `CYP2C9*2`, `CYP2D6*4` |
| `diplotype` | string | Optional | Diplotype | - | `CYP2C9*1/*2`, `CYP2D6*1/*4` |
| `phenotype` | string | Optional | Metabolizer phenotype classification | - | `Poor Metabolizer`, `Intermediate Metabolizer`, `Normal Metabolizer`, `Ultra-rapid Metabolizer` |
| `activity_value` | string | Optional | Functional impact of allele | - | `Decreased function`, `No function`, `Normal function` |
| `recommendation` | string | Optional | Clinical prescribing recommendation | - | `Consider alternative drug`, `Reduce dose by 50%` |
| `guideline_source` | enum | Optional | Source of clinical guideline | `CPIC`, `DPWG`, `FDA` | `CPIC` |
| `population_frequency` | object | Optional | Population allele frequencies by ethnicity | - | `{"European": 0.12, "East Asian": 0.05}` |

---

## Source Field Mappings

### DisGeNET Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `geneId` | `gene_id` |
| `geneSymbol` | `gene_symbol` |
| `diseaseId` | `disease_id` |
| `diseaseName` | `disease_name` |
| `score` | `score` |
| `EI` | `ei` |
| `DSI` | `dsi` |
| `DPI` | `dpi` |
| `associationType` | `association_type` |
| `snpId` | `snp_id` |
| `chromosome` | `chromosome` |
| `position` | `position` |
| `consequence` | `consequence` |
| `NofPmids` | `number_of_pmids` |
| `pmid` | `pmids` |

### Open Targets Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `targetId` | `target_id` |
| `approvedSymbol` | `gene_symbol` |
| `diseaseId` | `ot_disease_id` |
| `diseaseName` | `disease_name` |
| `score` | `score` |
| `datatypeScores` | `datatype_scores` |
| `evidenceCount` | `evidence_count` |
| `tractability` | `tractability` |
| `safety` | `safety` |
| `mechanismOfAction` | `mechanism_of_action` |
| `clinicalPhase` | `clinical_phase` |

### PharmGKB Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `pgkb_id` | `pgkb_id` |
| `Symbol` | `gene_symbol` |
| `Gene_ID` | `gene_id` |
| `Phenotype_Category` | `disease_name` |
| `Level_of_Evidence` | `evidence_level` |
| `Haplotype` | `haplotype` |
| `Diplotype` | `diplotype` |
| `Phenotype` | `phenotype` |
| `Activity_Value` | `activity_value` |
| `Recommendation` | `recommendation` |
| `Guideline_Source` | `guideline_source` |
| `Biogeographical_Frequency` | `population_frequency` |

---

## Evidence Level Classifications

### PharmGKB Evidence Levels

| Level | Description |
|-------|-------------|
| `1A` | Annotation supported by CPIC or PharmGKB clinical guideline |
| `1B` | Annotation supported by published clinical study |
| `2A` | Annotation supported by multiple studies showing association |
| `2B` | Annotation supported by single study |
| `3` | Annotation based on case reports or in vitro studies |
| `4` | Annotation based on conflicting evidence |

### DisGeNET Association Types

| Type | Description |
|------|-------------|
| `Therapeutic` | Gene is a drug target for the disease |
| `Biomarker` | Gene is a biomarker for the disease |
| `GeneticVariation` | Genetic variants in gene associated with disease |
| `Causal` | Gene is causally implicated in disease |
| `SomaticMutation` | Somatic mutations in gene contribute to disease |
| `GermlineMutation` | Germline mutations in gene cause disease |
| `Susceptibility` | Gene variants increase disease susceptibility |

---

## Metadata Fields

| Field Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `_source.database` | string | Name of the source database | `DisGeNET`, `Open_Targets`, `PharmGKB` |
| `_source.version` | string | Version of the source data | `v7.0`, `24.09`, `2024-01` |
| `_source.access_date` | date | Date the data was accessed | `2026-01-24` |
| `_source.original_id` | string | Original identifier in source | `C0020443` |
