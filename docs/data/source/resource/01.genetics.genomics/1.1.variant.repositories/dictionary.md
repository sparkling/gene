# 1.1 Variant Repositories - Data Dictionary

## Overview

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 1.1 |
| Subcategory Name | Variant Repositories |
| Data Sources | ClinVar, dbSNP, dbVar |
| Schema ID | https://gene.ai/schemas/1.1-variant-repositories.json |

---

## Unified Fields

Core fields that are mapped and normalized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| variant_id | string | Required | Unique identifier for the variant record | `VCV000000123.4`, `rs334`, `nsv1234567` | ClinVar (VCV/RCV/SCV), dbSNP (rsid), dbVar (nsv/nssv) |
| chromosome | string | Required | Chromosome location of the variant | `11`, `chr1`, `chrX` | ClinVar, dbSNP, dbVar |
| position | integer | Required | 1-based genomic coordinate position on reference assembly | `5227002`, `10001`, `140753336` | ClinVar, dbSNP, dbVar |
| reference_allele | string | Required | Reference genome allele (IUPAC nucleotides) | `A`, `T`, `AG` | ClinVar, dbSNP, dbVar |
| alternate_allele | string | Required | Observed variant allele (IUPAC nucleotides) | `T`, `A`, `A` | ClinVar, dbSNP, dbVar |
| variant_type | string | Required | Sequence Ontology term for variant type classification | `SNV`, `deletion`, `insertion`, `structural_variant` | ClinVar, dbSNP, dbVar |
| clinical_significance | array[string] | Optional | ACMG/AMP clinical interpretation of variant pathogenicity | `Pathogenic`, `Likely pathogenic`, `Uncertain significance`, `Benign` | ClinVar, dbSNP (via ClinVar), dbVar |
| gene_symbols | array[string] | Optional | HGNC official gene symbol(s) associated with variant | `HBB`, `BRCA1`, `TP53` | ClinVar, dbSNP, dbVar |
| hgvs_expressions | array[string] | Optional | Human Genome Variation Society standardized notation | `NM_000518.5:c.20A>T`, `NC_000011.10:g.5227002T>A` | ClinVar, dbSNP |
| assembly | string | Required | Genome Reference Consortium build identifier | `GRCh38`, `GRCh37`, `hg38` | ClinVar, dbSNP, dbVar |

---

## ClinVar-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| clinvar_allele_id | integer | Optional | ClinVar unique simple variant identifier | AlleleID |
| clinvar_variation_id | integer | Optional | ClinVar variation record ID | VariationID |
| review_status | string | Optional | Aggregate review status level (stars) | ReviewStatus, CLNREVSTAT |
| number_submitters | integer | Optional | Number of submitting organizations | NumberSubmitters |
| date_last_evaluated | date | Optional | Most recent evaluation date | DateLastEvaluated |
| phenotype_ids | array[string] | Optional | MedGen CUI identifiers for associated phenotypes | PhenotypeIDS |
| somatic_clinical_impact | string | Optional | Somatic clinical impact classification | SomaticClinicalImpact |
| oncogenicity | string | Optional | Oncogenicity classification | Oncogenicity |
| origin | string | Optional | Variant origin (germline, somatic, de novo) | Origin |
| disease_names | array[string] | Optional | Associated disease name(s) | CLNDN |
| molecular_consequence | string | Optional | Molecular consequence annotation | MC |

---

## dbSNP-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| dbsnp_refsnp_id | string | Optional | RefSNP identifier (without rs prefix) | refsnp_id |
| dbsnp_create_date | datetime | Optional | Initial submission date to dbSNP | create_date |
| dbsnp_last_update | datetime | Optional | Most recent modification date | last_update_date |
| mane_select_ids | array[string] | Optional | MANE Select transcript IDs | mane_select_ids |
| allele_count | integer | Optional | Number of allele observations | allele_count |
| total_count | integer | Optional | Total chromosomes sampled | total_count |
| citations | array[integer] | Optional | PubMed IDs for associated publications | citations |
| population_frequencies | object | Optional | Population-specific allele frequencies (ALFA) | ALFA_frequencies |

---

## dbVar-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| dbvar_outermost_start | integer | Optional | 1-based start position (outermost boundary) | outermost_start |
| dbvar_outermost_stop | integer | Optional | 1-based stop position (outermost boundary) | outermost_stop |
| dbvar_inner_start | integer | Optional | Minimum affected start position | inner_start |
| dbvar_inner_stop | integer | Optional | Minimum affected stop position | inner_stop |
| sv_variant_count | integer | Optional | Number of supporting variants (structural variants) | variant_count |
| detection_method | string | Optional | Detection method(s) for structural variants | method |
| analysis_type | string | Optional | Analysis type(s) | analysis |
| platform | string | Optional | Sequencing/array platform | platform |
| study_accession | string | Optional | Source study accession(s) | study |
| sv_type | string | Optional | Structural variant type | SVTYPE |
| sv_length | integer | Optional | SV length (negative for deletions) | SVLEN |
| is_imprecise | boolean | Optional | Imprecise structural variant flag | IMPRECISE |
| size_category | string | Optional | Size category (small, medium, large) | bin_size |

---

## Field Mapping Reference

### ClinVar Mappings

| Original Field | Unified Field |
|----------------|---------------|
| AlleleID | clinvar_allele_id |
| VariationID | clinvar_variation_id |
| ReviewStatus | review_status |
| CLNREVSTAT | review_status |
| NumberSubmitters | number_submitters |
| DateLastEvaluated | date_last_evaluated |
| PhenotypeIDS | phenotype_ids |
| SomaticClinicalImpact | somatic_clinical_impact |
| Oncogenicity | oncogenicity |
| Origin | origin |
| CLNDN | disease_names |
| MC | molecular_consequence |
| GeneSymbol | gene_symbols |
| ClinicalSignificance | clinical_significance |
| CLNSIG | clinical_significance |

### dbSNP Mappings

| Original Field | Unified Field |
|----------------|---------------|
| refsnp_id | dbsnp_refsnp_id |
| create_date | dbsnp_create_date |
| last_update_date | dbsnp_last_update |
| mane_select_ids | mane_select_ids |
| allele_count | allele_count |
| total_count | total_count |
| citations | citations |
| ALFA_frequencies | population_frequencies |
| CHROM | chromosome |
| POS | position |
| REF | reference_allele |
| ALT | alternate_allele |

### dbVar Mappings

| Original Field | Unified Field |
|----------------|---------------|
| outermost_start | dbvar_outermost_start |
| outermost_stop | dbvar_outermost_stop |
| inner_start | dbvar_inner_start |
| inner_stop | dbvar_inner_stop |
| variant_count | sv_variant_count |
| method | detection_method |
| analysis | analysis_type |
| platform | platform |
| study | study_accession |
| SVTYPE | sv_type |
| SVLEN | sv_length |
| IMPRECISE | is_imprecise |
| bin_size | size_category |

---

## Semantic Definitions

| Field | Semantic Definition |
|-------|---------------------|
| variant_id | Primary key for variant identification across databases |
| chromosome | Chromosome identifier (1-22, X, Y, MT) |
| position | 1-based genomic coordinate on reference assembly |
| reference_allele | Reference genome allele (IUPAC nucleotides) |
| alternate_allele | Observed variant allele (IUPAC nucleotides) |
| variant_type | Sequence Ontology term for variant type |
| clinical_significance | ACMG/AMP classification of variant pathogenicity |
| gene_symbols | HGNC official gene symbol(s) |
| hgvs_expressions | Human Genome Variation Society standardized notation |
| assembly | Genome Reference Consortium build identifier |

---

*Generated: 2026-01-24*
