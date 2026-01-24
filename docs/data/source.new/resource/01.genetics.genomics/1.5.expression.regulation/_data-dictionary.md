# 1.5 Expression & Regulation - Data Dictionary

## Overview

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 1.5 |
| Subcategory Name | Expression & Regulation |
| Data Sources | ENCODE, eQTLGen, GTEx, GWAS Catalog |
| Schema ID | https://gene.example.org/schemas/1.5-expression-regulation.json |

---

## Unified Fields

Core fields that are mapped and normalized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| gene_id | string | Required | Ensembl gene ID or HGNC symbol | `ENSG00000012048`, `BRCA1` | ENCODE, eQTLGen, GTEx, GWAS Catalog |
| variant_id | string | Optional | SNP or variant identifier | `rs12345`, `chr17_43044295_A_G_b38` | eQTLGen, GTEx, GWAS Catalog |
| p_value | float | Optional | Statistical significance of association | `1e-12`, `7e-8` | eQTLGen, GTEx, GWAS Catalog |
| effect_size | float | Optional | Effect size (beta/slope/OR) | `0.45`, `-0.23` | eQTLGen, GTEx, GWAS Catalog |

---

## ENCODE-Specific Fields

### Experiment Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| encode_accession | string | Optional | ENCODE accession (ENCSR, ENCFF, ENCBS) | accession |
| encode_assay_title | string | Optional | Assay name (ChIP-seq, ATAC-seq, RNA-seq, etc.) | assay_title |
| encode_target | object | Optional | Target protein/modification | target |
| biosample_ontology | object | Optional | Cell type/tissue information | biosample_ontology |

### File Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| file_format | string | Optional | File format: `bam`, `fastq`, `bigWig`, `bed`, `bigBed` | file_format |
| output_type | string | Optional | Output type: `peaks`, `signal`, `alignments` | output_type |
| assembly | string | Optional | Genome assembly: `GRCh38`, `hg19`, `mm10` | assembly |

### cCRE (Candidate Cis-Regulatory Elements)

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| ccre_accession | string | Optional | Candidate cis-regulatory element ID | cCRE_accession |
| ccre_classification | string | Optional | Classification: `PLS`, `pELS`, `dELS`, `DNase-H3K4me3`, `CTCF-only` | cCRE_classification |

### Genomic Coordinates

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| chromosome | string | Optional | Chromosome | chrom |
| start_position | integer | Optional | 0-based start position | chromStart |
| end_position | integer | Optional | 1-based end position | chromEnd |
| encode_score | integer | Optional | Normalized signal (0-1000) | score |

---

## eQTLGen-Specific Fields

### SNP Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| eqtlgen_snp | string | Optional | RS ID for eQTL SNP | SNP |
| eqtlgen_snp_chr | integer | Optional | SNP chromosome | SNPChr |
| eqtlgen_snp_pos | integer | Optional | SNP position (hg19) | SNPPos |

### Allele Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| assessed_allele | string | Optional | Effect allele | AssessedAllele |
| other_allele | string | Optional | Non-effect allele | OtherAllele |

### Association Statistics

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| z_score | float | Optional | Association Z-score | Zscore |
| fdr | float | Optional | False discovery rate | FDR |

### Gene Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gene_symbol | string | Optional | HGNC symbol | GeneSymbol |
| gene_chr | integer | Optional | Gene chromosome | GeneChr |
| gene_tss_position | integer | Optional | Gene TSS position | GenePos |

### Study Metadata

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| number_cohorts | integer | Optional | Number of contributing cohorts | NrCohorts |
| sample_size | integer | Optional | Total sample size | NrSamples |

---

## GTEx-Specific Fields

### Subject Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gtex_subject_id | string | Optional | Subject identifier (GTEX-XXXXX) | SUBJID |
| gtex_sex | integer | Optional | Sex: 1=Male, 2=Female | SEX |
| gtex_age | string | Optional | Age bracket (20-29, 30-39, etc.) | AGE |
| gtex_death_classification | integer | Optional | Hardy scale death classification (0-4) | DTHHRDY |

### Sample Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gtex_sample_id | string | Optional | Sample identifier | SAMPID |
| tissue_type_general | string | Optional | Tissue type (general) | SMTS |
| tissue_type_detailed | string | Optional | Tissue type (detailed) | SMTSD |
| rna_integrity_number | float | Optional | RNA Integrity Number (RIN) | SMRIN |
| analysis_freeze | string | Optional | Analysis freeze: `RNASEQ`, `WGS` | SMAFRZE |

### Expression Data

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| tpm | float | Optional | Transcripts Per Million | TPM |
| read_counts | integer | Optional | Raw read counts | read_counts |

### eQTL Statistics

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| tss_distance | integer | Optional | Distance to transcription start site | tss_distance |
| minor_allele_samples | integer | Optional | Number of samples with minor allele | ma_samples |
| minor_allele_count | integer | Optional | Minor allele count | ma_count |
| minor_allele_frequency | float | Optional | Minor allele frequency | maf |
| nominal_pvalue | float | Optional | Nominal p-value | pval_nominal |
| slope | float | Optional | Effect size (beta) | slope |
| slope_se | float | Optional | Standard error of slope | slope_se |
| beta_adjusted_pvalue | float | Optional | Beta-adjusted p-value | pval_beta |

### sQTL Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| sqtl_phenotype_id | string | Optional | sQTL intron cluster ID | phenotype_id |

---

## GWAS Catalog-Specific Fields

### Study Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gwas_study_accession | string | Optional | GCST study accession number | accessionId |
| disease_trait | string | Optional | Reported trait name | diseaseTrait |
| initial_sample_description | string | Optional | Sample description text | initialSampleSize |
| snp_count | integer | Optional | Number of SNPs analyzed | snpCount |
| is_imputed | boolean | Optional | Whether genotypes were imputed | imputed |

### Publication Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| pubmed_id | string | Optional | PubMed ID | pubmedId |
| publication_date | string | Optional | Publication date (ISO 8601) | publicationDate |
| author_info | object | Optional | First author information | author |

### Association Statistics

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| pvalue_description | string | Optional | P-value context (e.g., progression) | pvalueDescription |
| risk_allele_frequency | float | Optional | Risk allele frequency | riskFrequency |
| odds_ratio | float | Optional | Odds ratio per allele copy | orPerCopyNum |
| beta | float | Optional | Effect size for quantitative traits | betaNum |
| beta_unit | string | Optional | Unit of effect size | betaUnit |
| confidence_interval | string | Optional | Confidence interval notation | range |

### Variant Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| rsid | string | Optional | dbSNP RS identifier | rsId |
| functional_class | string | Optional | Sequence Ontology annotation | functionalClass |

### Trait Ontology

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| efo_trait | string | Optional | Experimental Factor Ontology trait | EFO_trait |
| efo_uri | string | Optional | EFO trait URI | EFO_uri |
| ancestral_groups | array[string] | Optional | Population ancestry groups | ancestralGroups |

---

## Field Mapping Reference

### ENCODE Mappings

| Original Field | Unified Field |
|----------------|---------------|
| accession | encode_accession |
| assay_title | encode_assay_title |
| target | encode_target |
| biosample_ontology | biosample_ontology |
| file_format | file_format |
| output_type | output_type |
| assembly | assembly |
| cCRE_accession | ccre_accession |
| cCRE_classification | ccre_classification |
| chrom | chromosome |
| chromStart | start_position |
| chromEnd | end_position |
| score | encode_score |

### eQTLGen Mappings

| Original Field | Unified Field |
|----------------|---------------|
| SNP | eqtlgen_snp |
| SNPChr | eqtlgen_snp_chr |
| SNPPos | eqtlgen_snp_pos |
| AssessedAllele | assessed_allele |
| OtherAllele | other_allele |
| Zscore | z_score |
| Gene | gene_id |
| GeneSymbol | gene_symbol |
| GeneChr | gene_chr |
| GenePos | gene_tss_position |
| NrCohorts | number_cohorts |
| NrSamples | sample_size |
| FDR | fdr |
| Pvalue | p_value |

### GTEx Mappings

| Original Field | Unified Field |
|----------------|---------------|
| SUBJID | gtex_subject_id |
| SEX | gtex_sex |
| AGE | gtex_age |
| DTHHRDY | gtex_death_classification |
| SAMPID | gtex_sample_id |
| SMTS | tissue_type_general |
| SMTSD | tissue_type_detailed |
| SMRIN | rna_integrity_number |
| SMAFRZE | analysis_freeze |
| TPM | tpm |
| read_counts | read_counts |
| tss_distance | tss_distance |
| ma_samples | minor_allele_samples |
| ma_count | minor_allele_count |
| maf | minor_allele_frequency |
| pval_nominal | nominal_pvalue |
| slope | slope |
| slope_se | slope_se |
| pval_beta | beta_adjusted_pvalue |
| phenotype_id | sqtl_phenotype_id |

### GWAS Catalog Mappings

| Original Field | Unified Field |
|----------------|---------------|
| accessionId | gwas_study_accession |
| diseaseTrait | disease_trait |
| initialSampleSize | initial_sample_description |
| snpCount | snp_count |
| imputed | is_imputed |
| pubmedId | pubmed_id |
| publicationDate | publication_date |
| author | author_info |
| pvalue | p_value |
| pvalueDescription | pvalue_description |
| riskFrequency | risk_allele_frequency |
| orPerCopyNum | odds_ratio |
| betaNum | beta |
| betaUnit | beta_unit |
| range | confidence_interval |
| rsId | rsid |
| chromosomeName | chromosome |
| chromosomePosition | position |
| functionalClass | functional_class |
| EFO_trait | efo_trait |
| EFO_uri | efo_uri |
| ancestralGroups | ancestral_groups |

---

## Interpretation Guide

### P-value Significance Thresholds

| Threshold | Description |
|-----------|-------------|
| 5x10^-8 | Genome-wide significance (GWAS) |
| 1x10^-5 | Suggestive significance |
| 0.05 | Nominal significance |

### cCRE Classifications

| Classification | Description |
|----------------|-------------|
| PLS | Promoter-like signature |
| pELS | Proximal enhancer-like signature |
| dELS | Distal enhancer-like signature |
| DNase-H3K4me3 | DNase with H3K4me3 |
| CTCF-only | CTCF binding only |

### Hardy Scale (Death Classification)

| Score | Description |
|-------|-------------|
| 0 | Ventilator case |
| 1 | Violent and fast death |
| 2 | Fast death of natural causes |
| 3 | Intermediate death |
| 4 | Slow death |

### Common GTEx Tissues

| Tissue Code | Description |
|-------------|-------------|
| Whole_Blood | Whole blood |
| Liver | Liver |
| Brain_Cortex | Brain cortex |
| Heart_Left_Ventricle | Heart left ventricle |
| Muscle_Skeletal | Skeletal muscle |
| Adipose_Subcutaneous | Subcutaneous adipose tissue |

---

*Generated: 2026-01-24*
