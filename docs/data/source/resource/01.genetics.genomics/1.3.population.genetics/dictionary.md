# 1.3 Population Genetics - Data Dictionary

## Overview

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 1.3 |
| Subcategory Name | Population Genetics |
| Data Sources | 1000 Genomes, gnomAD, TOPMed, UK Biobank |
| Schema ID | https://gene.ai/schemas/1.3-population-genetics.json |

---

## Unified Fields

Core fields that are mapped and normalized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| chromosome | string | Required | Chromosome identifier | `1`, `22`, `X` | 1000 Genomes, gnomAD, TOPMed, UK Biobank |
| position | integer | Required | 1-based genomic coordinate | `10177`, `5227002` | 1000 Genomes, gnomAD, TOPMed, UK Biobank |
| reference_allele | string | Required | Reference genome nucleotide(s) | `A`, `T` | 1000 Genomes, gnomAD, TOPMed, UK Biobank |
| alternate_allele | string | Required | Variant nucleotide(s) | `AC`, `G` | 1000 Genomes, gnomAD, TOPMed, UK Biobank |
| allele_frequency | float (0-1) | Optional | Global allele frequency (proportion of alternate allele) | `0.425`, `0.0108` | 1000 Genomes (AF), gnomAD (AF), TOPMed (allele_freq) |
| allele_count | integer | Optional | Count of alternate alleles observed | `2130`, `153036` | 1000 Genomes (AC), gnomAD (AC), TOPMed (allele_count) |
| allele_number | integer | Optional | Total alleles in called genotypes | `5008`, `360000` | 1000 Genomes (AN), gnomAD (AN), TOPMed (allele_num) |
| filter_status | string | Optional | Quality filter status (PASS or filter reason) | `PASS`, `AC0`, `RF` | gnomAD, TOPMed |

---

## Population-Specific Allele Frequencies

### Continental Population Frequencies

| Field Name | Data Type | Cardinality | Description | Source Mappings |
|------------|-----------|-------------|-------------|-----------------|
| af_african | float | Optional | African/African-American population allele frequency | 1000 Genomes (AFR_AF), gnomAD (AC_afr derived) |
| af_american | float | Optional | American/Latino population allele frequency | 1000 Genomes (AMR_AF), gnomAD (AC_amr derived) |
| af_east_asian | float | Optional | East Asian population allele frequency | 1000 Genomes (EAS_AF), gnomAD (AC_eas derived) |
| af_european | float | Optional | European population allele frequency | 1000 Genomes (EUR_AF), gnomAD (AC_nfe derived) |
| af_south_asian | float | Optional | South Asian population allele frequency | 1000 Genomes (SAS_AF), gnomAD (AC_sas derived) |

---

## gnomAD-Specific Fields

### Genotype Counts

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gnomad_homozygote_count | integer | Optional | Number of homozygotes | nhomalt |
| gnomad_ac_xx | integer | Optional | Allele count in XX individuals | AC_XX |
| gnomad_ac_xy | integer | Optional | Allele count in XY individuals | AC_XY |

### Population-Specific Frequencies

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gnomad_af_ashkenazi | float | Optional | Ashkenazi Jewish allele frequency | AC_asj (derived) |
| gnomad_af_finnish | float | Optional | Finnish allele frequency | AC_fin (derived) |

### Filtering Allele Frequency

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gnomad_faf95 | float | Optional | Filtering allele frequency (95% CI) | faf95 |
| gnomad_faf99 | float | Optional | Filtering allele frequency (99% CI) | faf99 |
| gnomad_fafmax | float | Optional | Maximum filtering allele frequency across populations | fafmax |

### Quality Metrics

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| qc_fisher_strand | float | Optional | Fisher strand bias (FS) | FS |
| qc_mapping_quality | float | Optional | Mapping quality (MQ) | MQ |
| qc_quality_by_depth | float | Optional | Quality by depth (QD) | QD |
| qc_vqslod | float | Optional | VQSR log-odds score | VQSLOD |

### Gene Constraint Metrics

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| constraint_expected_lof | float | Optional | Expected loss-of-function variants | exp_lof |
| constraint_observed_lof | integer | Optional | Observed loss-of-function variants | obs_lof |
| constraint_oe_lof | float | Optional | Observed/expected ratio for LoF | oe_lof |
| constraint_pli | float | Optional | Probability of LoF intolerance (0-1) | pLI |
| constraint_loeuf | float | Optional | LoF observed/expected upper bound fraction | LOEUF |

### LOFTEE Annotations

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| loftee_annotation | string | Optional | LOFTEE consequence annotation | lof |
| loftee_filter | string | Optional | LOFTEE filter status | lof_filter |
| loftee_flags | string | Optional | LOFTEE warning flags | lof_flags |

### Functional Annotations

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| annotation_cadd_phred | float | Optional | CADD PHRED-scaled score | cadd_phred |
| annotation_revel | float | Optional | REVEL ensemble score | revel |
| annotation_sift | string | Optional | SIFT prediction | sift |
| annotation_polyphen | string | Optional | PolyPhen prediction | polyphen |
| annotation_spliceai_max | float | Optional | Maximum SpliceAI delta score | spliceai_ds_max |

---

## TOPMed-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| topmed_het_count | integer | Optional | Heterozygote count | het_count |
| topmed_hom_count | integer | Optional | Homozygote count | hom_count |
| topmed_sample_id | string | Optional | TOPMed sample identifier (authorized access only) | NWD_ID |

---

## UK Biobank-Specific Fields

### Participant Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| ukb_participant_id | integer | Optional | Encrypted participant ID | eid |
| ukb_field_id | integer | Optional | Data-field identifier | field_id |
| ukb_instance | integer | Optional | Visit instance (0=baseline) | instance |

### Genotyping Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| ukb_genotyping_array | string | Optional | Genotyping array type used | genotyping_array |
| ukb_batch | string | Optional | Genotyping batch | batch |

### Sample Quality Metrics

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| ukb_sex_inferred | integer | Optional | Genetically inferred sex (1=male, 2=female) | sex_inferred |
| ukb_heterozygosity | float | Optional | Sample heterozygosity | heterozygosity |
| ukb_missingness | float | Optional | Genotype call rate | missingness |
| ukb_genetic_ancestry | string | Optional | Genetic ancestry classification | genetic_ethnic_grouping |

---

## Field Mapping Reference

### 1000 Genomes Mappings

| Original Field | Unified Field |
|----------------|---------------|
| AF | allele_frequency |
| AC | allele_count |
| AN | allele_number |
| AFR_AF | af_african |
| AMR_AF | af_american |
| EAS_AF | af_east_asian |
| EUR_AF | af_european |
| SAS_AF | af_south_asian |

### gnomAD Mappings

| Original Field | Unified Field |
|----------------|---------------|
| AF | allele_frequency |
| AC | allele_count |
| AN | allele_number |
| nhomalt | gnomad_homozygote_count |
| AC_XX | gnomad_ac_xx |
| AC_XY | gnomad_ac_xy |
| faf95 | gnomad_faf95 |
| faf99 | gnomad_faf99 |
| fafmax | gnomad_fafmax |
| FS | qc_fisher_strand |
| MQ | qc_mapping_quality |
| QD | qc_quality_by_depth |
| VQSLOD | qc_vqslod |
| exp_lof | constraint_expected_lof |
| obs_lof | constraint_observed_lof |
| oe_lof | constraint_oe_lof |
| pLI | constraint_pli |
| LOEUF | constraint_loeuf |
| lof | loftee_annotation |
| lof_filter | loftee_filter |
| lof_flags | loftee_flags |
| cadd_phred | annotation_cadd_phred |
| revel | annotation_revel |
| sift | annotation_sift |
| polyphen | annotation_polyphen |
| spliceai_ds_max | annotation_spliceai_max |

### TOPMed Mappings

| Original Field | Unified Field |
|----------------|---------------|
| allele_freq | allele_frequency |
| allele_count | allele_count |
| allele_num | allele_number |
| het_count | topmed_het_count |
| hom_count | topmed_hom_count |
| NWD_ID | topmed_sample_id |

### UK Biobank Mappings

| Original Field | Unified Field |
|----------------|---------------|
| eid | ukb_participant_id |
| field_id | ukb_field_id |
| instance | ukb_instance |
| genotyping_array | ukb_genotyping_array |
| batch | ukb_batch |
| sex_inferred | ukb_sex_inferred |
| heterozygosity | ukb_heterozygosity |
| missingness | ukb_missingness |
| genetic_ethnic_grouping | ukb_genetic_ancestry |

---

## Interpretation Guide

### Allele Frequency Thresholds

| Frequency | Classification |
|-----------|----------------|
| >5% | Common variant |
| 1-5% | Low-frequency variant |
| 0.1-1% | Rare variant |
| <0.1% | Very rare variant |
| <0.01% | Ultra-rare variant |

### Gene Constraint (pLI) Interpretation

| pLI Score | Interpretation |
|-----------|----------------|
| >0.9 | Extremely LoF intolerant (haploinsufficient) |
| 0.5-0.9 | Moderately LoF intolerant |
| <0.5 | LoF tolerant |

### LOEUF Interpretation

| LOEUF Score | Interpretation |
|-------------|----------------|
| <0.35 | Strong constraint (LoF intolerant) |
| 0.35-1.0 | Moderate constraint |
| >1.0 | Little or no constraint |

---

*Generated: 2026-01-24*
