# 1.2 Functional Prediction - Data Dictionary

## Overview

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 1.2 |
| Subcategory Name | Functional Prediction |
| Data Sources | AlphaMissense, CADD, dbNSFP, SpliceAI |
| Schema ID | https://gene.example.org/schemas/1.2-functional-prediction.json |

---

## Unified Fields

Core fields that are mapped and normalized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| chromosome | string | Required | Chromosome identifier | `chr1`, `1`, `22` | AlphaMissense, CADD, dbNSFP, SpliceAI |
| position | integer | Required | 1-based genomic position | `10001`, `5227002` | AlphaMissense, CADD, dbNSFP, SpliceAI |
| reference_allele | string | Required | Reference genome nucleotide(s) | `A`, `T`, `G` | AlphaMissense, CADD, dbNSFP, SpliceAI |
| alternate_allele | string | Required | Variant nucleotide(s) | `G`, `A`, `C` | AlphaMissense, CADD, dbNSFP, SpliceAI |
| gene_symbols | array[string] | Optional | HGNC gene symbol(s) | `BRCA1`, `TP53` | AlphaMissense, dbNSFP, SpliceAI |
| transcript_ids | array[string] | Optional | Ensembl or RefSeq transcript ID(s) | `ENST00000357654`, `NM_007294.4` | AlphaMissense, dbNSFP |

---

## AlphaMissense-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| alphamissense_genome | string | Optional | Reference genome version (hg38) | genome |
| alphamissense_uniprot_id | string | Optional | UniProt accession | uniprot_id |
| alphamissense_protein_variant | string | Optional | Amino acid change (e.g., E7V) | protein_variant |
| alphamissense_pathogenicity | float (0-1) | Optional | AlphaMissense pathogenicity score | am_pathogenicity |
| alphamissense_class | string | Optional | Classification: `likely_benign`, `ambiguous`, `likely_pathogenic` | am_class |
| alphamissense_gene_mean | float | Optional | Gene-level mean pathogenicity | mean_am_pathogenicity |

---

## CADD-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| cadd_raw_score | float | Optional | CADD untransformed SVM score | RawScore |
| cadd_phred | float | Optional | CADD PHRED-scaled score (-10*log10(rank/total)) | PHRED |

---

## dbNSFP-Specific Fields

### Amino Acid Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| dbnsfp_ref_amino_acid | string | Optional | Reference amino acid | aaref |
| dbnsfp_alt_amino_acid | string | Optional | Alternate amino acid | aaalt |

### Gene/Transcript Identifiers

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| ensembl_gene_id | string | Optional | Ensembl gene ID | Ensembl_geneid |
| ensembl_transcript_id | string | Optional | Ensembl transcript ID | Ensembl_transcriptid |
| ensembl_protein_id | string | Optional | Ensembl protein ID | Ensembl_proteinid |
| uniprot_accession | string | Optional | UniProt accession | Uniprot_acc |

### HGVS Annotations

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| hgvs_coding | string | Optional | HGVS coding annotation | HGVSc_snpEff |
| hgvs_protein | string | Optional | HGVS protein annotation | HGVSp_snpEff |

### Prediction Scores

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| sift_score | float | Optional | SIFT prediction score (lower = more damaging) | SIFT_score |
| sift_prediction | string | Optional | SIFT prediction (D=Deleterious/T=Tolerated) | SIFT_pred |
| polyphen2_hdiv_score | float | Optional | PolyPhen-2 HumDiv score | Polyphen2_HDIV_score |
| polyphen2_hdiv_prediction | string | Optional | PolyPhen-2 HumDiv prediction | Polyphen2_HDIV_pred |
| polyphen2_hvar_score | float | Optional | PolyPhen-2 HumVar score | Polyphen2_HVAR_score |
| polyphen2_hvar_prediction | string | Optional | PolyPhen-2 HumVar prediction | Polyphen2_HVAR_pred |
| lrt_score | float | Optional | LRT (Likelihood Ratio Test) score | LRT_score |
| mutationtaster_score | float | Optional | MutationTaster score | MutationTaster_score |
| mutationassessor_score | float | Optional | MutationAssessor score | MutationAssessor_score |
| fathmm_score | float | Optional | FATHMM score | FATHMM_score |
| provean_score | float | Optional | PROVEAN score | PROVEAN_score |

### Ensemble Scores

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| metasvm_score | float | Optional | MetaSVM ensemble score | MetaSVM_score |
| metalr_score | float | Optional | MetaLR ensemble score | MetaLR_score |
| revel_score | float | Optional | REVEL ensemble score (0-1, higher = more pathogenic) | REVEL_score |

### Conservation Scores

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| phylop_100way | float | Optional | phyloP conservation score (100 vertebrates) | phyloP100way_vertebrate |
| phastcons_100way | float | Optional | phastCons conservation score | phastCons100way_vertebrate |
| gerp_rs | float | Optional | GERP++ rejected substitution score | GERP++_RS |

### Population Frequencies

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| gnomad_exomes_af | float | Optional | gnomAD exome allele frequency | gnomAD_exomes_AF |
| gnomad_genomes_af | float | Optional | gnomAD genome allele frequency | gnomAD_genomes_AF |
| thousand_genomes_af | float | Optional | 1000 Genomes Phase 3 allele frequency | 1000Gp3_AF |
| topmed_af | float | Optional | TOPMed allele frequency | TOPMed_AF |

### ClinVar Integration

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| clinvar_significance | string | Optional | ClinVar clinical significance | clinvar_clnsig |
| clinvar_review | string | Optional | ClinVar review status | clinvar_review |

---

## SpliceAI-Specific Fields

### Delta Scores (Splice Effect Predictions)

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| spliceai_allele | string | Optional | SpliceAI alternate allele | ALLELE |
| spliceai_gene_symbol | string | Optional | SpliceAI gene symbol | SYMBOL |
| spliceai_ds_ag | float | Optional | Delta score for acceptor gain (0-1) | DS_AG |
| spliceai_ds_al | float | Optional | Delta score for acceptor loss (0-1) | DS_AL |
| spliceai_ds_dg | float | Optional | Delta score for donor gain (0-1) | DS_DG |
| spliceai_ds_dl | float | Optional | Delta score for donor loss (0-1) | DS_DL |

### Distance to Splice Site

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| spliceai_dp_ag | integer | Optional | Distance to acceptor gain position | DP_AG |
| spliceai_dp_al | integer | Optional | Distance to acceptor loss position | DP_AL |
| spliceai_dp_dg | integer | Optional | Distance to donor gain position | DP_DG |
| spliceai_dp_dl | integer | Optional | Distance to donor loss position | DP_DL |

---

## Field Mapping Reference

### AlphaMissense Mappings

| Original Field | Unified Field |
|----------------|---------------|
| genome | alphamissense_genome |
| uniprot_id | alphamissense_uniprot_id |
| protein_variant | alphamissense_protein_variant |
| am_pathogenicity | alphamissense_pathogenicity |
| am_class | alphamissense_class |
| mean_am_pathogenicity | alphamissense_gene_mean |

### CADD Mappings

| Original Field | Unified Field |
|----------------|---------------|
| RawScore | cadd_raw_score |
| PHRED | cadd_phred |

### dbNSFP Mappings

| Original Field | Unified Field |
|----------------|---------------|
| aaref | dbnsfp_ref_amino_acid |
| aaalt | dbnsfp_alt_amino_acid |
| genename | gene_symbols |
| Ensembl_geneid | ensembl_gene_id |
| Ensembl_transcriptid | ensembl_transcript_id |
| Ensembl_proteinid | ensembl_protein_id |
| Uniprot_acc | uniprot_accession |
| HGVSc_snpEff | hgvs_coding |
| HGVSp_snpEff | hgvs_protein |
| SIFT_score | sift_score |
| SIFT_pred | sift_prediction |
| Polyphen2_HDIV_score | polyphen2_hdiv_score |
| Polyphen2_HDIV_pred | polyphen2_hdiv_prediction |
| Polyphen2_HVAR_score | polyphen2_hvar_score |
| Polyphen2_HVAR_pred | polyphen2_hvar_prediction |
| LRT_score | lrt_score |
| MutationTaster_score | mutationtaster_score |
| MutationAssessor_score | mutationassessor_score |
| FATHMM_score | fathmm_score |
| PROVEAN_score | provean_score |
| MetaSVM_score | metasvm_score |
| MetaLR_score | metalr_score |
| REVEL_score | revel_score |
| CADD_phred | cadd_phred |
| AlphaMissense_score | alphamissense_pathogenicity |
| phyloP100way_vertebrate | phylop_100way |
| phastCons100way_vertebrate | phastcons_100way |
| GERP++_RS | gerp_rs |
| gnomAD_exomes_AF | gnomad_exomes_af |
| gnomAD_genomes_AF | gnomad_genomes_af |
| 1000Gp3_AF | thousand_genomes_af |
| TOPMed_AF | topmed_af |
| clinvar_clnsig | clinvar_significance |
| clinvar_review | clinvar_review |

### SpliceAI Mappings

| Original Field | Unified Field |
|----------------|---------------|
| ALLELE | spliceai_allele |
| SYMBOL | spliceai_gene_symbol |
| DS_AG | spliceai_ds_ag |
| DS_AL | spliceai_ds_al |
| DS_DG | spliceai_ds_dg |
| DS_DL | spliceai_ds_dl |
| DP_AG | spliceai_dp_ag |
| DP_AL | spliceai_dp_al |
| DP_DG | spliceai_dp_dg |
| DP_DL | spliceai_dp_dl |

---

## Score Interpretation Guide

| Score | Range | Interpretation |
|-------|-------|----------------|
| CADD PHRED | 0-50+ | >20 = top 1% most deleterious, >30 = top 0.1% |
| AlphaMissense | 0-1 | >0.564 = likely pathogenic, <0.340 = likely benign |
| REVEL | 0-1 | >0.5 = likely pathogenic |
| SIFT | 0-1 | <0.05 = deleterious |
| PolyPhen-2 | 0-1 | >0.85 = probably damaging |
| SpliceAI | 0-1 | >0.5 = high splice impact, >0.2 = moderate |

---

*Generated: 2026-01-24*
