# 1.4 Pharmacogenomics - Data Dictionary

## Overview

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 1.4 |
| Subcategory Name | Pharmacogenomics |
| Data Sources | CPIC, DPWG, PharmGKB, PharmVar |
| Schema ID | https://gene.ai/schemas/1.4-pharmacogenomics.json |

---

## Unified Fields

Core fields that are mapped and normalized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| gene | string | Required | HGNC gene symbol for pharmacogene | `CYP2D6`, `CYP2C19`, `SLCO1B1` | CPIC, DPWG, PharmGKB, PharmVar |
| drugs | array[string] | Optional | Generic drug name(s) | `codeine`, `clopidogrel`, `simvastatin` | CPIC, DPWG, PharmGKB |
| phenotype | string | Optional | Predicted metabolizer status or phenotype from diplotype | `Poor Metabolizer`, `Normal Metabolizer`, `Ultrarapid Metabolizer` | CPIC, DPWG, PharmGKB |
| recommendation | string | Optional | Dosing guidance or alternative therapy recommendation | `Avoid use`, `Use standard dose`, `Consider 50% dose reduction` | CPIC, DPWG, PharmGKB |
| evidence_level | string | Optional | Classification of evidence supporting recommendation | `A`, `B`, `Strong`, `Moderate`, `1A`, `1B` | CPIC, DPWG, PharmGKB |

---

## CPIC-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| cpic_guideline_id | string | Optional | CPIC unique guideline identifier | guideline_id |
| cpic_publication_doi | string | Optional | DOI reference for CPIC guideline | publication_doi |
| cpic_level | string | Optional | CPIC level: `A`, `B`, `C`, `D` | cpic_level |
| allele | string | Optional | Star allele (e.g., *1, *4) | allele |
| allele_function | string | Optional | Allele function: `Normal`, `Decreased`, `No function` | function |
| activity_score | float | Optional | Numeric activity value (0, 0.5, 1, 1.5, 2, etc.) | activity_score |
| recommendation_strength | string | Optional | Recommendation strength: `Strong`, `Moderate`, `Optional` | strength |
| clinical_implications | string | Optional | Clinical implications of phenotype | implications |

---

## DPWG-Specific Fields

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| dpwg_id | string | Optional | DPWG guideline identifier | dpwg_id |
| atc_code | string | Optional | WHO ATC drug classification code | atc_code |
| dpwg_urgency | string | Optional | Urgency classification: `Urgent`, `Non-urgent` | urgency |
| diplotype | string | Optional | Allele combination (e.g., *1/*4) | diplotype |

---

## PharmGKB-Specific Fields

### Identifiers

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| pharmgkb_accession | string | Optional | PharmGKB PA + integer identifier | PharmGKB_Accession |
| annotation_id | integer | Optional | Clinical annotation identifier | Annotation_ID |

### Gene Properties

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| is_vip_gene | boolean | Optional | Very Important Pharmacogene status | Is_VIP |
| has_variant_annotation | boolean | Optional | Has curated variants | Has_Variant_Annotation |

### Drug Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| phenotype_category | string | Optional | Category: `Efficacy`, `Toxicity`, `Dosage`, `Metabolism` | Phenotype_Category |
| trade_names | array[string] | Optional | Brand names for drug | Trade_Names |
| atc_codes | array[string] | Optional | ATC classification codes | ATC_Codes |
| rxnorm_ids | array[string] | Optional | RxNorm CUI identifiers | RxNorm_Identifiers |
| drugbank_id | string | Optional | DrugBank identifier (DB#####) | DrugBank_ID |
| has_pgx_on_label | boolean | Optional | FDA label contains PGx information | Has_PGx_On_Label |

### Variant Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| variant_haplotypes | array[string] | Optional | Associated haplotypes | Variant_Haplotypes |
| population_frequencies | object | Optional | Biogeographical population frequencies | Biogeographical_Groups |

### Pathway Information

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| pathway_id | string | Optional | Drug pathway identifier | Pathway_ID |
| pk_pd_type | string | Optional | Pharmacokinetics or pharmacodynamics | PK_PD |

---

## PharmVar-Specific Fields

### Allele Nomenclature

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| pharmvar_allele_name | string | Optional | Star notation (e.g., *4) | allele_name |
| pharmvar_core_allele | string | Optional | Core allele with suballele (e.g., *4.001) | core_allele |

### Variant Definition

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| defining_variants | array[string] | Optional | List of defining variants | defining_variants |
| reference_sequence_id | string | Optional | NG reference sequence | ref_seq_id |

### HGVS Nomenclature

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| hgvs_genomic | string | Optional | Genomic HGVS notation | hgvs_g |
| hgvs_coding | string | Optional | Coding HGVS notation | hgvs_c |
| hgvs_protein | string | Optional | Protein HGVS notation | hgvs_p |

### Variant Identifiers

| Field Name | Data Type | Cardinality | Description | Original Field |
|------------|-----------|-------------|-------------|----------------|
| rsid | string | Optional | dbSNP RS ID | rsid |
| variant_impact | string | Optional | Variant effect type | impact |

---

## Field Mapping Reference

### CPIC Mappings

| Original Field | Unified Field |
|----------------|---------------|
| guideline_id | cpic_guideline_id |
| publication_doi | cpic_publication_doi |
| cpic_level | cpic_level |
| allele | allele |
| function | allele_function |
| activity_score | activity_score |
| strength | recommendation_strength |
| implications | clinical_implications |

### DPWG Mappings

| Original Field | Unified Field |
|----------------|---------------|
| dpwg_id | dpwg_id |
| atc_code | atc_code |
| urgency | dpwg_urgency |
| diplotype | diplotype |
| activity_score | activity_score |

### PharmGKB Mappings

| Original Field | Unified Field |
|----------------|---------------|
| PharmGKB_Accession | pharmgkb_accession |
| Is_VIP | is_vip_gene |
| Has_Variant_Annotation | has_variant_annotation |
| Annotation_ID | annotation_id |
| Level_of_Evidence | evidence_level |
| Phenotype_Category | phenotype_category |
| Trade_Names | trade_names |
| ATC_Codes | atc_codes |
| RxNorm_Identifiers | rxnorm_ids |
| DrugBank_ID | drugbank_id |
| Has_PGx_On_Label | has_pgx_on_label |
| Variant_Haplotypes | variant_haplotypes |
| Biogeographical_Groups | population_frequencies |
| Pathway_ID | pathway_id |
| PK_PD | pk_pd_type |

### PharmVar Mappings

| Original Field | Unified Field |
|----------------|---------------|
| allele_name | pharmvar_allele_name |
| core_allele | pharmvar_core_allele |
| activity_value | activity_score |
| defining_variants | defining_variants |
| ref_seq_id | reference_sequence_id |
| hgvs_g | hgvs_genomic |
| hgvs_c | hgvs_coding |
| hgvs_p | hgvs_protein |
| rsid | rsid |
| impact | variant_impact |

---

## Evidence Level Reference

### CPIC Levels

| Level | Description |
|-------|-------------|
| A | Prescribing action recommended, high-quality evidence |
| B | Prescribing action recommended, moderate-quality evidence |
| C | No prescribing action at this time, moderate-quality evidence |
| D | Insufficient evidence, no prescribing recommendation |

### PharmGKB Evidence Levels

| Level | Description |
|-------|-------------|
| 1A | Annotation is based on variant-drug combination in a clinical PGx guideline |
| 1B | Strong evidence with clinical implementation |
| 2A | Moderate evidence, from multiple smaller studies |
| 2B | Moderate evidence, limited replication |
| 3 | Low-level evidence, single significant study |
| 4 | In vitro, case reports, or conflicting evidence |

### Phenotype Categories

| Phenotype | Activity Score (CYP2D6) | Description |
|-----------|------------------------|-------------|
| Ultrarapid Metabolizer | >2.25 | Very high enzyme activity |
| Normal Metabolizer | 1.25-2.25 | Normal enzyme activity |
| Intermediate Metabolizer | 0.25-1.25 | Reduced enzyme activity |
| Poor Metabolizer | 0 | No enzyme activity |

---

## Common Pharmacogenes

| Gene | Function | Common Drugs Affected |
|------|----------|----------------------|
| CYP2D6 | Drug metabolism | Codeine, tramadol, tamoxifen, antidepressants |
| CYP2C19 | Drug metabolism | Clopidogrel, PPIs, voriconazole |
| CYP2C9 | Drug metabolism | Warfarin, phenytoin, NSAIDs |
| CYP3A5 | Drug metabolism | Tacrolimus, cyclosporine |
| SLCO1B1 | Drug transport | Simvastatin, other statins |
| VKORC1 | Warfarin target | Warfarin |
| TPMT | Drug metabolism | Thiopurines (azathioprine, 6-MP) |
| DPYD | Drug metabolism | Fluoropyrimidines (5-FU, capecitabine) |
| UGT1A1 | Drug metabolism | Irinotecan, atazanavir |

---

*Generated: 2026-01-24*
