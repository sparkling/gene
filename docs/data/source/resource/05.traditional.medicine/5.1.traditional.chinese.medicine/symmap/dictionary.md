# SymMap - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | symmap |
| **Name** | SymMap TCM Symptom Mapping |
| **Total Fields** | 30 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| symptom_id | String | Yes | SymMap symptom identifier | SMSY000001 |
| symptom_name | String | Yes | TCM symptom name | Fatigue |
| symptom_name_cn | String | No | Chinese symptom name | 乏力 |
| herb_id | String | No | Associated herb ID | SMHB000001 |
| herb_name | String | No | Herb name | Ginseng |
| compound_id | String | No | Compound identifier | SMCP000001 |
| target_id | String | No | Target identifier | SMTG000001 |
| gene_symbol | String | No | Gene symbol | EGFR |

---

## TCM Symptom Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| tcm_syndrome | String | TCM syndrome/pattern | Qi deficiency |
| modern_disease | String | Western disease mapping | Chronic fatigue syndrome |
| symptom_category | String | Symptom classification | Constitutional |
| body_system | String | Affected body system | General |
| herb_efficacy | String | How herb treats symptom | Tonifies Qi |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Symptom ID | SMSY###### | Symptom identifier | SMSY000001 |
| Herb ID | SMHB###### | Herb identifier | SMHB000001 |
| Compound ID | SMCP###### | Compound identifier | SMCP000001 |
| Target ID | SMTG###### | Target identifier | SMTG000001 |

---

## Enumerations

### TCM Syndromes

| Value | Description |
|-------|-------------|
| Qi deficiency | Insufficient Qi energy |
| Blood stasis | Blood circulation impairment |
| Yin deficiency | Insufficient Yin fluids |
| Yang deficiency | Insufficient Yang warmth |
| Phlegm-dampness | Phlegm accumulation |
| Heat toxin | Pathogenic heat |

### Symptom Categories

| Value | Description |
|-------|-------------|
| Constitutional | General body symptoms |
| Respiratory | Lung/breathing symptoms |
| Digestive | GI tract symptoms |
| Cardiovascular | Heart/circulation symptoms |
| Neurological | Nervous system symptoms |
| Musculoskeletal | Bone/muscle symptoms |

---

## Entity Relationships

### Symptom to Herb
- **Cardinality:** N:M
- **Description:** Symptoms treated by multiple herbs
- **Key Fields:** symptom_id, herb_id

### Herb to Compound
- **Cardinality:** 1:N
- **Description:** Herbs contain compounds
- **Key Fields:** herb_id, compound_id

### Compound to Target
- **Cardinality:** N:M
- **Description:** Compounds affect targets
- **Key Fields:** compound_id, target_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| SymMap | Symptom Mapping | Database name |
| TCM | Traditional Chinese Medicine | Medical system |
| SMSY | SymMap Symptom | ID prefix |
| SMHB | SymMap Herb | ID prefix |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
