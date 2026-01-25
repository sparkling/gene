# ETCM - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | etcm |
| **Name** | Encyclopedia of Traditional Chinese Medicine |
| **Total Fields** | 28 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| herb_id | String | Yes | ETCM herb identifier | ETCM_H00001 |
| herb_name_en | String | Yes | English herb name | Ginseng |
| herb_name_cn | String | No | Chinese herb name | 人参 |
| herb_name_pinyin | String | No | Pinyin romanization | renshen |
| scientific_name | String | No | Latin botanical name | Panax ginseng |
| compound_id | String | No | ETCM compound identifier | ETCM_C00001 |
| compound_name | String | No | Compound name | Ginsenoside Rg1 |
| target_id | String | No | Target protein identifier | P04637 |
| gene_symbol | String | No | Gene symbol | TP53 |

---

## TCM Classification Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| tcm_class | String | Primary TCM classification | Qi-tonifying herbs |
| pharmacopoeia_status | String | Chinese Pharmacopoeia status | Official |
| plant_family | String | Botanical family | Araliaceae |
| used_part | String | Medicinal part used | Root |
| nature | String | TCM temperature nature | Warm |
| flavor | Array | TCM taste properties | ["Sweet", "Slightly bitter"] |
| channel_tropism | Array | Target meridians | ["Spleen", "Lung", "Heart"] |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Herb ID | ETCM_H##### | Herb identifier | ETCM_H00001 |
| Compound ID | ETCM_C##### | Compound identifier | ETCM_C00001 |
| Target ID | ETCM_T##### | Target identifier | ETCM_T00001 |
| PubChem CID | Numeric | Cross-reference | 441923 |

---

## Enumerations

### TCM Herb Classes

| Value | Description |
|-------|-------------|
| Qi-tonifying | Supplements Qi energy |
| Blood-tonifying | Nourishes blood |
| Yin-tonifying | Nourishes Yin |
| Yang-tonifying | Warms Yang |
| Heat-clearing | Clears pathogenic heat |
| Exterior-releasing | Releases exterior patterns |
| Dampness-draining | Removes dampness |

### Pharmacopoeia Status

| Value | Description |
|-------|-------------|
| Official | Listed in Chinese Pharmacopoeia |
| Not official | Not in Pharmacopoeia |
| Regional | Regional/provincial standard |

---

## Entity Relationships

### Herb to Compound
- **Cardinality:** 1:N
- **Description:** Each herb contains multiple compounds
- **Key Fields:** herb_id, compound_id

### Compound to Target
- **Cardinality:** N:M
- **Description:** Compounds interact with targets
- **Key Fields:** compound_id, target_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ETCM | Encyclopedia of Traditional Chinese Medicine | Database name |
| TCM | Traditional Chinese Medicine | Medical system |
| CP | Chinese Pharmacopoeia | Official standard |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
