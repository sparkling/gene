# BATMAN-TCM 2.0 - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | batman.tcm |
| **Name** | BATMAN-TCM 2.0 |
| **Total Fields** | 35 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| ingredient_id | String | Yes | Ingredient/compound identifier | BATMAN-TCM-ING-5679 |
| compound_name | String | No | Compound name | Ginsenoside Rg1 |
| target_id | String | Yes | Target protein identifier (UniProt) | P35354 |
| gene_symbol | String | No | Target gene symbol | PTGS2 |
| interaction_type | Enum | Yes | Known or predicted interaction | known, predicted |
| confidence_score | Number | No | Interaction confidence (0-1) | 0.85 |

---

## TCM-Specific Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| herb.herb_name_chinese | String | Chinese herb name | 人参 |
| herb.herb_name_pinyin | String | Pinyin romanization | renshen |
| herb.herb_name_latin | String | Latin botanical name | Panax ginseng |
| herb.properties_tcm.nature | Enum | TCM temperature nature | cold, cool, neutral, warm, hot |
| herb.properties_tcm.flavor | Array | TCM taste classifications | ["sweet", "bitter"] |
| herb.properties_tcm.meridian_tropism | Array | Target meridians | ["Lung", "Spleen"] |
| formula.formula_name_chinese | String | Formula Chinese name | 四君子汤 |
| formula.formula_name_pinyin | String | Formula Pinyin | sijunzitang |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Ingredient ID | BATMAN-TCM-ING-#### | Compound identifier | BATMAN-TCM-ING-5679 |
| TTI ID | BATMAN-TCM-####-###### | Interaction ID | BATMAN-TCM-2024-001234 |
| UniProt | Alphanumeric | Protein accession | P35354 |
| PubChem CID | Numeric | Compound ID | 441923 |

---

## Enumerations

### Interaction Type

| Value | Description |
|-------|-------------|
| known | Experimentally validated interaction |
| predicted | Computationally predicted interaction |

### TCM Nature (Si Qi)

| Value | Description |
|-------|-------------|
| cold | Han - Cold nature |
| cool | Liang - Cool nature |
| neutral | Ping - Neutral |
| warm | Wen - Warm nature |
| hot | Re - Hot nature |

### TCM Flavor (Wu Wei)

| Value | Description |
|-------|-------------|
| sweet | Gan - Sweet |
| bitter | Ku - Bitter |
| sour | Suan - Sour |
| pungent | Xin - Pungent/Acrid |
| salty | Xian - Salty |

---

## Entity Relationships

### Herb to Ingredient
- **Cardinality:** 1:N
- **Description:** Each herb contains multiple compounds
- **Key Fields:** herb.herb_id, ingredient_id

### Ingredient to Target
- **Cardinality:** N:M
- **Description:** Compounds interact with multiple targets
- **Key Fields:** ingredient_id, target_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| BATMAN | Bioinformatics Analysis Tool for Molecular mechANism | Database name |
| TCM | Traditional Chinese Medicine | Medical system |
| TTI | Target-TCM Interaction | Interaction type |
| Si Qi | Four Natures | TCM classification |
| Wu Wei | Five Flavors | TCM classification |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
