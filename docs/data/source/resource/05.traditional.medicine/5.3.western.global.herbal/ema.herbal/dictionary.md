# EMA Herbal Medicines - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | ema.herbal |
| **Name** | EMA Herbal Medicines |
| **Total Fields** | 28 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| monograph_id | String | No | EMA monograph reference | EMA/HMPC/150848/2006 |
| herbal_substance | String | Yes | Latin binomial name | Valeriana officinalis L. |
| plant_part | String | Yes | Standardized plant part | radix |
| common_name | String | No | Common English name | Valerian root |
| pharmacopoeia_name | String | No | European Pharmacopoeia name | Valerianae radix |
| status | Enum | No | Regulatory status | well-established use |

---

## Clinical Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| therapeutic_area | Array | Therapeutic indication areas |
| indications | Array | Approved indications with evidence |
| posology.preparation | String | Pharmaceutical preparation |
| posology.single_dose | String | Single dose recommendation |
| posology.daily_dose | String | Daily dose recommendation |
| posology.duration | String | Treatment duration |
| contraindications | Array | Contraindications |
| special_warnings | Array | Warnings and precautions |
| interactions | Array | Drug interactions |
| adverse_reactions | Array | Known adverse reactions |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Monograph ID | EMA/HMPC/######/#### | EMA reference | EMA/HMPC/150848/2006 |
| Plant Part | Latin term | Standardized part | radix, folium, flos |

---

## Enumerations

### Regulatory Status

| Value | Description |
|-------|-------------|
| well-established use | Clinical evidence supports efficacy |
| traditional use | 30+ years of traditional use |
| under assessment | Currently being evaluated |

### Plant Parts (Latin)

| Latin | English |
|-------|---------|
| radix | Root |
| folium | Leaf |
| flos | Flower |
| herba | Herb (aerial parts) |
| cortex | Bark |
| fructus | Fruit |
| semen | Seed |
| rhizoma | Rhizome |

### Evidence Levels

| Level | Description |
|-------|-------------|
| Clinical trials | Randomized controlled trials |
| Clinical studies | Non-randomized clinical evidence |
| Traditional use | Historical use documentation |
| Expert opinion | HMPC expert assessment |

---

## Entity Relationships

### Substance to Indications
- **Cardinality:** 1:N
- **Description:** Each substance has multiple indications
- **Key Fields:** herbal_substance, indications

### Substance to Interactions
- **Cardinality:** 1:N
- **Description:** Known drug interactions
- **Key Fields:** herbal_substance, interactions

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| EMA | European Medicines Agency | Regulatory body |
| HMPC | Committee on Herbal Medicinal Products | EMA committee |
| EP | European Pharmacopoeia | Quality standard |
| WEU | Well-Established Use | Regulatory category |
| TU | Traditional Use | Regulatory category |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
