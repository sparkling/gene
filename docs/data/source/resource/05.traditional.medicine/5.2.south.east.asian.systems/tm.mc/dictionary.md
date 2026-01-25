# TM-MC - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | tm.mc |
| **Name** | TM-MC Traditional Medicine Multi-Cultural |
| **Total Fields** | 32 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| plant_id | String | Yes | TM-MC plant identifier | TMMC_P_001234 |
| scientific_name | String | No | Latin botanical name | Curcuma longa |
| family | String | No | Botanical family | Zingiberaceae |
| system_code | Enum | No | Traditional system code | AYU, TCM, UNA |
| local_name | String | No | Name in local language | Haridra |
| compound_id | String | No | Compound identifier | TMMC_C_000123 |
| compound_name | String | No | Chemical name | Curcumin |

---

## Multi-System Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| traditional_system | Array | Systems using this plant |
| convergent_uses | Array | Shared therapeutic uses |
| cross_system_validation | Object | Multi-system validation data |
| traditional_terms | Array | System-specific terminology |

---

## Ayurveda-Specific Fields

| Field Name | Data Type | Description | Examples |
|------------|-----------|-------------|----------|
| ayurveda_properties.rasa | Array | Taste | ["tikta", "katu"] |
| ayurveda_properties.guna | Array | Quality | ["laghu", "ruksha"] |
| ayurveda_properties.vipaka | String | Post-digestive effect | katu |
| ayurveda_properties.virya | String | Potency | ushna |
| ayurveda_properties.dosha_effect | Array | Dosha balancing | ["kapha-vata shamaka"] |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Plant ID | TMMC_P_###### | Plant identifier | TMMC_P_001234 |
| Compound ID | TMMC_C_###### | Compound identifier | TMMC_C_000123 |
| Indication ID | TMMC_I_###### | Indication identifier | TMMC_I_000456 |
| PubChem CID | Numeric | Cross-reference | 969516 |

---

## Enumerations

### Traditional System Codes

| Code | System | Region |
|------|--------|--------|
| AYU | Ayurveda | India |
| TCM | Traditional Chinese Medicine | China |
| UNA | Unani | Middle East/South Asia |
| SID | Siddha | Tamil Nadu, India |
| KAM | Kampo | Japan |
| JAM | Jamu | Indonesia |

### Ayurvedic Rasa (Taste)

| Value | Translation |
|-------|-------------|
| madhura | Sweet |
| amla | Sour |
| lavana | Salty |
| katu | Pungent |
| tikta | Bitter |
| kashaya | Astringent |

### Ayurvedic Virya (Potency)

| Value | Translation |
|-------|-------------|
| ushna | Hot |
| sheeta | Cold |

---

## Entity Relationships

### Plant to Systems
- **Cardinality:** 1:N
- **Description:** Plants used in multiple traditional systems
- **Key Fields:** plant_id, traditional_system

### Plant to Compound
- **Cardinality:** 1:N
- **Description:** Plants contain compounds
- **Key Fields:** plant_id, compound_id

### Indication Convergence
- **Cardinality:** N:M
- **Description:** Shared uses across systems
- **Key Fields:** plant_id, convergent_uses

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| TM-MC | Traditional Medicine Multi-Cultural | Database name |
| AYU | Ayurveda | Indian system |
| TCM | Traditional Chinese Medicine | Chinese system |
| UNA | Unani | Greco-Arabic system |
| SID | Siddha | Tamil system |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
