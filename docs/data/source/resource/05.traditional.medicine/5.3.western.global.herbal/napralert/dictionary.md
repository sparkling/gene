# NAPRALERT - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | napralert |
| **Name** | NAPRALERT - Natural Products Alert |
| **Total Fields** | 30 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| organism_id | String | Yes | NAPRALERT organism ID | NAP_ORG_001234 |
| scientific_name | String | No | Latin binomial name | Ginkgo biloba |
| family | String | No | Taxonomic family | Ginkgoaceae |
| kingdom | Enum | No | Taxonomic kingdom | Plantae |
| compound_id | String | No | Compound identifier | NAP_CPD_001234 |
| compound_name | String | No | Chemical name | Ginkgolide B |
| cas_number | String | No | CAS Registry Number | 15291-77-7 |

---

## Ethnobotany Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| ethnobotany.traditional_use | String | Documented therapeutic use |
| ethnobotany.geographic_origin | String | Region/country of use |
| ethnobotany.preparation_method | String | Traditional preparation |
| ethnobotany.administration_route | String | Route of administration |
| ethnobotany.dosage | String | Traditional dosage |

---

## Pharmacology Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| pharmacology.activity_type | String | Type of biological activity |
| pharmacology.model_system | String | In vitro/in vivo model |
| pharmacology.activity_value | Number | Quantitative measurement |
| pharmacology.activity_unit | String | Unit of measurement |
| pharmacology.reference_standard | String | Positive control |

---

## Literature Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| literature.pubmed_id | Integer | PubMed reference |
| literature.napralert_ref_id | String | Internal reference ID |
| literature.authors | String | Author names |
| literature.title | String | Article title |
| literature.journal | String | Journal name |
| literature.year | Integer | Publication year |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Organism ID | NAP_ORG_###### | NAPRALERT organism | NAP_ORG_001234 |
| Compound ID | NAP_CPD_###### | NAPRALERT compound | NAP_CPD_001234 |
| CAS Number | ###-##-# | CAS Registry | 458-37-7 |
| PubMed ID | Numeric | Literature reference | 12345678 |

---

## Enumerations

### Kingdoms

| Value | Description |
|-------|-------------|
| Plantae | Plants |
| Fungi | Fungi and mushrooms |
| Bacteria | Bacteria |
| Animalia | Animals (marine, insects) |

### Data Categories

| Category | Description |
|----------|-------------|
| Ethnobotany | Traditional use data |
| Pharmacology | Activity data |
| Toxicology | Safety data |
| Chemistry | Structural data |

---

## Entity Relationships

### Organism to Compound
- **Cardinality:** 1:N
- **Description:** Organisms contain compounds
- **Key Fields:** organism_id, compound_id

### Compound to Activity
- **Cardinality:** 1:N
- **Description:** Compounds have activity data
- **Key Fields:** compound_id, pharmacology

### Data to Literature
- **Cardinality:** N:1
- **Description:** Data linked to citations
- **Key Fields:** *, literature.pubmed_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NAPRALERT | Natural Products Alert | Database name |
| CAS | Chemical Abstracts Service | Chemical registry |
| UIC | University of Illinois at Chicago | Host institution |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
