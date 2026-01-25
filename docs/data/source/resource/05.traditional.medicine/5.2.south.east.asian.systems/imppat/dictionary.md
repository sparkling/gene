# IMPPAT 2.0 - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | imppat |
| **Name** | IMPPAT 2.0 |
| **Total Fields** | 45 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| plant_id | String | Yes | IMPPAT plant identifier | IMPPAT_PLANT_03985 |
| scientific_name | String | No | Latin binomial name | Withania somnifera |
| family | String | No | Botanical family | Solanaceae |
| compound_id | String | No | Compound identifier | IMPPAT_CHEM_15234 |
| compound_name | String | No | Chemical name | Withanolide A |
| gene_symbol | String | No | Target gene symbol | NR3C1 |
| stitch_score | Integer | No | STITCH confidence (0-1000) | 850 |

---

## Vernacular Name Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| vernacular_names.hindi | String | Hindi name |
| vernacular_names.tamil | String | Tamil name |
| vernacular_names.telugu | String | Telugu name |
| vernacular_names.kannada | String | Kannada name |
| vernacular_names.malayalam | String | Malayalam name |
| vernacular_names.bengali | String | Bengali name |
| vernacular_names.gujarati | String | Gujarati name |
| vernacular_names.marathi | String | Marathi name |
| vernacular_names.punjabi | String | Punjabi name |
| vernacular_names.oriya | String | Oriya name |

---

## Drug-Likeness Fields

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| lipinski_rule_of_five.passes | Boolean | Passes Lipinski RO5 |
| lipinski_rule_of_five.violations | Integer | Number of violations |
| qed_score | Number | QED score (0-1) |
| ghose_filter.passes | Boolean | Passes Ghose filter |
| veber_filter.passes | Boolean | Passes Veber filter |

---

## ADMET Fields

| Category | Field | Description |
|----------|-------|-------------|
| Absorption | gi_absorption | GI absorption level |
| Absorption | pgp_substrate | P-glycoprotein substrate |
| Distribution | bbb_permeant | BBB permeability |
| Metabolism | cyp1a2_inhibitor | CYP1A2 inhibition |
| Metabolism | cyp3a4_inhibitor | CYP3A4 inhibition |
| Toxicity | ames_mutagenicity | Ames test result |
| Toxicity | hepatotoxicity | Liver toxicity prediction |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| Plant ID | IMPPAT_PLANT_##### | Plant identifier | IMPPAT_PLANT_03985 |
| Compound ID | IMPPAT_CHEM_##### | Compound identifier | IMPPAT_CHEM_15234 |
| InChIKey | 27-char string | Chemical hash | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| PubChem CID | Numeric | Cross-reference | 11294368 |

---

## Enumerations

### Traditional Systems

| Value | Description |
|-------|-------------|
| Ayurveda | Indian Ayurvedic medicine |
| Siddha | Tamil Siddha medicine |
| Unani | Greco-Arabic medicine |
| Homeopathy | Homeopathic medicine |

### IUCN Conservation Status

| Code | Status |
|------|--------|
| LC | Least Concern |
| NT | Near Threatened |
| VU | Vulnerable |
| EN | Endangered |
| CR | Critically Endangered |

---

## Entity Relationships

### Plant to Compound
- **Cardinality:** 1:N
- **Description:** Plants contain multiple phytochemicals
- **Key Fields:** plant_id, compound_id

### Compound to Target
- **Cardinality:** N:M
- **Description:** Compounds interact with targets (via STITCH)
- **Key Fields:** compound_id, gene_symbol

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IMPPAT | Indian Medicinal Plants, Phytochemistry And Therapeutics | Full name |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity | Property category |
| QED | Quantitative Estimate of Drug-likeness | Metric |
| BBB | Blood-Brain Barrier | Distribution barrier |
| STITCH | Search Tool for Interactions of Chemicals | Interaction source |

---

## See Also

- [schema.json](./schema.json) - Full schema documentation
- [sample.json](./sample.json) - Example records
