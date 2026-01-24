# 6.2 Dietary Supplements - Data Dictionary

## Overview

This data dictionary documents all fields for the Dietary Supplements subcategory, integrating data from ConsumerLab, DSLD (Dietary Supplement Label Database), and Natural Medicines.

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 6.2 |
| Subcategory Name | Dietary Supplements |
| Data Sources | ConsumerLab, DSLD, Natural Medicines |
| Schema ID | https://gene.taxonomy/schemas/6.2-dietary-supplements |

---

## Unified Fields

These fields have been harmonized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| product_id | string | Required (1:1) | Unique identifier for the supplement product | `CL-2024-VIT-D-001`, `123456`, `NM-001234` |
| product_name | string | Required (1:1) | Name of the supplement product | `Vitamin D3 5000 IU`, `Turmeric` |
| brand_name | string | Optional (1:1) | Brand name of the supplement | `Nature Made`, `Nature's Way` |
| upc | string | Optional (1:1) | Universal Product Code barcode | `012345678901`, `033674100233` |
| category | array[string] | Optional (1:N) | Category classification of the supplement | `["Vitamins & Minerals"]`, `["Herbs/Botanicals"]` |
| ingredients | array[object] | Optional (1:N) | List of active ingredients with amounts | See structure below |
| serving_size | string | Optional (1:1) | Amount per single serving | `1 softgel`, `2 capsules` |
| product_type | string | Optional (1:1) | Physical form of the supplement | `Softgel`, `Capsule`, `Tablet` |
| net_contents | string | Optional (1:1) | Net quantity/count in package | `120 Softgels`, `90 Tablets` |
| safety_grade | string | Optional (1:1) | Overall safety classification | `Likely Safe`, `Possibly Safe` |

### Ingredients Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| ingredient_name | string | Name of the active ingredient | `Vitamin D3 (as Cholecalciferol)` |
| amount | string | Amount of ingredient | `5000` |
| unit | string | Unit of measurement | `IU`, `mg`, `mcg` |
| daily_value_percent | string | Percentage of FDA recommended daily intake | `1250`, `100` |
| is_blend_component | boolean | Whether ingredient is part of a proprietary blend | `true`, `false` |
| blend_name | string | Name of proprietary blend if applicable | `Energy Blend` |
| source_name | string | Source organism or material | `Lanolin`, `Fish oil` |

---

## Source-Specific Fields

### ConsumerLab

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| approval_status | string | Optional (1:1) | Quality test result status | `Approved`, `Not Approved`, `Voluntarily Withdrawn` |
| test_results | object | Optional (1:1) | Comprehensive test result data including potency, purity, contaminants | See structure below |
| heavy_metals | object | Optional (1:1) | Heavy metal contamination test results (lead, arsenic, cadmium, mercury) | See structure below |
| disintegration | object | Optional (1:1) | Tablet/capsule breakdown test results | See structure below |
| evidence_summary | object | Optional (1:N) | Evidence ratings for health claims by condition | `{"bone_health": "Strong evidence"}` |
| drug_interactions | array[string] | Optional (1:N) | Known drug interaction information | `["Corticosteroids", "Orlistat"]` |

#### Test Results Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| claimed_amount | string | Amount stated on label | `5000 IU` |
| actual_amount | string | Amount found in testing | `5120 IU` |
| percent_of_claim | number | Percentage of claimed amount | `102.4` |

#### Heavy Metals Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| lead_mcg_daily | number/string | Lead content per daily serving | `0.3`, `ND` |
| arsenic_mcg_daily | number/string | Arsenic content per daily serving | `0.1`, `ND` |
| cadmium_mcg_daily | number/string | Cadmium content per daily serving | `0.05`, `ND` |
| mercury_mcg_daily | number/string | Mercury content per daily serving | `0.02`, `ND` |

#### Disintegration Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| time_minutes | number | Time to disintegrate | `18` |
| limit_minutes | number | Maximum allowed time | `30` |
| status | string | Pass/Fail result | `Pass`, `Fail` |

### DSLD

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| dsld_id | integer | Optional (1:1) | DSLD internal unique identifier | `123456` |
| dietary_claims | array[string] | Optional (1:N) | Health/dietary claims on label | `["Gluten Free", "No artificial colors"]` |
| nhp_id | string | Optional (1:1) | Canadian Natural Health Product ID if applicable | `NPN12345678` |
| is_blend_component | boolean | Optional (1:1) | Whether ingredient is part of a proprietary blend | `true`, `false` |
| blend_name | string | Optional (1:1) | Name of proprietary blend if applicable | `Energy Blend`, `Immune Support Complex` |
| source_name | string | Optional (1:1) | Source organism or material for ingredient | `Lanolin`, `Fish oil` |
| manufacturer | object | Optional (1:N) | Manufacturer/distributor information | See structure below |

#### Manufacturer Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| name | string | Company name | `Nature's Way Products` |
| city | string | City location | `Green Bay` |
| state | string | State/Province | `WI` |
| country | string | Country | `USA` |

### Natural Medicines

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| monograph_id | string | Optional (1:1) | Unique identifier for ingredient monograph | `NM-001234` |
| scientific_names | array[string] | Optional (1:N) | Taxonomic names for the ingredient source | `["Curcuma longa", "Curcuma domestica"]` |
| common_names | array[string] | Optional (1:N) | Alternative/regional names for the ingredient | `["Curcumin", "Indian Saffron", "Haridra"]` |
| effectiveness_rating | string | Optional (1:N) | Evidence-based efficacy assessment by condition | `Effective`, `Likely Effective`, `Possibly Effective` |
| evidence_level | string | Optional (1:1) | Quality of supporting research | `A`, `B`, `C`, `D` |
| mechanism_of_action | string | Optional (1:1) | Pharmacological basis of action | `Curcumin inhibits COX-2, LOX, NF-kB pathways` |
| pharmacokinetics | object | Optional (1:1) | Absorption, metabolism, and elimination data | See structure below |
| interaction_severity | string | Optional (1:1) | Clinical significance of drug-supplement interaction | `Contraindicated`, `Major`, `Moderate`, `Minor` |
| interaction_mechanism | string | Optional (1:1) | Mechanism of drug interaction | `CYP3A4 induction`, `CYP2C9 inhibition` |
| cas_number | string | Optional (1:1) | Chemical Abstracts Service registry number | `458-37-7` |
| unii | string | Optional (1:1) | FDA Unique Ingredient Identifier | `IT942ZTH98` |

#### Pharmacokinetics Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| absorption | string | Absorption characteristics | `Poor oral bioavailability` |
| half_life | string | Elimination half-life | `1-2 hours` |
| metabolism | string | Metabolic pathway | `Hepatic glucuronidation` |

---

## Field Mappings

### ConsumerLab to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| review_id | product_id |
| product_name | product_name |
| brand | brand_name |
| upc | upc |
| category | category |
| subcategory | category |
| ingredient | ingredients[].ingredient_name |
| claimed_amount | ingredients[].amount |
| approval_status | approval_status |
| test_results | test_results |
| heavy_metals | heavy_metals |
| disintegration | disintegration |
| evidence_summary | effectiveness_ratings |
| drug_interactions | drug_interactions |

### DSLD to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| dsld_id | product_id |
| product_name | product_name |
| brand_name | brand_name |
| upc | upc |
| ingredient_group | category |
| ingredient_name | ingredients[].ingredient_name |
| amount | ingredients[].amount |
| unit | ingredients[].unit |
| daily_value_percent | ingredients[].daily_value_percent |
| serving_size | serving_size |
| product_type | product_type |
| net_contents | net_contents |
| is_blend_component | ingredients[].is_blend_component |
| blend_name | ingredients[].blend_name |
| source_name | ingredients[].source_name |
| dietary_claims | dietary_claims |
| manufacturer | manufacturer |

### Natural Medicines to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| monograph_id | product_id |
| name | product_name |
| category | category |
| scientific_names | scientific_names |
| common_names | common_names |
| safety_grade | safety_grade |
| effectiveness_rating | effectiveness_ratings |
| evidence_level | effectiveness_ratings[].evidence_level |
| mechanism_of_action | mechanism_of_action |
| pharmacokinetics | pharmacokinetics |
| interaction_severity | drug_interactions[].severity |
| interaction_mechanism | drug_interactions[].mechanism |
| cas_number | cas_number |
| unii | unii |

---

## Cardinality Legend

| Symbol | Meaning |
|--------|---------|
| 1:1 | Exactly one value per record |
| 1:N | One or more values per record |
| Required | Field must be present |
| Optional | Field may be null or absent |

---

## Data Quality Notes

1. **product_id** - Format varies by source: ConsumerLab uses review IDs, DSLD uses numeric IDs, Natural Medicines uses monograph IDs
2. **approval_status** - Only available from ConsumerLab quality testing
3. **safety_grade** - Uses standardized scale: `Likely Safe`, `Possibly Safe`, `Possibly Unsafe`, `Likely Unsafe`
4. **effectiveness_rating** - Natural Medicines uses evidence-based scale: `Effective`, `Likely Effective`, `Possibly Effective`, `Insufficient Evidence`
5. **evidence_level** - Letter grades A-D indicating quality of supporting research
6. **heavy_metals** - Values may be numeric (in mcg) or `ND` (Not Detected)
7. **drug_interactions** - ConsumerLab provides list format; Natural Medicines provides severity and mechanism details
