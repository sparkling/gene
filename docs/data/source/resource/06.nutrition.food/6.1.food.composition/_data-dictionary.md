# 6.1 Food Composition - Data Dictionary

## Overview

This data dictionary documents all fields for the Food Composition subcategory, integrating data from Open Food Facts and FooDB.

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 6.1 |
| Subcategory Name | Food Composition |
| Data Sources | Open Food Facts, FooDB |
| Schema ID | https://gene.taxonomy/schemas/6.1-food-composition |

---

## Unified Fields

These fields have been harmonized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| food_id | string | Required (1:1) | Primary identifier for a food item | `3017620425035`, `FOOD00234` |
| food_name | string | Required (1:1) | Common name of the food product | `Nutella hazelnut spread`, `Apple, raw, unpeeled` |
| scientific_name | string | Optional (1:1) | Scientific/Latin name of the food source organism | `Malus domestica`, `Angelica keiskei` |
| food_category | array[string] | Optional (N:M) | Category classification of the food | `["en:breakfast-cereals"]`, `["Herbs and Spices"]` |
| brand | array[string] | Optional (1:N) | Brand name(s) of the product | `["Ferrero"]`, `["Coca-Cola"]` |
| nutrient_values | object | Optional (1:N) | Nutritional content per 100g | `{"energy_kcal": 539, "protein_g": 6.3}` |
| compound_content | array[object] | Optional (1:N) | Chemical compound content in the food | `[{"compound_name": "Quercetin", "content": 4.01, "unit": "mg/100g"}]` |

### Nutrient Values Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| energy_kcal | number | Energy content in kilocalories | `539` |
| protein_g | number | Protein content in grams | `6.3` |
| fat_g | number | Fat content in grams | `30.7` |
| carbs_g | number | Carbohydrate content in grams | `57.5` |
| fiber_g | number | Fiber content in grams | `3.4` |
| sugar_g | number | Sugar content in grams | `56.3` |
| sodium_mg | number | Sodium content in milligrams | `41` |

### Compound Content Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| compound_name | string | Name of the chemical compound | `Quercetin` |
| compound_id | string | Unique identifier for the compound | `FDB000001` |
| content | number | Concentration amount | `4.01` |
| unit | string | Unit of measurement | `mg/100g` |

---

## Source-Specific Fields

### Open Food Facts

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| nutriscore_grade | string | Optional (1:1) | Nutri-Score nutritional grade from A (best) to E (worst) | `a`, `b`, `c`, `d`, `e` |
| nova_group | integer | Optional (1:1) | NOVA food processing classification (1-4, 4=ultra-processed) | `1`, `2`, `3`, `4` |
| ecoscore_grade | string | Optional (1:1) | Environmental impact grade from A (best) to E (worst) | `a`, `b`, `c`, `d`, `e` |
| quantity | string | Optional (1:1) | Package quantity with unit | `400 g`, `1 L` |

### FooDB

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| compound_id | string | Optional (1:1) | Unique identifier for chemical compounds (FDB prefix) | `FDB000001` |
| moldb_smiles | string | Optional (1:1) | SMILES notation representing molecular structure | `CC(=O)Oc1ccccc1C(=O)O` |
| moldb_inchikey | string | Optional (1:1) | InChI key - unique chemical structure identifier | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` |
| cas_number | string | Optional (1:1) | Chemical Abstracts Service registry number | `7084-24-4` |
| itis_id | integer | Optional (1:1) | Integrated Taxonomic Information System ID for species | `18098` |
| export_status | integer | Optional (1:1) | Compound status indicator (1=active, 2=inactive) | `1`, `2` |
| orig_content | decimal | Optional (1:1) | Measured or average concentration of compound in food | `25.3` |
| orig_unit | string | Optional (1:1) | Unit of measurement for compound concentration | `mg/100g`, `mcg/100g` |

---

## Field Mappings

### Open Food Facts to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| code | food_id |
| barcode | food_id |
| product_name | food_name |
| categories_tags | food_category |
| brands | brand |
| nutriscore_grade | nutriscore_grade |
| nova_group | nova_group |
| ecoscore_grade | ecoscore_grade |
| quantity | quantity |
| nutriments | nutrient_values |

### FooDB to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| food_id | food_id |
| public_id | food_id |
| name | food_name |
| name_scientific | scientific_name |
| food_group | food_category |
| compound_id | compound_content[].compound_id |
| moldb_smiles | moldb_smiles |
| moldb_inchikey | moldb_inchikey |
| cas_number | cas_number |
| itis_id | itis_id |
| orig_content | compound_content[].content |
| orig_unit | compound_content[].unit |

---

## Cardinality Legend

| Symbol | Meaning |
|--------|---------|
| 1:1 | Exactly one value per record |
| 1:N | One or more values per record |
| N:M | Many-to-many relationship |
| Required | Field must be present |
| Optional | Field may be null or absent |

---

## Data Quality Notes

1. **food_id** - Open Food Facts uses barcode/EAN codes; FooDB uses internal identifiers with `FOOD` prefix
2. **nutriscore_grade** and **ecoscore_grade** - Lowercase letters a-e only from Open Food Facts
3. **nova_group** - Integer values 1-4 only; 4 indicates ultra-processed foods
4. **compound_content** - FooDB provides extensive compound data; Open Food Facts does not
5. **scientific_name** - Only available from FooDB for biological food sources
