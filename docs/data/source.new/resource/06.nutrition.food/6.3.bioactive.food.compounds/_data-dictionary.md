# 6.3 Bioactive Food Compounds - Data Dictionary

## Overview

This data dictionary documents all fields for the Bioactive Food Compounds subcategory, containing data from the eBASIS (Bioactive Substances in Food Information System) database.

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 6.3 |
| Subcategory Name | Bioactive Food Compounds |
| Data Sources | eBASIS |
| Schema ID | https://gene.taxonomy/schemas/6.3-bioactive-food-compounds |

---

## Unified Fields

These fields represent the core data model for bioactive compound information.

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| component_id | integer | Required (1:1) | Unique identifier for bioactive compound | `1234` |
| component_name | string | Required (1:1) | Name of the bioactive compound | `Quercetin`, `Beta-carotene`, `Resveratrol` |
| component_group | string | Optional (1:1) | Major category of bioactive compounds | `Flavonoids`, `Carotenoids`, `Phenolic Acids` |
| component_subgroup | string | Optional (1:1) | Subcategory within component group | `Flavonols`, `Flavones`, `Anthocyanidins` |
| food_id | integer | Required (1:1) | Unique identifier for the food source | `101` |
| food_name | string | Required (1:1) | Common name of the food | `Apple, raw, with skin`, `Red onion, raw` |
| scientific_name | string | Optional (1:1) | Scientific name of the food source organism | `Malus domestica`, `Allium cepa` |
| composition_value | decimal | Required (1:1) | Concentration value of compound in food | `4.01`, `233.84` |
| composition_unit | string | Required (1:1) | Unit of measurement for concentration | `mg/100g`, `mcg/100g` |
| data_quality_score | integer | Optional (1:1) | EuroFIR quality rating (1-5, 5=excellent) | `1`, `2`, `3`, `4`, `5` |

---

## Source-Specific Fields (eBASIS)

### Chemical Identification

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| cas_number | string | Optional (1:1) | Chemical Abstracts Service registry number | `117-39-5` |
| pubchem_cid | integer | Optional (1:1) | PubChem compound identifier | `5280343` |
| molecular_formula | string | Optional (1:1) | Chemical formula of the compound | `C15H10O7` |
| molecular_weight | decimal | Optional (1:1) | Molecular weight in g/mol | `302.236` |
| iupac_name | string | Optional (1:1) | IUPAC systematic chemical name | `2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxychromen-4-one` |
| inchi | string | Optional (1:1) | International Chemical Identifier string | `InChI=1S/C15H10O7/c16-7-4-10(19)...` |
| inchikey | string | Optional (1:1) | Hashed InChI identifier | `REFJWTPEDVJJIY-UHFFFAOYSA-N` |
| synonyms | array[string] | Optional (1:N) | Alternative names for the compound | `["Sophoretin", "Meletin", "Quercetine"]` |

### Food Classification

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| langual_codes | array[string] | Optional (1:N) | LanguaL food description thesaurus codes | `["A0148", "B1001"]` |
| eurofir_code | string | Optional (1:1) | EuroFIR food code | `FR001234` |
| edible_portion | decimal | Optional (1:1) | Percentage of food that is consumable | `95.0` |

### Measurement Statistics

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| value_type | string | Optional (1:1) | Type of measurement value | `mean`, `median`, `single` |
| min_value | decimal | Optional (1:1) | Minimum observed concentration | `2.10` |
| max_value | decimal | Optional (1:1) | Maximum observed concentration | `6.45` |
| n_samples | integer | Optional (1:1) | Number of samples analyzed | `12` |
| std_dev | decimal | Optional (1:1) | Standard deviation of measurements | `1.23` |
| analytical_method | string | Optional (1:1) | Analytical technique used | `HPLC-DAD`, `GC-MS` |

### References

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| reference_id | integer | Optional (1:1) | Reference to source publication | `5678` |
| pubmed_id | integer | Optional (1:1) | PubMed reference ID | `32145678` |
| doi | string | Optional (1:1) | Digital Object Identifier for publication | `10.1016/j.jfca.2020.103567` |

---

## Field Mappings

### eBASIS to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| component_id | component_id |
| component_name | component_name |
| component_group | component_group |
| component_subgroup | component_subgroup |
| food_id | food_id |
| food_name | food_name |
| scientific_name | scientific_name |
| value | composition_value |
| unit | composition_unit |
| data_quality_score | data_quality_score |
| cas_number | cas_number |
| pubchem_cid | pubchem_cid |
| molecular_formula | molecular_formula |
| molecular_weight | molecular_weight |
| iupac_name | iupac_name |
| inchi | inchi |
| inchikey | inchikey |
| synonyms | synonyms |
| langual_code | langual_codes |
| eurofir_code | eurofir_code |
| edible_portion | edible_portion |
| value_type | value_type |
| min_value | min_value |
| max_value | max_value |
| n_samples | n_samples |
| std_dev | std_dev |
| analytical_method | analytical_method |
| reference_id | reference_id |
| pubmed_id | pubmed_id |
| doi | doi |

---

## Cardinality Legend

| Symbol | Meaning |
|--------|---------|
| 1:1 | Exactly one value per record |
| 1:N | One or more values per record |
| Required | Field must be present |
| Optional | Field may be null or absent |

---

## Component Groups Reference

| Group | Description | Example Compounds |
|-------|-------------|-------------------|
| Flavonoids | Polyphenolic compounds with flavone backbone | Quercetin, Kaempferol, Catechins |
| Carotenoids | Tetraterpenoid pigments | Beta-carotene, Lutein, Lycopene |
| Phenolic Acids | Aromatic acids with hydroxyl groups | Caffeic acid, Ferulic acid |
| Glucosinolates | Sulfur-containing compounds | Sulforaphane, Indole-3-carbinol |
| Phytosterols | Plant sterols | Beta-sitosterol, Campesterol |
| Terpenes | Hydrocarbon compounds | Limonene, Pinene |

## Component Subgroups Reference

| Parent Group | Subgroups |
|--------------|-----------|
| Flavonoids | Flavonols, Flavones, Flavanones, Flavanols, Anthocyanidins, Isoflavones |
| Carotenoids | Carotenes, Xanthophylls |
| Phenolic Acids | Hydroxybenzoic acids, Hydroxycinnamic acids |

---

## Data Quality Score Interpretation

| Score | Rating | Description |
|-------|--------|-------------|
| 5 | Excellent | Complete documentation, validated methods, multiple sources |
| 4 | Good | Well-documented, standard analytical methods |
| 3 | Acceptable | Adequate documentation, recognized methods |
| 2 | Limited | Incomplete documentation or limited validation |
| 1 | Minimal | Basic data only, minimal quality assurance |

---

## Analytical Methods Reference

| Method | Full Name | Typical Use |
|--------|-----------|-------------|
| HPLC-DAD | High-Performance Liquid Chromatography with Diode Array Detection | Flavonoids, phenolic acids |
| HPLC-MS | High-Performance Liquid Chromatography-Mass Spectrometry | Complex mixtures |
| GC-MS | Gas Chromatography-Mass Spectrometry | Volatile compounds, terpenes |
| LC-MS/MS | Liquid Chromatography-Tandem Mass Spectrometry | Trace compounds |
| UV-Vis | UV-Visible Spectrophotometry | Total phenolics, carotenoids |

---

## Data Quality Notes

1. **component_id** - Internal eBASIS identifier; unique per compound
2. **food_id** - Internal eBASIS identifier; unique per food item
3. **composition_value** - Primary value reported; check `value_type` for interpretation
4. **data_quality_score** - EuroFIR standardized quality assessment (1-5 scale)
5. **langual_codes** - LanguaL (Langua aLimentaria) international food description system
6. **analytical_method** - Important for data comparability across studies
7. **n_samples** - Sample size affects reliability of mean values
8. **inchikey** - Enables cross-referencing with other chemical databases
