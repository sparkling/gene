# 2.4 Food Compounds and Nutrients - Data Dictionary

## Overview

This data dictionary documents the unified schema for food compound and nutrient data from three major databases: Phenol-Explorer, PhytoHub, and USDA FoodData Central.

**Subcategory ID:** 2.4
**Subcategory Name:** Food Compounds and Nutrients
**Data Sources:** Phenol-Explorer, PhytoHub, USDA FoodData

---

## Unified Fields

### Core Identification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| compound_id | string/integer | 1:1, Required | Primary identifier for compound | `1`, `PHUB000001`, `167512` | Phenol-Explorer: compound_id, PhytoHub: phytohub_id, USDA: fdc_id |
| name | string | 1:1, Required | Compound or nutrient name | `Quercetin`, `Vitamin C` | Phenol-Explorer, PhytoHub, USDA FoodData |

### Structure Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| smiles | string | 1:1, Optional | SMILES structure notation | `OC1=CC(O)=C2C(=O)C(O)=C(OC2=C1)C1=CC=C(O)C(O)=C1` | Phenol-Explorer, PhytoHub |
| inchi_key | string | 1:1, Optional | InChI Key identifier | `REFJWTPEDVJJIY-UHFFFAOYSA-N` | PhytoHub |
| molecular_weight | decimal | 1:1, Optional | Molecular weight in Daltons | `302.24` | Phenol-Explorer, PhytoHub |

### Food Content Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| food_sources | object[] | 1:N, Optional | Foods containing the compound | `[{"food_name": "Onion", "food_id": "11282", "content_value": 20.3, "content_unit": "mg/100g"}]` | Phenol-Explorer, PhytoHub, USDA FoodData |
| content_value | decimal | 1:1, Optional | Amount in food | `20.3` | Phenol-Explorer, PhytoHub, USDA FoodData |
| content_unit | string | 1:1, Optional | Unit of measurement | `mg/100g`, `mg/100ml`, `ug/100g` | Phenol-Explorer, PhytoHub, USDA FoodData |

---

## Source-Specific Fields

### Phenol-Explorer

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| compound_class | string | 1:1, Optional | Polyphenol class | `Flavonoids`, `Phenolic acids`, `Stilbenes`, `Lignans` |
| subclass | string | 1:1, Optional | Polyphenol subclass | `Flavonols`, `Flavones`, `Isoflavones`, `Anthocyanins` |
| pubchem_cid | integer | 1:1, Optional | PubChem compound ID | `5280343` |
| cas_number | string | 1:1, Optional | CAS registry number | `117-39-5` |
| metabolites | object[] | 1:N, Optional | Human metabolites (Phase I/II, microbial) | `[{"name": "Quercetin-3-glucuronide", "type": "Phase II", "phase": "glucuronidation"}]` |
| pharmacokinetics | object | 1:1, Optional | Cmax, Tmax, AUC values | `{"cmax": 0.5, "tmax": 1.5, "auc": 2.3}` |

### PhytoHub

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| phytohub_id | string | 1:1, Required | PhytoHub identifier (PHUB prefix) | `PHUB000001` |
| compound_type | string | 1:1, Optional | Compound type classification | `parent`, `phase1`, `phase2`, `microbial` |
| parent_id | string | 1:1, Optional | Parent compound for metabolites | `PHUB000042` |
| transformation | string | 1:1, Optional | Metabolic transformation type | `glucuronidation`, `sulfation`, `methylation` |
| ms_spectra | object[] | 1:N, Optional | MS/MS reference spectra | `[{"mz": 301.035, "intensity": 100}, {"mz": 151.003, "intensity": 45}]` |
| external_ids | object[] | 1:N, Optional | Cross-references to PubChem, HMDB, ChEBI | `[{"database": "PubChem", "id": "5280343"}, {"database": "HMDB", "id": "HMDB0005794"}]` |

### USDA FoodData Central

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| fdc_id | integer | 1:1, Required | FoodData Central ID | `167512` |
| data_type | string | 1:1, Optional | Data type classification | `Foundation`, `SR Legacy`, `Branded`, `FNDDS` |
| nutrient_id | integer | 1:1, Optional | Nutrient identifier | `1162` |
| nutrient_nbr | string | 1:1, Optional | USDA nutrient number | `401` |
| food_portions | object[] | 1:N, Optional | Serving sizes with gram weights | `[{"measure_unit": "cup", "gram_weight": 160, "amount": 1}]` |
| gtin_upc | string | 1:1, Optional | Product barcode (branded foods) | `0018787001007` |

---

## Metadata Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| _source.database | string | 1:1, Required | Source database name | `Phenol-Explorer`, `PhytoHub`, `USDA FoodData` |
| _source.version | string | 1:1, Optional | Database version | `3.6`, `1.0`, `2024-04` |
| _source.access_date | date | 1:1, Optional | Date data was accessed | `2026-01-24` |
| _source.original_id | string | 1:1, Optional | Original identifier from source | - |

---

## Field Mappings by Source

### Phenol-Explorer to Unified Schema

| Phenol-Explorer Field | Unified Field |
|-----------------------|---------------|
| compound_id | compound_id |
| name | name |
| smiles | smiles |
| molecular_weight | molecular_weight |
| class | compound_class |
| subclass | subclass |
| pubchem_cid | pubchem_cid |
| cas_number | cas_number |
| metabolites | metabolites |
| pharmacokinetics | pharmacokinetics |
| food_content | food_sources |

### PhytoHub to Unified Schema

| PhytoHub Field | Unified Field |
|----------------|---------------|
| phytohub_id | phytohub_id |
| name | name |
| smiles | smiles |
| inchi_key | inchi_key |
| molecular_weight | molecular_weight |
| compound_type | compound_type |
| parent_id | parent_id |
| transformation | transformation |
| ms_spectra | ms_spectra |
| external_ids | external_ids |
| food_content | food_sources |

### USDA FoodData Central to Unified Schema

| USDA Field | Unified Field |
|------------|---------------|
| fdcId | fdc_id |
| description | name |
| dataType | data_type |
| nutrient_id | nutrient_id |
| nutrientNumber | nutrient_nbr |
| foodPortions | food_portions |
| gtinUpc | gtin_upc |
| amount | content_value |
| unitName | content_unit |

---

## Data Quality Notes

- **Required Field:** Only `name` is strictly required across all sources
- **Structure Data:** Chemical structure fields primarily from Phenol-Explorer and PhytoHub (phytochemicals); USDA focuses on nutrients without structures
- **Food Content:** All three sources provide food-compound associations with quantitative values
- **Metabolites:** Phenol-Explorer and PhytoHub provide metabolite data; USDA does not track metabolites
- **Pharmacokinetics:** Only Phenol-Explorer provides absorption kinetics data
- **MS Spectra:** Only PhytoHub provides mass spectrometry reference spectra

### Polyphenol Classes (Phenol-Explorer)

| Class | Description |
|-------|-------------|
| Flavonoids | Largest class; includes flavonols, flavones, isoflavones, flavanones, anthocyanins, flavanols |
| Phenolic acids | Hydroxybenzoic and hydroxycinnamic acids |
| Stilbenes | Includes resveratrol and derivatives |
| Lignans | Plant-derived polyphenols metabolized by gut microbiota |
| Other polyphenols | Tyrosols, alkylphenols, curcuminoids, etc. |

### USDA Data Types

| Data Type | Description |
|-----------|-------------|
| Foundation | Nutrient and food component values from analytical data |
| SR Legacy | Standard Reference data (historical) |
| Branded | Branded food products with label data |
| FNDDS | Food and Nutrient Database for Dietary Studies |
| Survey (FNDDS) | Foods reported in dietary surveys |
| Experimental | Research-based data |

### Metabolite Types (PhytoHub)

| Type | Description |
|------|-------------|
| parent | Original dietary compound |
| phase1 | Phase I metabolism products (oxidation, reduction, hydrolysis) |
| phase2 | Phase II metabolism products (conjugation: glucuronidation, sulfation, methylation) |
| microbial | Gut microbiota metabolites |

### Common Content Units

| Unit | Description |
|------|-------------|
| mg/100g | Milligrams per 100 grams (solid foods) |
| mg/100ml | Milligrams per 100 milliliters (beverages) |
| ug/100g | Micrograms per 100 grams |
| mg/serving | Milligrams per standard serving |
| mg/L | Milligrams per liter |
