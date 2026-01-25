# Phenol-Explorer - Data Dictionary

## Overview

This data dictionary documents the schema for Phenol-Explorer polyphenol database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | phenol.explorer |
| **Name** | Phenol-Explorer |
| **Parent** | 2.4.food.compounds.nutrients |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| compound_id | integer | 1:1 | Yes | Primary identifier | 1 |
| name | string | 1:1 | Yes | Compound name | Quercetin |
| class | string | 1:1 | Yes | Polyphenol class | Flavonoids |
| subclass | string | 1:1 | No | Subclass | Flavonols |
| pubchem_cid | integer | 1:1 | No | PubChem compound ID | 5280343 |
| cas_number | string | 1:1 | No | CAS registry number | 117-39-5 |
| molecular_formula | string | 1:1 | No | Chemical formula | C15H10O7 |
| molecular_weight | decimal | 1:1 | No | MW in Daltons | 302.24 |
| smiles | string | 1:1 | No | SMILES structure | O=c1c(O)c(-c2ccc(O)c(O)c2)... |

### Food Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| food_id | integer | 1:1 | Yes | Primary identifier | 101 |
| food_name | string | 1:1 | Yes | Food name | Onion, red, raw |
| food_group | string | 1:1 | No | Food category | Vegetables |
| scientific_name | string | 1:1 | No | Plant species | Allium cepa |
| langual_code | string | 1:1 | No | LanguaL food code | A0148 |

### Food Content Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| content_id | integer | 1:1 | Yes | Primary identifier | 1001 |
| content_value | decimal | 1:1 | Yes | Amount in food | 39.21 |
| content_unit | string | 1:1 | Yes | Unit of measurement | mg/100g |
| min_value | decimal | 1:1 | No | Minimum reported | 20.3 |
| max_value | decimal | 1:1 | No | Maximum reported | 50.6 |
| data_points | integer | 1:1 | No | Number of samples | 12 |
| analytical_method | string | 1:1 | No | Analysis method | HPLC-DAD |
| reference_id | integer | 1:1 | No | Literature source | 101 |

### Metabolite Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| metabolite_id | integer | 1:1 | Yes | Primary identifier | 201 |
| metabolite_name | string | 1:1 | Yes | Metabolite name | Quercetin-3-O-glucuronide |
| metabolite_type | string | 1:1 | Yes | Metabolism phase | Phase II |
| transformation | string | 1:1 | No | Transformation type | Glucuronidation |
| biofluid | string | 1:1 | No | Detection matrix | plasma, urine |

### Pharmacokinetics

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| pk_id | integer | 1:1 | Yes | Primary identifier | 301 |
| dose | decimal | 1:1 | No | Administered dose | 150 |
| dose_unit | string | 1:1 | No | Dose unit | mg |
| food_matrix | string | 1:1 | No | Source food/pure | Onion |
| cmax | decimal | 1:1 | No | Maximum concentration | 1.5 |
| cmax_unit | string | 1:1 | No | Concentration unit | umol/L |
| tmax | decimal | 1:1 | No | Time to Cmax | 0.7 |
| tmax_unit | string | 1:1 | No | Time unit | hours |
| auc | decimal | 1:1 | No | Area under curve | 3.2 |
| half_life | decimal | 1:1 | No | Elimination half-life | 3.8 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Compound ID | Integer | 1 | Internal identifier |
| PubChem CID | Integer | 5280343 | Chemical database |
| CAS Number | XX-XX-X | 117-39-5 | Chemical registry |
| LanguaL Code | A + digits | A0148 | Food descriptor |

---

## Enumerations

### Polyphenol Classes

| Class | Subclasses |
|-------|------------|
| Flavonoids | Flavonols, Flavones, Flavan-3-ols, Anthocyanins, Flavanones, Isoflavones |
| Phenolic acids | Hydroxybenzoic acids, Hydroxycinnamic acids |
| Stilbenes | Resveratrol and derivatives |
| Lignans | Secoisolariciresinol, Matairesinol |
| Other polyphenols | Curcuminoids, Ellagitannins, Tyrosols |

### Flavonoid Subclasses

| Subclass | Examples |
|----------|----------|
| Flavonols | Quercetin, Kaempferol, Myricetin |
| Flavones | Luteolin, Apigenin |
| Flavan-3-ols | Catechin, Epicatechin, EGCG |
| Anthocyanins | Cyanidin, Delphinidin, Malvidin |
| Flavanones | Naringenin, Hesperetin |
| Isoflavones | Genistein, Daidzein |

### Metabolite Types

| Type | Description | Examples |
|------|-------------|----------|
| Phase I | Oxidation/reduction | Hydroxylated metabolites |
| Phase II | Conjugation | Glucuronides, sulfates |
| Microbial | Gut bacteria | Ring fission products |

### Food Groups

| Group | Examples |
|-------|----------|
| Fruits | Apples, berries, citrus |
| Vegetables | Onions, broccoli, spinach |
| Beverages | Tea, coffee, wine |
| Grains | Whole wheat, oats |
| Legumes | Soybeans, lentils |
| Nuts | Walnuts, pecans |

### Analytical Methods

| Method | Description |
|--------|-------------|
| HPLC-DAD | High-performance liquid chromatography with diode array |
| HPLC-MS | LC with mass spectrometry |
| HPLC-UV | LC with UV detection |
| GC-MS | Gas chromatography-MS |

---

## Entity Relationships

### Compound to Food Contents
- **Cardinality:** 1:N
- **Description:** One compound measured in multiple foods
- **Key Fields:** compound_id, food_id

### Food to Contents
- **Cardinality:** 1:N
- **Description:** One food contains multiple polyphenols
- **Key Fields:** food_id, compound_id

### Compound to Metabolites
- **Cardinality:** 1:N
- **Description:** One parent compound produces multiple metabolites
- **Key Fields:** compound_id, metabolite_id

### Compound to Pharmacokinetics
- **Cardinality:** 1:N
- **Description:** Multiple PK studies per compound
- **Key Fields:** compound_id, pk_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HPLC | High-Performance Liquid Chromatography | Analysis method |
| DAD | Diode Array Detector | HPLC detector |
| MS | Mass Spectrometry | Detection method |
| GC | Gas Chromatography | Analysis method |
| PK | Pharmacokinetics | Bioavailability |
| Cmax | Maximum Concentration | PK parameter |
| Tmax | Time to Maximum | PK parameter |
| AUC | Area Under the Curve | PK parameter |
| UGT | UDP-glucuronosyltransferase | Phase II enzyme |
| SULT | Sulfotransferase | Phase II enzyme |
| COMT | Catechol-O-methyltransferase | Phase II enzyme |
| EGCG | Epigallocatechin Gallate | Tea polyphenol |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubChem | CID | Chemical data |
| ChEBI | ChEBI ID | Chemical ontology |
| HMDB | HMDB ID | Metabolite data |
| FooDB | FooDB ID | Food composition |
| PhytoHub | PhytoHub ID | Metabolomics |

---

## Data Quality Notes

1. **Comprehensive Coverage:** 500+ polyphenol compounds in 450+ foods
2. **Quantitative Data:** 35,000+ content data points with ranges
3. **Metabolite Focus:** 380+ metabolites for bioavailability studies
4. **Analytical Methods:** Data includes method information
5. **Literature Linked:** 1,500+ publications referenced
6. **Standard Units:** Content in mg/100g fresh weight
7. **Open Access:** Free for academic and research use
