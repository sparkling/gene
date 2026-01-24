---
id: schema-ebasis
title: "eBASIS - Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, bioactive, phytochemicals, food, composition]
---

# eBASIS - Bioactive Substances in Food Information System Schema Documentation

**Document ID:** SCHEMA-EBASIS
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://ebasis.eurofir.org/

---

## TL;DR

eBASIS (Bioactive Substances in Food Information System) is a EuroFIR database containing quality-evaluated composition data on bioactive compounds in plant-based foods. It focuses on compounds with potential health benefits beyond basic nutrition, including polyphenols, carotenoids, phytosterols, and glucosinolates.

---

## Database Statistics

| Metric | Count |
|--------|-------|
| Bioactive Compounds | 300+ |
| Foods | 500+ |
| Data Points | 25,000+ |
| Component Groups | 15+ |
| Quality Evaluated | Yes (EuroFIR criteria) |

---

## Database Schema

### Core Tables

#### 1. Components (Bioactive Compounds)
```sql
CREATE TABLE components (
    component_id INTEGER PRIMARY KEY,
    component_name VARCHAR(255),      -- Compound name
    component_group VARCHAR(100),     -- Category (e.g., Flavonoids)
    component_subgroup VARCHAR(100),  -- Subcategory (e.g., Flavonols)
    cas_number VARCHAR(20),           -- CAS registry number
    pubchem_cid INTEGER,              -- PubChem compound ID
    molecular_formula VARCHAR(50),
    molecular_weight DECIMAL(10,4),
    iupac_name VARCHAR(500),
    inchi VARCHAR(1000),
    inchikey VARCHAR(50),
    synonyms TEXT,                    -- Alternative names
    INDEX idx_group (component_group),
    INDEX idx_cas (cas_number)
);
```

#### 2. Foods
```sql
CREATE TABLE foods (
    food_id INTEGER PRIMARY KEY,
    food_name VARCHAR(255),           -- Common name
    scientific_name VARCHAR(255),     -- Latin binomial
    food_group VARCHAR(100),          -- Category (e.g., Fruits)
    langual_code VARCHAR(20),         -- LanguaL classification
    eurofir_code VARCHAR(20),         -- EuroFIR food code
    edible_portion DECIMAL(5,2),      -- Percent edible
    INDEX idx_group (food_group)
);
```

#### 3. Composition Data
```sql
CREATE TABLE composition (
    id INTEGER PRIMARY KEY,
    food_id INTEGER REFERENCES foods(food_id),
    component_id INTEGER REFERENCES components(component_id),
    value DECIMAL(15,6),              -- Concentration value
    unit VARCHAR(20),                 -- Unit (mg/100g, mcg/100g)
    value_type VARCHAR(20),           -- "mean", "median", "single"
    min_value DECIMAL(15,6),          -- Minimum observed
    max_value DECIMAL(15,6),          -- Maximum observed
    n_samples INTEGER,                -- Number of samples
    std_dev DECIMAL(15,6),            -- Standard deviation
    data_quality_score INTEGER,       -- EuroFIR quality rating (1-5)
    analytical_method VARCHAR(100),   -- HPLC, GC-MS, etc.
    reference_id INTEGER,             -- Source reference
    INDEX idx_food (food_id),
    INDEX idx_component (component_id),
    INDEX idx_quality (data_quality_score)
);
```

#### 4. References
```sql
CREATE TABLE references (
    reference_id INTEGER PRIMARY KEY,
    citation TEXT,                    -- Full citation
    pubmed_id INTEGER,                -- PubMed ID if available
    doi VARCHAR(100),                 -- DOI if available
    year INTEGER,
    authors TEXT,
    journal VARCHAR(255),
    title TEXT
);
```

---

## Component Groups (Bioactive Categories)

| Group | Subgroups | Examples |
|-------|-----------|----------|
| Flavonoids | Flavonols, Flavones, Flavanones, Isoflavones, Anthocyanidins, Flavanols | Quercetin, Kaempferol |
| Carotenoids | Carotenes, Xanthophylls | Beta-carotene, Lutein |
| Phenolic Acids | Hydroxybenzoic, Hydroxycinnamic | Gallic acid, Caffeic acid |
| Stilbenes | - | Resveratrol |
| Lignans | - | Secoisolariciresinol |
| Phytosterols | Sterols, Stanols | Beta-sitosterol |
| Glucosinolates | Aliphatic, Indole, Aromatic | Glucoraphanin, Sinigrin |
| Terpenes | Monoterpenes, Diterpenes | Limonene |
| Organosulfur | - | Allicin, Alliin |
| Alkaloids | - | Caffeine, Theobromine |

---

## Data Quality Framework (EuroFIR)

### Quality Score Scale

| Score | Description | Criteria |
|-------|-------------|----------|
| 5 | Excellent | Full documentation, validated method |
| 4 | Good | Good documentation, appropriate method |
| 3 | Acceptable | Adequate documentation |
| 2 | Poor | Limited documentation |
| 1 | Very Poor | Minimal information |

### Quality Evaluation Criteria

| Criterion | Weight | Assessment |
|-----------|--------|------------|
| Food description | High | LanguaL coding, botanical ID |
| Sampling procedure | High | Representative, documented |
| Analytical method | High | Validated, appropriate |
| Analytical quality control | Medium | Standards, replicates |
| Number of samples | Medium | Statistical adequacy |
| Component identification | High | Confirmed identity |

---

## JSON Schemas

### Component Object
```json
{
  "component_id": 1234,
  "component_name": "Quercetin",
  "component_group": "Flavonoids",
  "component_subgroup": "Flavonols",
  "cas_number": "117-39-5",
  "pubchem_cid": 5280343,
  "molecular_formula": "C15H10O7",
  "molecular_weight": 302.236,
  "iupac_name": "2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxychromen-4-one",
  "inchikey": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
  "synonyms": ["3,3',4',5,7-Pentahydroxyflavone", "Sophoretin"]
}
```

### Composition Record Object
```json
{
  "id": 5678,
  "food": {
    "food_id": 101,
    "food_name": "Apple, raw, with skin",
    "scientific_name": "Malus domestica"
  },
  "component": {
    "component_id": 1234,
    "component_name": "Quercetin",
    "component_group": "Flavonoids"
  },
  "value": 4.01,
  "unit": "mg/100g",
  "value_type": "mean",
  "min_value": 2.10,
  "max_value": 6.45,
  "n_samples": 12,
  "std_dev": 1.23,
  "data_quality_score": 4,
  "analytical_method": "HPLC-DAD",
  "reference": {
    "pubmed_id": 12345678,
    "citation": "Smith et al. (2020)"
  }
}
```

### Food Object
```json
{
  "food_id": 101,
  "food_name": "Apple, raw, with skin",
  "scientific_name": "Malus domestica",
  "food_group": "Fruits",
  "langual_code": "A0148",
  "eurofir_code": "FR001234",
  "edible_portion": 95.0,
  "components": [
    {
      "component_name": "Quercetin",
      "value": 4.01,
      "unit": "mg/100g"
    },
    {
      "component_name": "Catechin",
      "value": 8.50,
      "unit": "mg/100g"
    }
  ]
}
```

---

## Data Access

### Web Interface

| Feature | Access |
|---------|--------|
| Browse foods | Registration required |
| Browse components | Registration required |
| Query composition | Registration required |
| Download data | Contact EuroFIR |

### API Access

- **Status:** Not publicly available
- **Contact:** EuroFIR for data access agreements

### Bulk Downloads

- **Status:** Restricted
- **Contact:** info@eurofir.org
- **License:** EuroFIR data use agreement required

---

## Relationships

### Entity Relationship Diagram
```
foods (1) ----< composition >---- (N) components
                    |
                    v
               references
```

### Foreign Key Relationships

| Table | Column | References |
|-------|--------|------------|
| composition | food_id | foods.food_id |
| composition | component_id | components.component_id |
| composition | reference_id | references.reference_id |

---

## Use Cases

### 1. Find Bioactive Content in a Food
```
Query: What are the flavonoid contents in apples?
Filter: food_name = "Apple", component_group = "Flavonoids"
Result: Quercetin 4.01 mg/100g, Catechin 8.50 mg/100g, etc.
```

### 2. Find Food Sources of a Compound
```
Query: Which foods contain quercetin?
Filter: component_name = "Quercetin"
Sort by: value DESC
Result: Onions (20.3 mg/100g), Apples (4.01 mg/100g), etc.
```

### 3. Dietary Intake Assessment
```
Query: Calculate total polyphenol intake from diet
Method: Sum (food_amount * component_value) for all foods
Use quality scores to weight reliability
```

---

## License

- **License Type:** EuroFIR Data Use Agreement
- **Academic Use:** Free with agreement
- **Commercial Use:** Requires license
- **Attribution:** Required
- **Contact:** info@eurofir.org

---

## Citation

```
eBASIS - Bioactive Substances in Food Information System.
EuroFIR AISBL.
https://ebasis.eurofir.org/
[Access date]

Kiely M, Black LJ, Plumb J, et al. (2010).
EuroFIR eBASIS: Application for health claims submissions.
Nutrition Bulletin, 35(2), 170-174.
```

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Bioactive compound | Non-nutrient food component with biological activity | Quercetin |
| Component group | Major category of bioactive compounds | Flavonoids |
| Data quality score | EuroFIR rating of data reliability (1-5) | 4 |
| LanguaL code | Food description thesaurus code | A0148 |
| mg/100g | Milligrams per 100 grams edible portion | 4.01 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Polyphenol | Plant compound with multiple phenol groups | Flavonoids, phenolic acids |
| Flavonoid | Class of plant polyphenols with flavone backbone | Component group |
| Carotenoid | Terpenoid pigments in plants | Component group |
| HPLC-DAD | High-performance liquid chromatography with diode array detection | Analytical method |
| Edible portion | Percentage of food that is consumable | Food description |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| eBASIS | Bioactive Substances in Food Information System | Database name |
| EuroFIR | European Food Information Resource | Network |
| CAS | Chemical Abstracts Service | Identifier system |
| LanguaL | Langua aLimentaria | Food thesaurus |
| HPLC | High-Performance Liquid Chromatography | Analytical method |
| DAD | Diode Array Detector | Detection method |
| GC-MS | Gas Chromatography-Mass Spectrometry | Analytical method |

---

## Sample Data

### Example Composition Record (Full Detail)
```json
{
  "id": 12345,
  "food": {
    "food_id": 101,
    "food_name": "Apple, raw, with skin",
    "scientific_name": "Malus domestica",
    "cultivar": "Gala",
    "food_group": "Fruits",
    "food_subgroup": "Pome fruits",
    "langual_codes": ["A0148", "B1001", "C0131"],
    "eurofir_code": "FR001234",
    "edible_portion_percent": 95.0,
    "preparation": "Raw, unprocessed"
  },
  "component": {
    "component_id": 1234,
    "component_name": "Quercetin",
    "component_group": "Flavonoids",
    "component_subgroup": "Flavonols",
    "cas_number": "117-39-5",
    "pubchem_cid": 5280343,
    "molecular_formula": "C15H10O7",
    "molecular_weight": 302.236,
    "inchikey": "REFJWTPEDVJJIY-UHFFFAOYSA-N"
  },
  "composition_data": {
    "value": 4.01,
    "unit": "mg/100g",
    "basis": "Fresh weight, edible portion",
    "value_type": "mean",
    "min_value": 2.10,
    "max_value": 6.45,
    "n_samples": 12,
    "n_data_points": 36,
    "std_dev": 1.23,
    "cv_percent": 30.7,
    "confidence_interval_95": [3.42, 4.60]
  },
  "quality_assessment": {
    "overall_score": 4,
    "food_description_score": 5,
    "sampling_plan_score": 4,
    "sample_handling_score": 4,
    "analytical_method_score": 4,
    "analytical_qc_score": 3,
    "component_id_score": 5
  },
  "analytical_details": {
    "method": "HPLC-DAD",
    "method_reference": "AOAC 2005.02",
    "detection_limit": 0.1,
    "quantification_limit": 0.5,
    "extraction_solvent": "Methanol/water 70:30",
    "hydrolysis": "None (free form measured)"
  },
  "reference": {
    "reference_id": 5678,
    "authors": "Smith J, Jones M, Williams K",
    "year": 2020,
    "title": "Flavonoid content of apple cultivars",
    "journal": "Journal of Food Composition and Analysis",
    "volume": 85,
    "pages": "103567",
    "doi": "10.1016/j.jfca.2020.103567",
    "pubmed_id": 32145678
  },
  "metadata": {
    "date_entered": "2020-05-15",
    "last_reviewed": "2024-01-10",
    "reviewed_by": "EuroFIR Quality Team"
  }
}
```

### Example Component Record (Full Detail)
```json
{
  "component_id": 1234,
  "component_name": "Quercetin",
  "component_group": "Flavonoids",
  "component_subgroup": "Flavonols",
  "systematic_name": "3,3',4',5,7-Pentahydroxyflavone",
  "iupac_name": "2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxychromen-4-one",
  "synonyms": [
    "Sophoretin",
    "Meletin",
    "Quercetine",
    "Xanthaurine"
  ],
  "identifiers": {
    "cas_number": "117-39-5",
    "pubchem_cid": 5280343,
    "chebi_id": "CHEBI:16243",
    "kegg_id": "C00389",
    "hmdb_id": "HMDB0005794"
  },
  "chemical_properties": {
    "molecular_formula": "C15H10O7",
    "molecular_weight": 302.236,
    "exact_mass": 302.0427,
    "inchi": "InChI=1S/C15H10O7/c16-7-4-10(19)12-11(5-7)22-15(14(21)13(12)20)6-1-2-8(17)9(18)3-6/h1-5,16-19,21H",
    "inchikey": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
    "smiles": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
    "logp": 1.54
  },
  "biological_activity": {
    "antioxidant": true,
    "anti_inflammatory": true,
    "orac_value_umol_te_umol": 4.7
  },
  "dietary_sources": [
    {
      "food": "Capers, raw",
      "mean_content_mg_100g": 233.84
    },
    {
      "food": "Red onion, raw",
      "mean_content_mg_100g": 39.21
    },
    {
      "food": "Cranberries, raw",
      "mean_content_mg_100g": 14.02
    },
    {
      "food": "Apple, raw with skin",
      "mean_content_mg_100g": 4.01
    }
  ],
  "related_compounds": [
    {
      "name": "Quercetin-3-glucoside",
      "relationship": "Glycoside"
    },
    {
      "name": "Isorhamnetin",
      "relationship": "Methylated derivative"
    },
    {
      "name": "Rutin",
      "relationship": "Rutinoside"
    }
  ]
}
```

### Example Food Record (Full Detail)
```json
{
  "food_id": 101,
  "food_name": "Apple, raw, with skin",
  "scientific_name": "Malus domestica",
  "family": "Rosaceae",
  "common_names": {
    "en": "Apple",
    "de": "Apfel",
    "fr": "Pomme",
    "es": "Manzana"
  },
  "classification": {
    "food_group": "Fruits",
    "food_subgroup": "Pome fruits",
    "langual_codes": {
      "product_type": "A0148",
      "food_source": "B1001",
      "part_of_plant": "C0131",
      "physical_state": "E0101",
      "treatment": "H0100"
    },
    "eurofir_code": "FR001234"
  },
  "nutritional_info": {
    "water_g_100g": 85.56,
    "energy_kcal_100g": 52,
    "edible_portion_percent": 95.0
  },
  "bioactive_summary": {
    "total_flavonoids_mg_100g": 22.3,
    "total_phenolics_mg_gae_100g": 296.3,
    "component_groups_present": [
      "Flavonols",
      "Flavanols",
      "Hydroxycinnamic acids",
      "Dihydrochalcones"
    ],
    "major_compounds": [
      {"name": "Epicatechin", "mg_100g": 7.20},
      {"name": "Procyanidin B2", "mg_100g": 6.10},
      {"name": "Quercetin", "mg_100g": 4.01},
      {"name": "Chlorogenic acid", "mg_100g": 3.50}
    ]
  },
  "data_quality": {
    "n_composition_records": 45,
    "avg_quality_score": 3.8,
    "last_data_update": "2024-06-15"
  }
}
```

**Note:** The sample data above represents the structure of data available through eBASIS. Actual data requires registration and a EuroFIR data use agreement. Values shown are illustrative and based on publicly available information about polyphenol content in foods.

---

## Related Documents

- [Download Instructions](./download.md)
- [FooDB](../../6.1.food.composition/foodb/_index.md) - Food compound database
- [Phenol-Explorer](../phenol.explorer/) - Polyphenol database
