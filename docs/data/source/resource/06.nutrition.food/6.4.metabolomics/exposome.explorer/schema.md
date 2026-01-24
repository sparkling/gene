---
id: schema-exposome-explorer
title: "Exposome-Explorer - Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, exposome, biomarkers, diet, environment, metabolomics]
---

# Exposome-Explorer - Schema Documentation

**Document ID:** SCHEMA-EXPOSOME-EXPLORER
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** http://exposome-explorer.iarc.fr/

---

## TL;DR

Exposome-Explorer is a manually curated database of biomarkers of dietary and environmental exposures. It systematically catalogs published associations between biomarkers measured in human biological samples and specific exposures including foods, nutrients, pollutants, and other environmental factors.

---

## Database Statistics

| Metric | Count |
|--------|-------|
| Biomarkers | 900+ |
| Dietary Exposures | 600+ |
| Publications | 1,500+ |
| Specimen Types | 10+ |
| Exposure Categories | 20+ |

---

## Database Schema

### Core Tables

#### 1. Biomarkers Table
```sql
CREATE TABLE biomarkers (
    biomarker_id INTEGER PRIMARY KEY,
    biomarker_name VARCHAR(255),      -- Common name
    pubchem_cid INTEGER,              -- PubChem compound ID
    hmdb_id VARCHAR(20),              -- HMDB identifier
    chebi_id VARCHAR(20),             -- ChEBI identifier
    cas_number VARCHAR(20),           -- CAS registry number
    molecular_formula VARCHAR(50),
    biomarker_type VARCHAR(50),       -- "metabolite", "nutrient", etc.
    description TEXT,
    INDEX idx_pubchem (pubchem_cid),
    INDEX idx_hmdb (hmdb_id)
);
```

#### 2. Exposures Table
```sql
CREATE TABLE exposures (
    exposure_id INTEGER PRIMARY KEY,
    exposure_name VARCHAR(255),       -- Exposure name
    exposure_category VARCHAR(100),   -- "Dietary", "Environmental", etc.
    exposure_subcategory VARCHAR(100),-- "Fruits", "Vegetables", etc.
    food_group VARCHAR(100),          -- If dietary
    description TEXT,
    INDEX idx_category (exposure_category)
);
```

#### 3. Associations Table
```sql
CREATE TABLE associations (
    association_id INTEGER PRIMARY KEY,
    biomarker_id INTEGER REFERENCES biomarkers(biomarker_id),
    exposure_id INTEGER REFERENCES exposures(exposure_id),
    specimen_type VARCHAR(50),        -- "plasma", "urine", "serum", etc.
    association_type VARCHAR(50),     -- "positive", "negative", "correlation"
    effect_size DECIMAL(10,4),        -- Correlation coefficient or OR
    effect_unit VARCHAR(50),          -- Unit or metric type
    p_value DECIMAL(20,15),           -- Statistical significance
    confidence_interval VARCHAR(50),   -- 95% CI
    sample_size INTEGER,              -- Study sample size
    population VARCHAR(100),          -- Study population
    reference_id INTEGER,             -- Publication reference
    INDEX idx_biomarker (biomarker_id),
    INDEX idx_exposure (exposure_id),
    INDEX idx_specimen (specimen_type)
);
```

#### 4. Concentrations Table
```sql
CREATE TABLE concentrations (
    concentration_id INTEGER PRIMARY KEY,
    biomarker_id INTEGER REFERENCES biomarkers(biomarker_id),
    specimen_type VARCHAR(50),        -- "plasma", "urine", "serum"
    mean_concentration DECIMAL(15,6),
    unit VARCHAR(30),                 -- "umol/L", "ng/mL", etc.
    min_concentration DECIMAL(15,6),
    max_concentration DECIMAL(15,6),
    std_dev DECIMAL(15,6),
    n_subjects INTEGER,
    population VARCHAR(100),
    reference_id INTEGER
);
```

#### 5. References Table
```sql
CREATE TABLE references (
    reference_id INTEGER PRIMARY KEY,
    pubmed_id INTEGER,
    doi VARCHAR(100),
    authors TEXT,
    title TEXT,
    journal VARCHAR(255),
    year INTEGER,
    study_type VARCHAR(50),           -- "cross-sectional", "cohort", etc.
    study_design VARCHAR(100)
);
```

---

## Exposure Categories

### Dietary Exposures

| Category | Examples |
|----------|----------|
| Fruits | Citrus fruits, berries, apples |
| Vegetables | Cruciferous, leafy greens, allium |
| Beverages | Coffee, tea, wine, beer |
| Animal products | Red meat, fish, dairy |
| Grains | Whole grains, refined grains |
| Legumes | Soy, beans, lentils |
| Nuts and seeds | Tree nuts, peanuts |
| Fats and oils | Olive oil, vegetable oils |
| Sweeteners | Sugar, artificial sweeteners |
| Additives | Preservatives, colorants |

### Environmental Exposures

| Category | Examples |
|----------|----------|
| Pollutants | PAHs, PCBs, dioxins |
| Heavy metals | Lead, mercury, arsenic, cadmium |
| Pesticides | Organophosphates, pyrethroids |
| Tobacco | Active smoking, passive smoke |
| Occupational | Industrial chemicals, solvents |

---

## Specimen Types

| Type | Description | Common Biomarkers |
|------|-------------|-------------------|
| Plasma | Blood plasma | Fatty acids, carotenoids |
| Serum | Blood serum | Vitamins, minerals |
| Urine | Spot or 24h collection | Polyphenol metabolites |
| Whole blood | Complete blood | Heavy metals |
| Erythrocytes | Red blood cells | Long-term markers |
| Adipose tissue | Fat tissue biopsy | Lipophilic compounds |
| Hair | Head hair | Chronic exposures |
| Nails | Fingernails/toenails | Selenium, arsenic |

---

## JSON Schemas

### Biomarker Object
```json
{
  "biomarker_id": 123,
  "biomarker_name": "4-O-methylgallic acid",
  "pubchem_cid": 13159,
  "hmdb_id": "HMDB0002122",
  "cas_number": "3934-84-7",
  "molecular_formula": "C8H8O5",
  "biomarker_type": "Phenolic acid metabolite",
  "description": "Urinary metabolite of gallic acid from dietary polyphenols"
}
```

### Association Object
```json
{
  "association_id": 456,
  "biomarker": {
    "biomarker_id": 123,
    "biomarker_name": "4-O-methylgallic acid"
  },
  "exposure": {
    "exposure_id": 789,
    "exposure_name": "Tea consumption",
    "exposure_category": "Dietary"
  },
  "specimen_type": "urine",
  "association_type": "positive",
  "effect_size": 0.42,
  "effect_unit": "Pearson correlation",
  "p_value": 0.00001,
  "sample_size": 450,
  "population": "European adults",
  "reference": {
    "pubmed_id": 23456789,
    "year": 2019
  }
}
```

### Concentration Object
```json
{
  "biomarker": "Quercetin",
  "specimen_type": "plasma",
  "mean_concentration": 0.065,
  "unit": "umol/L",
  "min_concentration": 0.01,
  "max_concentration": 0.35,
  "n_subjects": 238,
  "population": "Free-living adults, US"
}
```

---

## Data Access

### Web Interface

| Feature | URL |
|---------|-----|
| Browse biomarkers | http://exposome-explorer.iarc.fr/biomarkers |
| Browse exposures | http://exposome-explorer.iarc.fr/exposures |
| Search associations | http://exposome-explorer.iarc.fr/search |

### Downloads

| Format | URL |
|--------|-----|
| TSV exports | http://exposome-explorer.iarc.fr/downloads |

### API Access

- **Status:** No formal API
- **Alternative:** Web scraping (with rate limits)

---

## Relationships

### Entity Relationship Diagram
```
biomarkers (N) ----< associations >---- (M) exposures
     |                    |
     v                    v
concentrations       references
```

### Foreign Key Relationships

| Table | Column | References |
|-------|--------|------------|
| associations | biomarker_id | biomarkers.biomarker_id |
| associations | exposure_id | exposures.exposure_id |
| associations | reference_id | references.reference_id |
| concentrations | biomarker_id | biomarkers.biomarker_id |

---

## Use Cases

### 1. Find Biomarkers for a Food
```
Query: What are validated biomarkers of coffee consumption?
Filter: exposure_name LIKE '%coffee%'
Result: Trigonelline, caffeine metabolites, etc.
```

### 2. Find Dietary Sources of a Biomarker
```
Query: What foods are associated with urinary hippuric acid?
Filter: biomarker_name = 'hippuric acid'
Result: Tea, wine, fruits, whole grains
```

### 3. Assess Biomarker Reliability
```
Query: Evaluate plasma carotenoids as fruit/vegetable biomarkers
Check: Effect sizes, consistency across studies, specimen types
```

---

## Study Design Information

| Study Type | Description |
|------------|-------------|
| Cross-sectional | Single time point measurement |
| Cohort | Longitudinal follow-up |
| Intervention | Controlled feeding study |
| Case-control | Cases vs. controls comparison |

---

## License

- **License Type:** Creative Commons Attribution 4.0 (CC BY 4.0)
- **Commercial Use:** Yes
- **Attribution:** Required
- **Redistribution:** Permitted with attribution

---

## Citation

```
Neveu V, Moussy A, Rber H, et al.
Exposome-Explorer: a manually-curated database on biomarkers of exposure to dietary and environmental factors.
Nucleic Acids Research. 2017;45(D1):D979-D984.
doi: 10.1093/nar/gkw980

Exposome-Explorer. IARC, WHO.
http://exposome-explorer.iarc.fr/
```

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Biomarker | Measurable indicator of exposure | Urinary hippuric acid |
| Exposure | Dietary or environmental factor | Coffee consumption |
| Specimen | Biological sample type | Plasma, urine |
| Association | Statistical relationship | Correlation r=0.42 |
| Effect size | Magnitude of association | Pearson r, odds ratio |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Exposome | Total environmental exposures over lifetime | Study framework |
| Dietary biomarker | Indicator of food intake | Exposure assessment |
| Metabolite | Breakdown product of metabolism | Biomarker type |
| Controlled feeding | Study with provided diet | Gold standard validation |
| Free-living | Normal dietary conditions | Study population |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IARC | International Agency for Research on Cancer | Database host |
| WHO | World Health Organization | Parent organization |
| HMDB | Human Metabolome Database | Cross-reference |
| PAH | Polycyclic Aromatic Hydrocarbons | Environmental pollutants |
| PCB | Polychlorinated Biphenyls | Environmental pollutants |
| OR | Odds Ratio | Effect size measure |
| CI | Confidence Interval | Statistical uncertainty |

---

## Related Documents

- [Download Instructions](./download.md)
- [HMDB](../hmdb/_index.md) - Human Metabolome Database
- [FooDB](../../6.1.food.composition/foodb/_index.md) - Food compound database
