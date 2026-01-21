# FooDB - Schema Documentation

**Document ID:** SCHEMA-FOODB
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://foodb.ca/

---

## Overview

FooDB is the world's largest and most comprehensive resource on food constituents, chemistry, and biology. It provides information on both macronutrients and micronutrients, including many of the constituents that give foods their flavor, color, taste, texture, and aroma.

---

## Database Statistics

| Metric | Count |
|--------|-------|
| Foods | 778 |
| Compounds | 70,926 |
| Compound-Food Associations | Multiple per food |
| Total Pages (Compounds) | 2,838 pages (25/page) |
| Total Pages (Foods) | 32 pages (25/page) |

---

## Database Schema

### Core Tables

#### 1. Foods Table
```sql
CREATE TABLE foods (
    food_id INTEGER PRIMARY KEY,
    public_id VARCHAR(20),        -- e.g., "FOOD00001"
    name VARCHAR(255),            -- Common food name
    name_scientific VARCHAR(255), -- Scientific nomenclature
    description TEXT,             -- Manual documentation
    food_group VARCHAR(100),      -- Category (e.g., "Herbs and Spices")
    food_subgroup VARCHAR(100),   -- Subcategory (e.g., "Herbs")
    food_type VARCHAR(50),        -- Classification type
    itis_id INTEGER,              -- Integrated Taxonomic Information System ID
    wikipedia_id VARCHAR(255),    -- Wikipedia reference
    picture_file_name VARCHAR(255),
    picture_content_type VARCHAR(100),
    legacy_id INTEGER,
    created_at TIMESTAMP,
    updated_at TIMESTAMP
);
```

#### 2. Compounds Table
```sql
CREATE TABLE compounds (
    id INTEGER PRIMARY KEY,
    public_id VARCHAR(20),        -- e.g., "FDB000001"
    name VARCHAR(255),            -- Compound designation
    moldb_iupac VARCHAR(500),     -- IUPAC name
    moldb_smiles TEXT,            -- SMILES notation
    moldb_inchi TEXT,             -- InChI string
    moldb_inchikey VARCHAR(50),   -- InChI key
    moldb_formula VARCHAR(100),   -- Molecular formula
    moldb_average_mass DECIMAL(15,6),  -- Average molecular mass
    moldb_mono_mass DECIMAL(15,6),     -- Monoisotopic mass
    cas_number VARCHAR(20),       -- CAS registry number
    description TEXT,             -- Documentation
    export INTEGER,               -- Status: 1=active, 2=inactive
    state VARCHAR(20),            -- Physical state
    annotation_quality VARCHAR(20),
    created_at TIMESTAMP,
    updated_at TIMESTAMP
);
```

#### 3. Nutrients Table
```sql
CREATE TABLE nutrients (
    id INTEGER PRIMARY KEY,
    public_id VARCHAR(20),        -- Nutrient identifier
    name VARCHAR(255),            -- Nutrient designation
    description TEXT,             -- Nutrient documentation
    created_at TIMESTAMP,
    updated_at TIMESTAMP
);
```

#### 4. Contents Table (Association/Junction)
```sql
CREATE TABLE contents (
    id INTEGER PRIMARY KEY,
    food_id INTEGER REFERENCES foods(food_id),
    source_id INTEGER,            -- References compound or nutrient
    source_type VARCHAR(20),      -- "Nutrient" OR "Compound"
    orig_food_common_name VARCHAR(255),  -- Original food nomenclature
    orig_food_scientific_name VARCHAR(255),
    orig_food_part VARCHAR(100),
    orig_content DECIMAL(15,6),   -- Average or direct measurement
    orig_min DECIMAL(15,6),       -- Minimum concentration
    orig_max DECIMAL(15,6),       -- Maximum concentration
    orig_unit VARCHAR(20),        -- Unit (e.g., "mg/100g")
    orig_citation TEXT,
    orig_citation_type VARCHAR(50),
    orig_method VARCHAR(100),
    orig_unit_expression VARCHAR(50),
    standard_content DECIMAL(15,6),
    preparation_type VARCHAR(100),
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    FOREIGN KEY (food_id) REFERENCES foods(food_id),
    INDEX idx_source (source_id, source_type),
    INDEX idx_food (food_id)
);
```

---

## Food Categories (Groups and Subgroups)

### Main Food Groups

| Group | Example Subgroups |
|-------|-------------------|
| Herbs and Spices | Herbs, Spices |
| Vegetables | Cabbages, Onion-family, Root vegetables, Shoot vegetables, Leaf vegetables |
| Fruits | Tropical fruits, Berries, Citrus, Pome fruits, Stone fruits |
| Nuts | Tree nuts, Seeds |
| Cereals and Cereal Products | Cereals, Grains |
| Legumes | Beans, Lentils, Peas |
| Dairy | Milk, Cheese, Yogurt |
| Meat | Beef, Pork, Poultry |
| Fish and Seafood | Fish, Shellfish |
| Beverages | Alcoholic, Non-alcoholic |
| Confectioneries | Sweets, Chocolate |

---

## Compound Classification

### Status Codes

| Status | Value | Description |
|--------|-------|-------------|
| Quantified | 0 | Detected and quantified |
| Detected | 1 | Detected but not quantified |
| Expected | 2 | Expected but not quantified |
| Predicted | 3 | Predicted to exist |

### Compound Categories (Examples)

- Polyphenols (Flavonoids, Phenolic acids)
- Terpenoids
- Alkaloids
- Lipids and fatty acids
- Amino acids and peptides
- Carbohydrates
- Vitamins and cofactors
- Organic acids
- Aldehydes and ketones
- Alcohols
- Esters

---

## JSON Schemas

### Food Object
```json
{
  "public_id": "FOOD00001",
  "name": "Angelica",
  "name_scientific": "Angelica keiskei",
  "description": "A genus comprising approximately 60 species of biennial and perennial herbs in the Apiaceae family. Native to temperate and subarctic Northern Hemisphere regions.",
  "food_group": "Herbs and Spices",
  "food_subgroup": "Herbs",
  "taxonomy": {
    "superkingdom": "Eukaryota",
    "kingdom": "Viridiplantae",
    "phylum": "Streptophyta",
    "class": "Magnoliopsida",
    "order": "Apiales",
    "family": "Apiaceae",
    "genus": "Angelica",
    "species": "keiskei"
  },
  "external_links": {
    "wikipedia": "https://en.wikipedia.org/wiki/Angelica",
    "phenol_explorer": "phenol-explorer.eu/foods/XXX"
  }
}
```

### Compound Object
```json
{
  "public_id": "FDB000001",
  "name": "Cyanidin 3-O-glucoside",
  "cas_number": "7084-24-4",
  "moldb_formula": "C21H21O11+",
  "moldb_average_mass": 449.386,
  "moldb_mono_mass": 449.1084,
  "status": "Detected and quantified",
  "description": "A anthocyanin glycoside found in various berries and fruits.",
  "structure_image": "/compounds/FDB000001/structure.svg",
  "foods_containing": [
    {
      "food_id": "FOOD00234",
      "food_name": "Blueberry",
      "content": 25.3,
      "unit": "mg/100g"
    }
  ]
}
```

### Contents/Association Object
```json
{
  "food_id": 234,
  "food_public_id": "FOOD00234",
  "food_name": "Blueberry",
  "source_id": 1,
  "source_type": "Compound",
  "source_public_id": "FDB000001",
  "source_name": "Cyanidin 3-O-glucoside",
  "orig_content": 25.3,
  "orig_min": 15.2,
  "orig_max": 42.8,
  "orig_unit": "mg/100g",
  "orig_citation": "Rothwell et al., 2013",
  "preparation_type": "Raw"
}
```

### Nutrient Object
```json
{
  "public_id": "NUTR00001",
  "name": "Water",
  "description": "Total water content",
  "unit": "g/100g",
  "category": "Macronutrient"
}
```

---

## Data Access

### Web Interface
- Browse foods: https://foodb.ca/foods
- Browse compounds: https://foodb.ca/compounds
- Search: https://foodb.ca/unearth

### API Access
- **Status:** Contact required for API access
- **Contact:** Through FooDB website

### Bulk Downloads

| Format | Size (Compressed) | Content |
|--------|-------------------|---------|
| CSV | 952.52 MB | Compounds, proteins, contents, nutrients |
| XML | 6,438.08 MB | Full structured data |
| JSON | 86.66 MB | Structured data |
| MySQL Dump | 172.73 MB | Complete database |

**Download URL:** https://foodb.ca/downloads

### Spectral Data Downloads

| Data Type | Size |
|-----------|------|
| C-MS spectra (experimental) | Varies |
| C-MS spectra (predicted) | Varies |
| MS-MS spectra (experimental) | Varies |
| MS-MS spectra (predicted) | 1,198.74 MB |
| NMR spectra | Varies |
| NMR FID files | Varies |
| Peak lists | Varies |
| Images | 171.15 MB |

---

## Data Sources Integrated

FooDB integrates data from multiple sources:

| Source | Type | Data Contributed |
|--------|------|------------------|
| Phenol-Explorer | Database | Polyphenol content |
| KNApSAcK | Database | Metabolite relationships |
| HMDB | Database | Human metabolome data |
| DrugBank | Database | Compound information |
| PubChem | Database | Chemical structures |
| USDA SR | Database | Nutrient composition |
| Literature | Publications | Experimental values |

---

## Relationships

### Entity Relationship Diagram (Simplified)
```
foods (1) ----< contents >---- (N) compounds
                    |
                    v
                nutrients (via source_type)
```

### Foreign Key Relationships

| Table | Column | References |
|-------|--------|------------|
| contents | food_id | foods.food_id |
| contents | source_id + source_type | compounds.id OR nutrients.id |

---

## Use Cases

### 1. Find Compounds in a Food
```sql
SELECT c.name, co.orig_content, co.orig_unit
FROM contents co
JOIN compounds c ON co.source_id = c.id AND co.source_type = 'Compound'
WHERE co.food_id = (SELECT food_id FROM foods WHERE name = 'Blueberry');
```

### 2. Find Foods Containing a Compound
```sql
SELECT f.name, co.orig_content, co.orig_unit
FROM contents co
JOIN foods f ON co.food_id = f.food_id
WHERE co.source_id = (SELECT id FROM compounds WHERE public_id = 'FDB000001')
  AND co.source_type = 'Compound';
```

### 3. Get Nutrient Content for a Food
```sql
SELECT n.name, co.orig_content, co.orig_unit
FROM contents co
JOIN nutrients n ON co.source_id = n.id AND co.source_type = 'Nutrient'
WHERE co.food_id = (SELECT food_id FROM foods WHERE name = 'Apple');
```

---

## License

- **License Type:** Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
- **Academic Use:** Free with attribution
- **Commercial Use:** Requires permission
- **Contact:** Through FooDB website

---

## Citation

```
Scalbert, A., Andres-Lacueva, C., Arita, M., et al. (2011).
Databases on Food Phytochemicals and Their Health-Promoting Effects.
Journal of Agricultural and Food Chemistry, 59(1), 159-170.

FooDB Version 1.0. https://foodb.ca/
```

---

## Related Documents

- [43-79-SLEEP-LONGEVITY-NUTRI.md](../43-79-SLEEP-LONGEVITY-NUTRI.md) - References FooDB
- [43-52-NATURAL-PRODUCTS.md](../43-52-NATURAL-PRODUCTS.md) - Natural product databases
- [USDA-FOODDATA-CENTRAL.md](./USDA-FOODDATA-CENTRAL.md) - Related food database
