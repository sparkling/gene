---
title: "Nutrition Databases"
parent: ../_index.md
world: 3
last_updated: 2026-01-22
status: draft
---

# Nutrition Databases (World 3)

Databases for nutritional composition, food chemistry, and phytochemical content.

## Database Catalog

| Database | Type | Tier | Coverage | Access Method | Size |
|----------|------|------|----------|---------------|------|
| **FooDB** | Food Composition | 1 | 800+ foods, 70,000+ compounds | JSON API, Bulk | Comprehensive |
| **USDA FoodData Central** | Nutritional Data | 1 | 350,000+ foods | REST API, CSV | Official US database |
| **Phenol-Explorer** | Polyphenols | 2 | 500+ foods, 500+ polyphenols | Web, Download | Specialized |
| **PhytoHub** | Phytochemicals | 2 | 1,800+ metabolites | Web, Download | Metabolite-focused |
| **HMDB** | Human Metabolome | 2 | 114,000+ metabolites | XML, API | Comprehensive |
| **FRIDA** | Food-Risk Component | 2 | 300+ foods, toxins | Web | Safety-focused |
| **FoodOmicsGR** | Greek Foods | 3 | Greek diet, omics data | Web | Regional |
| **EuroFIR** | European Foods | 3 | 30+ national databases | Platform | European focus |
| **AUSNUT** | Australian Foods | 3 | Australian foods | Download | Regional |
| **UK Food Composition** | UK Foods | 3 | UK foods | Download | Regional |

## Primary Use Cases

### Tier 1 (MVP) Focus

#### FooDB
- **Purpose**: Comprehensive food composition database
- **Content**:
  - 800+ foods
  - 70,000+ food compounds
  - Chemical structures
  - Nutritional data
- **Use**: Link genes to food compounds
- **Access**: https://foodb.ca/
  - REST API: `https://foodb.ca/api/v1/compounds`
  - JSON export available
- **Integration**: Links to PubChem, HMDB, ChEBI

#### USDA FoodData Central
- **Purpose**: Official US nutritional database
- **Content**:
  - 350,000+ foods
  - Foundation Foods: nutrient composition
  - SR Legacy: Standard Reference
  - Branded Foods: packaged products
- **Use**: Nutritional recommendations, dietary analysis
- **Access**: https://fdc.nal.usda.gov/
  - REST API: `https://api.nal.usda.gov/fdc/v1/`
  - Bulk CSV downloads
- **API Key**: Required (free registration)

## Data Integration Workflow

```
Food Item
    ↓
USDA FoodData (nutritional profile)
    ↓
FooDB (compound composition)
    ↓
┌──────────┬──────────┐
│          │          │
Phenol-    PhytoHub   HMDB
Explorer   │          │
│          │          │
Polyphenols Phyto-   Metabolites
           chemicals
│          │          │
└──────────┴──────────┘
         ↓
   Gene Targets
         ↓
    Pathways
```

## Tier 2 Databases

### Phenol-Explorer
- **Purpose**: Polyphenol content in foods
- **Content**:
  - 500+ foods
  - 500+ polyphenols
  - Bioavailability data
  - Metabolism information
- **Use**: Polyphenol-gene interactions
- **Access**: http://phenol-explorer.eu/
  - Download CSV files
  - No official API

### PhytoHub
- **Purpose**: Dietary phytochemical metabolism
- **Content**:
  - 1,800+ phytochemical metabolites
  - Food sources
  - Transformation pathways
- **Use**: Track phytochemical metabolism
- **Access**: http://phytohub.eu/
  - Download SDF, CSV
  - Web interface

### HMDB (Human Metabolome Database)
- **Purpose**: Comprehensive metabolome reference
- **Content**:
  - 114,000+ metabolites
  - Food metabolites
  - Drug metabolites
  - Endogenous metabolites
- **Use**: Metabolite identification, pathway analysis
- **Access**: https://hmdb.ca/
  - XML download
  - REST API
  - SPARQL endpoint

## Access Methods

### FooDB
```bash
# API access (no key required)
curl "https://foodb.ca/api/v1/compounds/FDB012345"
curl "https://foodb.ca/api/v1/foods/FOOD00001"

# Bulk download
wget https://foodb.ca/downloads/compounds.json.zip
wget https://foodb.ca/downloads/foods.json.zip
```

### USDA FoodData Central
```bash
# API key required (register at https://fdc.nal.usda.gov/api-key-signup.html)
API_KEY="your_api_key_here"

# Search foods
curl "https://api.nal.usda.gov/fdc/v1/foods/search?query=apple&api_key=$API_KEY"

# Get food details
curl "https://api.nal.usda.gov/fdc/v1/food/12345?api_key=$API_KEY"

# Bulk download (no API key needed)
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_csv_2023-10-26.zip
```

### Phenol-Explorer
```bash
# Download database files
wget http://phenol-explorer.eu/contents/download
# Manual download required - no direct links
```

### PhytoHub
```bash
# Download compound data
wget http://phytohub.eu/files/PhytoHub_metabolites.sdf
wget http://phytohub.eu/files/PhytoHub_foods.csv
```

## Data Schema Examples

### FooDB Food Entry
```json
{
  "id": "FOOD00001",
  "name": "Apple",
  "scientific_name": "Malus domestica",
  "compounds": [
    {
      "id": "FDB012345",
      "name": "Quercetin",
      "content_mg_per_100g": 4.3,
      "pubchem_id": "5280343"
    }
  ]
}
```

### USDA Food Entry
```json
{
  "fdcId": 171688,
  "description": "Apple, raw",
  "foodNutrients": [
    {
      "nutrientId": 1003,
      "nutrientName": "Protein",
      "value": 0.26,
      "unitName": "g"
    }
  ]
}
```

## Nutrient-Gene Mapping Strategy

### Direct Mapping
```
Food (USDA) → Compound (FooDB) → Gene Target → Pathway
```

### Example: Apple → Quercetin → Genes
```
Apple (USDA)
    ↓
Quercetin (FooDB, Phenol-Explorer)
    ↓
Target Genes (KEAP1, NRF2, SIRT1)
    ↓
Pathways (Antioxidant response, inflammation)
    ↓
Health Benefits (Anti-inflammatory, neuroprotective)
```

## Regional Databases (Tier 3)

### EuroFIR
- **Coverage**: 30+ European national food databases
- **Access**: Platform subscription required
- **Use**: European dietary analysis

### AUSNUT
- **Coverage**: Australian foods
- **Access**: https://www.foodstandards.gov.au/science-data/food-composition-databases
- **Format**: Excel, CSV

### UK Food Composition
- **Coverage**: UK foods
- **Access**: https://www.gov.uk/government/publications/composition-of-foods-integrated-dataset-cofid
- **Format**: Excel, CSV

## Update Strategy

### Automated Updates (Tier 1)
- **USDA FoodData Central**: Quarterly releases (check monthly)
- **FooDB**: Check for updates biannually

### Manual Updates (Tier 2)
- **Phenol-Explorer**: Check annually
- **PhytoHub**: Check annually
- **HMDB**: Major releases every 1-2 years

## Storage Estimates

| Database | Storage Required | Format |
|----------|------------------|--------|
| FooDB | 1-2 GB | JSON, XML |
| USDA FoodData Central | 500 MB - 1 GB | CSV, JSON |
| Phenol-Explorer | 50-100 MB | CSV, Excel |
| PhytoHub | 100-200 MB | SDF, CSV |
| HMDB | 5-10 GB | XML, CSV |

## Data Quality Considerations

### Standardization Challenges
- Food name variations (regional, brand names)
- Portion size differences
- Preparation method effects
- Seasonal variations in compound content

### Integration Considerations
- Mapping food items across databases
- Handling compound concentration ranges
- Accounting for bioavailability
- Regional dietary patterns

## Compound Coverage

| Database | Macronutrients | Micronutrients | Phytochemicals | Metabolites |
|----------|----------------|----------------|----------------|-------------|
| USDA | ✓✓✓ | ✓✓✓ | ✓ | - |
| FooDB | ✓✓ | ✓✓✓ | ✓✓✓ | ✓✓ |
| Phenol-Explorer | - | - | ✓✓✓ | ✓✓ |
| PhytoHub | - | - | ✓✓✓ | ✓✓✓ |
| HMDB | ✓ | ✓✓ | ✓✓ | ✓✓✓ |

## Navigation

- **Parent**: [Database Sources](../_index.md)
- **Related**: [Compounds](../compounds/_index.md), [Traditional Medicine](../traditional/_index.md)
