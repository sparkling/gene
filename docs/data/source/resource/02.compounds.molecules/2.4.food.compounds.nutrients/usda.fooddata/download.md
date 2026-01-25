---
id: download-usda.fooddata
title: "USDA FoodData Central Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# USDA FoodData Central Download Instructions

## Quick Start

```bash
# Get API key (free)
# https://api.data.gov/signup/

# Query the API
curl "https://api.nal.usda.gov/fdc/v1/foods/search?api_key=YOUR_KEY&query=apple" | jq

# Or download bulk data
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_csv_2024-04.zip
```

## Prerequisites

- **API key** (free registration at api.data.gov)
- **wget** or **curl** for downloads
- ~5-10 GB disk space for full dataset
- JSON parser (jq) or CSV tools

## Registration

Free API key required: https://api.data.gov/signup/

## Download Methods

### Method 1: REST API (Recommended for Targeted Access)

```bash
# Get API key first
API_KEY="your_api_key_here"

# Search foods
curl "https://api.nal.usda.gov/fdc/v1/foods/search?api_key=$API_KEY&query=apple" | jq

# Get specific food by FDC ID
curl "https://api.nal.usda.gov/fdc/v1/food/167512?api_key=$API_KEY" | jq

# Search with filters
curl "https://api.nal.usda.gov/fdc/v1/foods/search?api_key=$API_KEY&query=chicken&dataType=Foundation&pageSize=25" | jq

# List foods with pagination
curl "https://api.nal.usda.gov/fdc/v1/foods/list?api_key=$API_KEY&pageNumber=1&pageSize=50" | jq
```

### Method 2: Bulk CSV Downloads

```bash
# Navigate to download page
# https://fdc.nal.usda.gov/download-datasets.html

# Download full dataset (all data types)
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_csv_2024-04.zip
unzip FoodData_Central_csv_2024-04.zip -d fooddata

# Download specific data types
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_foundation_food_csv_2024-04.zip
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_sr_legacy_food_csv_2024-04.zip
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_branded_food_csv_2024-04.zip
```

### Method 3: JSON Downloads

```bash
# JSON format downloads (same content, different format)
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_foundation_food_json_2024-04.zip
unzip FoodData_Central_foundation_food_json_2024-04.zip -d foundation_json
```

### Method 4: Supporting Data Files

```bash
# Download supporting files
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_Supporting_Data_csv_2024-04.zip
# Contains:
# - nutrient.csv (nutrient definitions)
# - food_category.csv (categories)
# - measure_unit.csv (portion units)
```

## File Inventory

### CSV Downloads

| File | Size | Description |
|------|------|-------------|
| food.csv | ~50 MB | Core food records |
| food_nutrient.csv | ~2 GB | Nutrient values |
| nutrient.csv | ~50 KB | Nutrient definitions |
| food_portion.csv | ~20 MB | Serving sizes |
| branded_food.csv | ~500 MB | Branded products |
| foundation_food.csv | ~1 MB | Foundation foods |

### JSON Downloads

| File | Size | Description |
|------|------|-------------|
| FoodData_Central_*.json | ~3 GB | Complete dataset |
| foundation_food.json | ~200 MB | Foundation only |
| sr_legacy_food.json | ~400 MB | SR Legacy only |

## Post-Download Processing

```bash
# Preview CSV structure
head -5 food.csv
head -5 food_nutrient.csv

# Count foods by data type
cut -d',' -f3 food.csv | sort | uniq -c

# Load into SQLite
sqlite3 fooddata.db << 'EOF'
.mode csv
.import food.csv food
.import food_nutrient.csv food_nutrient
.import nutrient.csv nutrient

CREATE INDEX idx_food_id ON food(fdc_id);
CREATE INDEX idx_fn_food ON food_nutrient(fdc_id);
CREATE INDEX idx_fn_nutrient ON food_nutrient(nutrient_id);
EOF

# Query for specific food
sqlite3 fooddata.db "SELECT * FROM food WHERE description LIKE '%apple%' LIMIT 10;"

# Get nutrient profile
sqlite3 fooddata.db << 'EOF'
SELECT n.name, fn.amount, n.unit_name
FROM food_nutrient fn
JOIN nutrient n ON fn.nutrient_id = n.id
WHERE fn.fdc_id = 167512
ORDER BY n.rank;
EOF

# Calculate macros for a food
sqlite3 fooddata.db << 'EOF'
SELECT
  SUM(CASE WHEN nutrient_id=1008 THEN amount END) as kcal,
  SUM(CASE WHEN nutrient_id=1003 THEN amount END) as protein_g,
  SUM(CASE WHEN nutrient_id=1004 THEN amount END) as fat_g,
  SUM(CASE WHEN nutrient_id=1005 THEN amount END) as carbs_g
FROM food_nutrient
WHERE fdc_id = 167512;
EOF
```

## API Query Examples

```bash
# Search branded foods only
curl "https://api.nal.usda.gov/fdc/v1/foods/search?api_key=$API_KEY&query=yogurt&dataType=Branded" | jq

# Get nutrients for multiple foods
curl -X POST "https://api.nal.usda.gov/fdc/v1/foods?api_key=$API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"fdcIds": [167512, 168193, 169148]}' | jq

# Search by GTIN/UPC barcode
curl "https://api.nal.usda.gov/fdc/v1/foods/search?api_key=$API_KEY&query=049000042566" | jq
```

## Verification

```bash
# Check CSV integrity
wc -l food.csv
wc -l food_nutrient.csv

# Verify column count
head -1 food.csv | tr ',' '\n' | wc -l

# Check for required fields
head -1 food_nutrient.csv

# Verify FDC ID coverage
cut -d',' -f1 food.csv | sort | uniq | wc -l
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| v13.6 | 2025-11-20 | ~5 GB | Current |
| v13.5 | 2025-09-18 | ~5 GB | Archived |
| v12.x | 2025-01-23 | ~5 GB | Archived |

### Version Notes

FoodData Central v13.6 (latest) includes:
- Foundation Foods (high quality reference data)
- SR Legacy (USDA Standard Reference)
- FNDDS (What We Eat in America survey)
- Branded Foods (USDA Global Branded Food Products Database)
- FNS Child Nutrition Database now hosted on FDC

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://api.nal.usda.gov/fdc/v1` |
| Rate Limit | 1000 req/hour (registered) |
| Auth Required | Yes (free API key) |
| Documentation | https://fdc.nal.usda.gov/api-guide.html |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Full releases | Monthly |
| API data | Continuous |
| Foundation Foods | Periodic |
| Branded Foods | Continuous |

## API Rate Limits

| Tier | Limit |
|------|-------|
| Default | 1,000 requests/hour |
| Demo key | 30 requests/hour |

## Common Issues

- **Large files**: Use streaming for food_nutrient.csv (2GB+)
- **CSV parsing**: Some fields contain commas; use proper CSV parser
- **Missing values**: Not all nutrients measured for all foods
- **Legacy NDB numbers**: Use FDC ID as primary key

## Nutrient ID Reference

| ID | Nutrient |
|----|----------|
| 1008 | Energy (kcal) |
| 1003 | Protein |
| 1004 | Total fat |
| 1005 | Carbohydrate |
| 1079 | Fiber |
| 1087 | Calcium |
| 1089 | Iron |
| 1104 | Vitamin A |
| 1162 | Vitamin C |

## License

Public Domain (US Government) - Free for any use
