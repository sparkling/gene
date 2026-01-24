---
id: download-foodb
title: "FooDB Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# FooDB Download Instructions

## Quick Start

```bash
# Download FooDB CSV export
wget https://foodb.ca/public/system/downloads/foodb_2020_04_07_csv.tar.gz
tar -xzf foodb_2020_04_07_csv.tar.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **tar** for extraction
- Approximately 2GB disk space

## No Registration Required

FooDB data is freely available for academic and non-commercial use.

## Download Methods

### Method 1: CSV Export (Recommended)

```bash
# Download latest CSV export
wget https://foodb.ca/public/system/downloads/foodb_2020_04_07_csv.tar.gz

# Extract
tar -xzf foodb_2020_04_07_csv.tar.gz
cd foodb_2020_04_07_csv/

# List files
ls -la *.csv
```

### Method 2: XML Export

```bash
# Download XML format
wget https://foodb.ca/public/system/downloads/foodb_2020_04_07_xml.tar.gz

# Extract
tar -xzf foodb_2020_04_07_xml.tar.gz
```

### Method 3: SDF Structures

```bash
# Download compound structures
wget https://foodb.ca/public/system/downloads/foodb_2020_04_07_sdf.tar.gz

# Extract
tar -xzf foodb_2020_04_07_sdf.tar.gz
```

### Method 4: JSON Export

```bash
# Download JSON format
wget https://foodb.ca/public/system/downloads/foodb_2020_04_07_json.tar.gz

# Extract
tar -xzf foodb_2020_04_07_json.tar.gz
```

### Method 5: Web Scraping API

```bash
# Get food details
curl "https://foodb.ca/foods/FOOD00001.json" -o apple.json

# Get compound details
curl "https://foodb.ca/compounds/FDB000001.json" -o water.json

# Note: API has rate limits; be respectful
```

## File Inventory

### CSV Files

| File | Size | Description |
|------|------|-------------|
| Food.csv | ~2 MB | Food items |
| Compound.csv | ~20 MB | Chemical compounds |
| Content.csv | ~50 MB | Food-compound relationships |
| Nutrient.csv | ~1 MB | Nutrient definitions |
| FoodNutrient.csv | ~10 MB | Food-nutrient values |
| FoodTaxonomy.csv | ~500 KB | Food classification |
| CompoundExternalDescriptor.csv | ~5 MB | External IDs |
| Reference.csv | ~10 MB | Literature references |
| CompoundFlavor.csv | ~1 MB | Flavor associations |
| CompoundHealthEffect.csv | ~2 MB | Health effects |

### Other Formats

| File | Size | Description |
|------|------|-------------|
| foodb_*.xml.tar.gz | ~500 MB | XML format |
| foodb_*.sdf.tar.gz | ~100 MB | SDF structures |
| foodb_*.json.tar.gz | ~800 MB | JSON format |

## Post-Download Processing

```bash
# Parse food data
python3 << 'EOF'
import pandas as pd

# Load foods
foods = pd.read_csv('Food.csv')
print(f"Foods: {len(foods)}")
print(foods.columns.tolist())

# Load compounds
compounds = pd.read_csv('Compound.csv')
print(f"Compounds: {len(compounds)}")

# Load food-compound relationships
contents = pd.read_csv('Content.csv')
print(f"Food-compound entries: {len(contents)}")

# Find compounds in a specific food (e.g., Apple)
apple = foods[foods['name'].str.contains('Apple', case=False, na=False)].iloc[0]
apple_compounds = contents[contents['food_id'] == apple['id']]
print(f"Compounds in apple: {len(apple_compounds)}")
EOF

# Create food-compound matrix
python3 << 'EOF'
import pandas as pd

foods = pd.read_csv('Food.csv')
compounds = pd.read_csv('Compound.csv')
contents = pd.read_csv('Content.csv')

# Merge for readable names
merged = contents.merge(
    foods[['id', 'name']],
    left_on='food_id',
    right_on='id',
    suffixes=('', '_food')
).merge(
    compounds[['id', 'name']],
    left_on='source_id',
    right_on='id',
    suffixes=('_food', '_compound')
)

# Save simplified version
merged[['name_food', 'name_compound', 'orig_content', 'orig_unit']].to_csv(
    'food_compound_simple.tsv', sep='\t', index=False
)
EOF

# Extract nutrient data
python3 << 'EOF'
import pandas as pd

foods = pd.read_csv('Food.csv')
nutrients = pd.read_csv('Nutrient.csv')
food_nutrients = pd.read_csv('FoodNutrient.csv')

# Merge for complete view
nutrient_data = food_nutrients.merge(
    foods[['id', 'name']],
    left_on='food_id',
    right_on='id'
).merge(
    nutrients[['id', 'name']],
    left_on='nutrient_id',
    right_on='id',
    suffixes=('_food', '_nutrient')
)

# Pivot to wide format
wide = nutrient_data.pivot(
    index='name_food',
    columns='name_nutrient',
    values='orig_content'
)
wide.to_csv('food_nutrient_matrix.csv')
EOF
```

## Verification

```bash
# Check CSV structure
head -5 Food.csv

# Count records
wc -l *.csv

# Check for specific food
grep -i "apple" Food.csv

# Verify relationships
awk -F',' 'NR==1 || $2==1' Content.csv | head
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major updates | Every 1-2 years |
| Bug fixes | As reported |

## Common Issues

- **Large Content.csv**: Contains millions of food-compound relationships
- **Null values**: Many compounds lack concentration data
- **Unit inconsistency**: Various units used; normalize carefully
- **Encoding**: UTF-8 encoding; handle special characters
- **ID mapping**: Use CompoundExternalDescriptor for cross-references

## Data Categories

| Category | Description |
|----------|-------------|
| Foods | 892 food items |
| Compounds | 70,000+ compounds |
| Nutrients | 100+ nutrients |
| Health Effects | 300+ health associations |
| Flavors | Flavor profiles |

## Food Classification

```bash
# Extract food taxonomy
python3 << 'EOF'
import pandas as pd

taxonomy = pd.read_csv('FoodTaxonomy.csv')
print(taxonomy[['id', 'name', 'parent_id']].head(20))

# Build hierarchy
foods = pd.read_csv('Food.csv')
merged = foods.merge(taxonomy, left_on='food_type', right_on='id', how='left')
print(merged[['name_x', 'name_y']].value_counts())
EOF
```

## Integration Examples

```bash
# Map to PubChem
python3 << 'EOF'
import pandas as pd

xrefs = pd.read_csv('CompoundExternalDescriptor.csv')
pubchem = xrefs[xrefs['external_id_source'] == 'PubChem']
pubchem[['compound_id', 'external_id']].to_csv('foodb_pubchem_map.tsv', sep='\t', index=False)
EOF

# Map to HMDB
python3 << 'EOF'
import pandas as pd

xrefs = pd.read_csv('CompoundExternalDescriptor.csv')
hmdb = xrefs[xrefs['external_id_source'] == 'HMDB']
hmdb[['compound_id', 'external_id']].to_csv('foodb_hmdb_map.tsv', sep='\t', index=False)
EOF
```

## Related Resources

- [HMDB](../../6.4.metabolomics/hmdb/) - Human metabolome
- [PubChem](../../02.compounds.molecules/2.6.chemical.ontology.classification/pubchem/) - Chemical database
- [USDA FoodData](../../2.4.food.compounds.nutrients/usda.fooddata/) - Nutrient database
