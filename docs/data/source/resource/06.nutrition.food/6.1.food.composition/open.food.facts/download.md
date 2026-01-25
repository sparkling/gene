---
id: download-open-food-facts
title: "Open Food Facts Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Open Food Facts Download Instructions

## Quick Start

```bash
# Download CSV export (compressed)
wget https://static.openfoodfacts.org/data/en.openfoodfacts.org.products.csv.gz
gunzip en.openfoodfacts.org.products.csv.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for extraction
- Approximately 5-10GB disk space for uncompressed data

## No Registration Required

Open Food Facts data is freely available under the Open Database License (ODbL).

## Download Methods

### Method 1: CSV Export (Recommended for Analysis)

```bash
# Download CSV (all products, flat format)
wget https://static.openfoodfacts.org/data/en.openfoodfacts.org.products.csv.gz

# Extract
gunzip en.openfoodfacts.org.products.csv.gz

# Check file size and structure
wc -l en.openfoodfacts.org.products.csv
head -1 en.openfoodfacts.org.products.csv | tr '\t' '\n' | nl
```

### Method 2: JSONL Export (Full Data)

```bash
# Download JSONL (complete product data)
wget https://static.openfoodfacts.org/data/openfoodfacts-products.jsonl.gz

# Extract
gunzip openfoodfacts-products.jsonl.gz

# Process line by line (memory efficient)
head -1 openfoodfacts-products.jsonl | jq '.'
```

### Method 3: Parquet (ML/Analytics)

```bash
# Using Hugging Face datasets
pip install datasets

python3 << 'EOF'
from datasets import load_dataset

# Load dataset (streams from HuggingFace)
ds = load_dataset("openfoodfacts/product-database", split="train")
print(f"Total products: {len(ds)}")

# Convert to pandas for analysis
df = ds.to_pandas()
df.to_parquet("openfoodfacts.parquet")
EOF
```

### Method 4: RDF/Linked Data

```bash
# Download RDF dump
wget https://static.openfoodfacts.org/data/openfoodfacts-products.rdf.gz

# Extract
gunzip openfoodfacts-products.rdf.gz
```

### Method 5: MongoDB Dump

```bash
# Download MongoDB export
wget https://static.openfoodfacts.org/data/openfoodfacts-mongodbdump.tar.gz

# Extract
tar -xzf openfoodfacts-mongodbdump.tar.gz

# Restore to local MongoDB
mongorestore dump/off
```

### Method 6: API for Specific Products

```bash
# Get single product by barcode
curl "https://world.openfoodfacts.org/api/v2/product/3017620422003" -o nutella.json

# Search products
curl "https://world.openfoodfacts.org/api/v2/search?categories_tags=en:breakfast-cereals&page_size=100" -o cereals.json

# Note: For bulk data, use downloads instead of API
```

### Method 7: Delta Updates

```bash
# Download daily delta (new/modified products since last full export)
wget https://static.openfoodfacts.org/data/delta/products.$(date +%Y-%m-%d).jsonl.gz
```

## File Inventory

### Primary Exports

| File | Size | Description |
|------|------|-------------|
| en.openfoodfacts.org.products.csv.gz | ~800 MB | CSV format |
| openfoodfacts-products.jsonl.gz | ~1.2 GB | JSONL format |
| openfoodfacts-mongodbdump.tar.gz | ~1.5 GB | MongoDB dump |

### Alternative Formats

| File | Size | Description |
|------|------|-------------|
| openfoodfacts-products.rdf.gz | ~5 GB | RDF/Linked data |
| HuggingFace Parquet | ~600 MB | ML-optimized |

### Country-Specific Exports

| Country | URL Pattern |
|---------|-------------|
| United States | us.openfoodfacts.org.products.csv.gz |
| France | fr.openfoodfacts.org.products.csv.gz |
| Germany | de.openfoodfacts.org.products.csv.gz |

## Post-Download Processing

```bash
# Load into Python for analysis
python3 << 'EOF'
import pandas as pd

# Load CSV
df = pd.read_csv('en.openfoodfacts.org.products.csv',
                 sep='\t',
                 low_memory=False,
                 encoding='utf-8')

print(f"Total products: {len(df)}")
print(f"Columns: {len(df.columns)}")

# Check Nutri-Score distribution
print(df['nutriscore_grade'].value_counts())

# Filter by country
us_products = df[df['countries_tags'].str.contains('en:united-states', na=False)]
print(f"US products: {len(us_products)}")

# Save subset
us_products[['code', 'product_name', 'brands', 'nutriscore_grade', 'nova_group']].to_csv(
    'us_products_subset.tsv', sep='\t', index=False
)
EOF

# Process JSONL (streaming)
python3 << 'EOF'
import json
import gzip

nutriscore_counts = {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}

with gzip.open('openfoodfacts-products.jsonl.gz', 'rt', encoding='utf-8') as f:
    for line in f:
        product = json.loads(line)
        grade = product.get('nutriscore_grade', '')
        if grade in nutriscore_counts:
            nutriscore_counts[grade] += 1

print(nutriscore_counts)
EOF

# Extract specific fields using jq
zcat openfoodfacts-products.jsonl.gz | head -1000 | \
  jq -r '[.code, .product_name, .nutriscore_grade] | @tsv' > sample.tsv
```

## Verification

```bash
# Check CSV structure
head -5 en.openfoodfacts.org.products.csv | cut -f1-10

# Count records
wc -l en.openfoodfacts.org.products.csv

# Check for specific product
grep "3017620422003" en.openfoodfacts.org.products.csv

# Verify JSONL
zcat openfoodfacts-products.jsonl.gz | head -1 | jq 'keys'
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Full exports | Daily (02:00 UTC) |
| Delta exports | Daily |
| Real-time API | Continuous |

## Common Issues

- **Large file size**: Use streaming/chunked processing
- **UTF-8 encoding**: Some product names contain special characters
- **Missing fields**: Not all products have complete data
- **Tab-separated**: CSV uses tabs, not commas
- **NULL values**: Empty strings, not NULL keyword
- **Duplicate barcodes**: Some products have multiple entries

## Key Columns in CSV

| Column | Description |
|--------|-------------|
| code | Product barcode |
| product_name | Product name |
| brands | Brand name |
| categories_tags | Product categories |
| countries_tags | Countries where sold |
| nutriscore_grade | Nutri-Score (a-e) |
| nova_group | NOVA classification (1-4) |
| ecoscore_grade | Eco-Score (a-e) |
| energy-kcal_100g | Calories per 100g |
| fat_100g | Fat per 100g |
| sugars_100g | Sugars per 100g |
| proteins_100g | Protein per 100g |
| salt_100g | Salt per 100g |

## API Integration Examples

```bash
# Build product lookup
python3 << 'EOF'
import pandas as pd
import sqlite3

# Load CSV
df = pd.read_csv('en.openfoodfacts.org.products.csv', sep='\t', low_memory=False)

# Create SQLite database
conn = sqlite3.connect('openfoodfacts.db')
df[['code', 'product_name', 'brands', 'categories_tags',
    'nutriscore_grade', 'nova_group', 'energy-kcal_100g']].to_sql(
    'products', conn, if_exists='replace', index=False
)

# Create index
conn.execute('CREATE INDEX IF NOT EXISTS idx_code ON products(code)')
conn.execute('CREATE INDEX IF NOT EXISTS idx_nutriscore ON products(nutriscore_grade)')
conn.close()

print("SQLite database created: openfoodfacts.db")
EOF

# Query SQLite database
sqlite3 openfoodfacts.db "SELECT COUNT(*) FROM products WHERE nutriscore_grade = 'a'"
```

## Related Resources

- [Schema Documentation](./schema.md)
- [FooDB](../foodb/) - Food compound database
- [USDA FoodData](../usda.fooddata/) - US nutrient database
