---
id: download-dsld
title: "DSLD Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# DSLD (Dietary Supplement Label Database) Download Instructions

## Quick Start

```bash
# Get API key first: https://api.ods.od.nih.gov/signup

# Search for products
curl -H "X-Api-Key: YOUR_API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/browse?q=vitamin+d&pagesize=100" \
  -o vitamin_d_products.json
```

## Prerequisites

- **API Key:** Free registration required
- **curl** or Python requests library
- **jq** for JSON processing (optional)
- Disk space: ~1-5GB for full dataset

## Registration Required

1. Visit: https://api.ods.od.nih.gov/signup
2. Register with email address
3. Receive API key via email
4. Use key in X-Api-Key header

## Download Methods

### Method 1: API Access (Recommended)

```bash
# Set your API key
export DSLD_API_KEY="your-api-key-here"

# Search products
curl -H "X-Api-Key: $DSLD_API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/browse?pagesize=100" \
  -o products_page1.json

# Get specific product details
curl -H "X-Api-Key: $DSLD_API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/label/123456" \
  -o product_123456.json

# Get ingredients for a product
curl -H "X-Api-Key: $DSLD_API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/ingredients/123456" \
  -o ingredients_123456.json
```

### Method 2: Bulk Download (Web Interface)

```bash
# Navigate to bulk download page
# https://dsld.od.nih.gov/download

# Download options available:
# - Full database export (CSV/JSON)
# - Ingredient dictionary
# - Brand list
```

### Method 3: Paginated API Download

```bash
# Download all products in batches
python3 << 'EOF'
import requests
import json
import time

API_KEY = "your-api-key-here"
BASE_URL = "https://api.ods.od.nih.gov/dsld/v8"
HEADERS = {"X-Api-Key": API_KEY}

def download_all_products():
    all_products = []
    offset = 0
    page_size = 100

    while True:
        url = f"{BASE_URL}/browse?pagesize={page_size}&offset={offset}"
        response = requests.get(url, headers=HEADERS)
        data = response.json()

        products = data.get('results', [])
        if not products:
            break

        all_products.extend(products)
        print(f"Downloaded {len(all_products)} products...")

        offset += page_size
        time.sleep(0.5)  # Rate limiting

        # Optional: limit for testing
        if offset >= 1000:
            break

    return all_products

products = download_all_products()

with open('dsld_products.json', 'w') as f:
    json.dump(products, f, indent=2)

print(f"Total products: {len(products)}")
EOF
```

### Method 4: Download by Category

```bash
# Download products by ingredient category
python3 << 'EOF'
import requests
import json
import time

API_KEY = "your-api-key-here"
BASE_URL = "https://api.ods.od.nih.gov/dsld/v8"
HEADERS = {"X-Api-Key": API_KEY}

categories = [
    "vitamin d",
    "vitamin c",
    "omega-3",
    "probiotics",
    "turmeric",
    "melatonin"
]

for category in categories:
    url = f"{BASE_URL}/browse?ingredient={category.replace(' ', '%20')}&pagesize=100"
    response = requests.get(url, headers=HEADERS)
    data = response.json()

    filename = f"dsld_{category.replace(' ', '_')}.json"
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"{category}: {data.get('totalCount', 0)} products")
    time.sleep(1)
EOF
```

### Method 5: Download Full Product Details

```bash
# Download complete label data for products
python3 << 'EOF'
import requests
import json
import time

API_KEY = "your-api-key-here"
BASE_URL = "https://api.ods.od.nih.gov/dsld/v8"
HEADERS = {"X-Api-Key": API_KEY}

# First get product IDs
browse_url = f"{BASE_URL}/browse?pagesize=50"
response = requests.get(browse_url, headers=HEADERS)
products = response.json().get('results', [])

# Then get full details for each
detailed_products = []
for product in products[:10]:  # Limit for testing
    dsld_id = product.get('dsldId')
    detail_url = f"{BASE_URL}/label/{dsld_id}"
    detail_response = requests.get(detail_url, headers=HEADERS)

    if detail_response.status_code == 200:
        detailed_products.append(detail_response.json())
        print(f"Downloaded: {product.get('productName', 'Unknown')[:50]}")

    time.sleep(0.5)

with open('dsld_detailed.json', 'w') as f:
    json.dump(detailed_products, f, indent=2)
EOF
```

## File Inventory

### API Response Fields

| Field | Description |
|-------|-------------|
| dsldId | Unique product identifier |
| productName | Full product name |
| brandName | Brand name |
| netContents | Package size |
| servingSize | Serving size |
| ingredients | List of ingredients |
| manufacturer | Company information |

### Bulk Download Files

| File | Description |
|------|-------------|
| products.csv | All products basic info |
| ingredients.csv | All ingredients |
| brands.csv | Brand directory |
| manufacturers.csv | Manufacturer directory |

## Post-Download Processing

```bash
# Convert JSON to CSV
python3 << 'EOF'
import json
import pandas as pd

with open('dsld_products.json') as f:
    products = json.load(f)

# Flatten to DataFrame
df = pd.json_normalize(products)
df.to_csv('dsld_products.csv', index=False)

print(f"Converted {len(df)} products to CSV")
print(f"Columns: {list(df.columns)}")
EOF

# Process ingredients separately
python3 << 'EOF'
import json
import pandas as pd

with open('dsld_detailed.json') as f:
    products = json.load(f)

ingredients = []
for product in products:
    dsld_id = product.get('dsldId')
    for ing in product.get('ingredients', []):
        ing['dsldId'] = dsld_id
        ingredients.append(ing)

df = pd.DataFrame(ingredients)
df.to_csv('dsld_ingredients.csv', index=False)

print(f"Extracted {len(df)} ingredient records")
EOF

# Create SQLite database
python3 << 'EOF'
import json
import sqlite3
import pandas as pd

conn = sqlite3.connect('dsld.db')

# Load products
with open('dsld_products.json') as f:
    products = json.load(f)

df_products = pd.json_normalize(products)
df_products.to_sql('products', conn, if_exists='replace', index=False)

# Create indexes
conn.execute('CREATE INDEX IF NOT EXISTS idx_brand ON products(brandName)')
conn.execute('CREATE INDEX IF NOT EXISTS idx_dsld ON products(dsldId)')

conn.close()
print("SQLite database created: dsld.db")
EOF
```

## Verification

```bash
# Check API response
curl -H "X-Api-Key: $DSLD_API_KEY" \
  "https://api.ods.od.nih.gov/dsld/v8/browse?pagesize=1" | jq

# Verify downloaded data
jq 'length' dsld_products.json

# Check specific product
jq '.[] | select(.dsldId == 123456)' dsld_products.json
```

## Update Schedule

| Content | Frequency |
|---------|-----------|
| New products | Monthly |
| Label updates | As received |
| Database refresh | Continuous |

## Rate Limits

| Limit Type | Value |
|------------|-------|
| Requests per second | ~2 |
| Requests per day | ~10,000 |
| Bulk download | No limit |

## Common Issues

- **401 Unauthorized**: Check API key is valid and properly formatted
- **429 Too Many Requests**: Add delays between requests
- **Missing fields**: Not all products have complete data
- **Encoding**: UTF-8 with some special characters
- **Null values**: Some fields may be null or empty strings

## API Endpoints Reference

| Endpoint | Method | Description |
|----------|--------|-------------|
| /browse | GET | Search/list products |
| /label/{id} | GET | Full label details |
| /ingredients/{id} | GET | Ingredient list |
| /ingredient/browse | GET | Search ingredients |

## Integration Examples

```bash
# Map to PubChem
python3 << 'EOF'
import json
import requests

with open('dsld_ingredients.csv') as f:
    # Extract unique ingredient names
    ingredients = set()
    for line in f:
        parts = line.strip().split(',')
        if len(parts) > 1:
            ingredients.add(parts[1])

# Save for cross-reference
with open('ingredient_names.txt', 'w') as f:
    for ing in sorted(ingredients):
        f.write(ing + '\n')

print(f"Unique ingredients: {len(ingredients)}")
EOF
```

## Related Resources

- [Schema Documentation](./schema.md)
- [ConsumerLab](../consumerlab/) - Product testing
- [Natural Medicines](../natural.medicines/) - Evidence database
- [NIH ODS](https://ods.od.nih.gov/) - Office of Dietary Supplements
