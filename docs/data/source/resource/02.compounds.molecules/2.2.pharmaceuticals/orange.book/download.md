---
id: download-orange.book
title: "FDA Orange Book Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# FDA Orange Book Download Instructions

## Quick Start

```bash
# Download all Orange Book data files
wget https://www.fda.gov/media/76860/download -O orangebook.zip
unzip orangebook.zip
```

## Prerequisites

- **wget** or **curl** for downloads
- **unzip** for extraction
- ~50 MB disk space
- Text processing tools (awk, sed) or database software

## No Registration Required

Data is public domain (US Government) and freely downloadable.

## Download Methods

### Method 1: FDA Data Files (Recommended)

```bash
# Download ZIP archive with all data files
wget https://www.fda.gov/media/76860/download -O orangebook.zip
unzip orangebook.zip

# Archive contains:
# - products.txt (approved products)
# - patent.txt (patent listings)
# - exclusivity.txt (exclusivity data)
```

### Method 2: Individual Files via FDA Website

```bash
# Navigate to:
# https://www.fda.gov/drugs/drug-approvals-and-databases/orange-book-data-files

# Download individual files:
# - Approved Drug Products
# - Patent and Exclusivity Data
# - Discontinued Drug Products
```

### Method 3: openFDA API

```bash
# Search products via openFDA
curl "https://api.fda.gov/drug/drugsfda.json?search=products.brand_name:LIPITOR" | jq

# Get application information
curl "https://api.fda.gov/drug/drugsfda.json?search=application_number:NDA020702" | jq

# Search by active ingredient
curl "https://api.fda.gov/drug/drugsfda.json?search=products.active_ingredients.name:ATORVASTATIN" | jq
```

### Method 4: Web Search Interface

```bash
# Interactive search at:
# https://www.accessdata.fda.gov/scripts/cder/ob/

# Search by:
# - Active Ingredient
# - Proprietary Name
# - Patent Number
# - Applicant
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| products.txt | ~15 MB | Approved drug products |
| patent.txt | ~8 MB | Patent listings |
| exclusivity.txt | ~2 MB | Exclusivity data |
| discontinued.txt | ~5 MB | Discontinued products |

## File Format

Files are tilde (~) delimited text:

```
# products.txt header
Ingredient~DF;Route~Trade_Name~Applicant~Strength~Appl_Type~Appl_No~Product_No~TE_Code~Approval_Date~RLD~RS~Type

# Example record
ACETAMINOPHEN~TABLET;ORAL~TYLENOL~JOHNSON & JOHNSON~325MG~N~020032~001~AA~Approved Prior to Jan 1, 1982~Yes~No~RX
```

## Post-Download Processing

```bash
# Preview file structure
head -5 products.txt

# Count products
wc -l products.txt

# Convert to CSV
sed 's/~/,/g' products.txt > products.csv

# Extract unique active ingredients
cut -d'~' -f1 products.txt | sort | uniq > ingredients.txt

# Filter RLD (Reference Listed Drugs)
awk -F'~' '$11=="Yes"' products.txt > rld_products.txt

# Find products with AB rating
awk -F'~' '$9=="AB"' products.txt > ab_rated.txt

# Load into SQLite
sqlite3 orangebook.db << 'EOF'
CREATE TABLE products (
  ingredient TEXT,
  df_route TEXT,
  trade_name TEXT,
  applicant TEXT,
  strength TEXT,
  appl_type TEXT,
  appl_no TEXT,
  product_no TEXT,
  te_code TEXT,
  approval_date TEXT,
  rld TEXT,
  rs TEXT,
  type TEXT
);
.mode csv
.separator ~
.import products.txt products
DELETE FROM products WHERE ingredient='Ingredient';
EOF

# Query for specific drugs
sqlite3 orangebook.db "SELECT * FROM products WHERE trade_name LIKE '%LIPITOR%';"
```

## Patent Analysis

```bash
# Patents expiring this year
YEAR=$(date +%Y)
awk -F'~' -v year="$YEAR" '$4 ~ year' patent.txt

# Count patents per application
cut -d'~' -f2 patent.txt | sort | uniq -c | sort -rn | head -20

# Find substance vs product patents
awk -F'~' '$5=="Y"' patent.txt | wc -l  # Drug substance
awk -F'~' '$6=="Y"' patent.txt | wc -l  # Drug product
```

## Exclusivity Analysis

```bash
# NCE exclusivities
awk -F'~' '$1=="NCE"' exclusivity.txt

# Exclusivities expiring soon
awk -F'~' '{print $3}' exclusivity.txt | sort | head -20

# Count by exclusivity type
cut -d'~' -f1 exclusivity.txt | sort | uniq -c
```

## Verification

```bash
# Check file integrity
file products.txt
# Should be ASCII text

# Validate field count
head -1 products.txt | tr '~' '\n' | wc -l

# Check for required fields
awk -F'~' 'NF!=13 {print NR": "NF" fields"}' products.txt | head
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Data files | Monthly |
| Web interface | Daily |
| openFDA API | Near real-time |

## Common Issues

- **Delimiter handling**: Fields may contain special characters; use proper parsing
- **Date formats**: Mix of formats (some "Approved Prior to..." text)
- **Discontinued products**: In separate file; may need merging
- **TE code changes**: Track version history for changes

## Integration with Other Sources

```bash
# Cross-reference with DailyMed via NDC
# Cross-reference with DrugBank via application number
# Cross-reference with RxNorm via drug name

# Example: Find DailyMed label for Orange Book product
DRUG="LIPITOR"
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name=$DRUG" | jq
```

## License

Public Domain (US Government) - Free for any use
