---
id: download-exposome-explorer
title: "Exposome-Explorer Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Exposome-Explorer Download Instructions

## Quick Start

```bash
# Download biomarker data
wget http://exposome-explorer.iarc.fr/downloads/biomarkers.tsv

# Download exposure associations
wget http://exposome-explorer.iarc.fr/downloads/associations.tsv
```

## Prerequisites

- **wget** or **curl** for downloads
- Minimal disk space (~50 MB)

## No Registration Required

Exposome-Explorer data is freely available under Creative Commons Attribution 4.0 license.

## Download Methods

### Method 1: Direct TSV Downloads

```bash
# Download all available data files
mkdir exposome_explorer && cd exposome_explorer

# Biomarkers
wget http://exposome-explorer.iarc.fr/downloads/biomarkers.tsv

# Exposures
wget http://exposome-explorer.iarc.fr/downloads/exposures.tsv

# Associations
wget http://exposome-explorer.iarc.fr/downloads/associations.tsv

# Concentrations
wget http://exposome-explorer.iarc.fr/downloads/concentrations.tsv

# References
wget http://exposome-explorer.iarc.fr/downloads/references.tsv
```

### Method 2: Web Interface Export

```
1. Navigate to: http://exposome-explorer.iarc.fr/
2. Browse or search data
3. Use export buttons where available
4. Download as TSV/CSV
```

### Method 3: Programmatic Access

```bash
# Python script to download and process
python3 << 'EOF'
import pandas as pd
import requests
import os

BASE_URL = "http://exposome-explorer.iarc.fr/downloads"

files = [
    "biomarkers.tsv",
    "exposures.tsv",
    "associations.tsv",
    "concentrations.tsv",
    "references.tsv"
]

os.makedirs("exposome_data", exist_ok=True)

for filename in files:
    url = f"{BASE_URL}/{filename}"
    print(f"Downloading {filename}...")

    response = requests.get(url)
    if response.ok:
        filepath = f"exposome_data/{filename}"
        with open(filepath, 'wb') as f:
            f.write(response.content)
        print(f"  Saved to {filepath}")
    else:
        print(f"  Failed: {response.status_code}")
EOF
```

### Method 4: Web Scraping (Supplementary)

```bash
# For data not in bulk downloads
python3 << 'EOF'
import requests
from bs4 import BeautifulSoup
import time

BASE_URL = "http://exposome-explorer.iarc.fr"

# Get biomarker list
response = requests.get(f"{BASE_URL}/biomarkers")
soup = BeautifulSoup(response.content, 'html.parser')

biomarkers = []
for link in soup.find_all('a', href=True):
    if '/biomarker/' in link['href']:
        biomarkers.append({
            'name': link.text.strip(),
            'url': BASE_URL + link['href']
        })

print(f"Found {len(biomarkers)} biomarkers")

# Respectful scraping with delays
time.sleep(1)
EOF
```

## File Inventory

### Available Downloads

| File | Size | Description |
|------|------|-------------|
| biomarkers.tsv | ~1 MB | Biomarker definitions |
| exposures.tsv | ~500 KB | Exposure definitions |
| associations.tsv | ~5 MB | Biomarker-exposure associations |
| concentrations.tsv | ~2 MB | Reference concentrations |
| references.tsv | ~1 MB | Publication references |

### File Structure

#### biomarkers.tsv
| Column | Description |
|--------|-------------|
| biomarker_id | Unique identifier |
| name | Biomarker name |
| pubchem_cid | PubChem compound ID |
| hmdb_id | HMDB identifier |
| cas_number | CAS registry number |
| type | Biomarker category |

#### associations.tsv
| Column | Description |
|--------|-------------|
| association_id | Unique identifier |
| biomarker_id | Foreign key to biomarkers |
| exposure_id | Foreign key to exposures |
| specimen_type | Sample type |
| effect_size | Association strength |
| p_value | Statistical significance |
| sample_size | Study size |
| reference_id | Publication reference |

## Post-Download Processing

```bash
# Load and analyze data
python3 << 'EOF'
import pandas as pd

# Load datasets
biomarkers = pd.read_csv('exposome_data/biomarkers.tsv', sep='\t')
exposures = pd.read_csv('exposome_data/exposures.tsv', sep='\t')
associations = pd.read_csv('exposome_data/associations.tsv', sep='\t')

print(f"Biomarkers: {len(biomarkers)}")
print(f"Exposures: {len(exposures)}")
print(f"Associations: {len(associations)}")

# Find dietary biomarkers
dietary = exposures[exposures['category'] == 'Dietary']
print(f"Dietary exposures: {len(dietary)}")

# Top biomarkers by number of associations
top_biomarkers = associations.groupby('biomarker_id').size().sort_values(ascending=False).head(10)
print("\nTop biomarkers by associations:")
print(top_biomarkers)
EOF

# Create SQLite database
python3 << 'EOF'
import pandas as pd
import sqlite3

conn = sqlite3.connect('exposome_explorer.db')

# Load and import data
for table in ['biomarkers', 'exposures', 'associations', 'concentrations', 'references']:
    try:
        df = pd.read_csv(f'exposome_data/{table}.tsv', sep='\t')
        df.to_sql(table, conn, if_exists='replace', index=False)
        print(f"Imported {table}: {len(df)} rows")
    except Exception as e:
        print(f"Error with {table}: {e}")

# Create indexes
conn.execute('CREATE INDEX IF NOT EXISTS idx_bio_id ON associations(biomarker_id)')
conn.execute('CREATE INDEX IF NOT EXISTS idx_exp_id ON associations(exposure_id)')
conn.commit()
conn.close()

print("\nSQLite database created: exposome_explorer.db")
EOF

# Query specific biomarker
python3 << 'EOF'
import sqlite3

conn = sqlite3.connect('exposome_explorer.db')

# Find coffee biomarkers
query = """
SELECT DISTINCT b.name as biomarker, e.name as exposure,
       a.specimen_type, a.effect_size
FROM associations a
JOIN biomarkers b ON a.biomarker_id = b.biomarker_id
JOIN exposures e ON a.exposure_id = e.exposure_id
WHERE e.name LIKE '%coffee%'
ORDER BY a.effect_size DESC
"""

results = conn.execute(query).fetchall()
print("Coffee biomarkers:")
for row in results[:10]:
    print(f"  {row[0]}: r={row[3]:.3f} ({row[2]})")

conn.close()
EOF
```

## Verification

```bash
# Check file integrity
wc -l exposome_data/*.tsv

# Check structure
head -5 exposome_data/biomarkers.tsv

# Count associations per exposure type
cut -f4 exposome_data/exposures.tsv | sort | uniq -c | sort -rn

# Verify biomarker-exposure links
awk -F'\t' 'NR>1 {print $2}' exposome_data/associations.tsv | sort -u | wc -l
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | As literature published |
| Major revisions | Annually |

## Common Issues

- **TSV format**: Tab-separated, not comma-separated
- **Encoding**: UTF-8
- **Missing values**: Empty strings or "NA"
- **ID cross-references**: Some external IDs may be outdated

## Integration Examples

```bash
# Map biomarkers to HMDB
python3 << 'EOF'
import pandas as pd

biomarkers = pd.read_csv('exposome_data/biomarkers.tsv', sep='\t')

# Extract HMDB mappings
hmdb_map = biomarkers[biomarkers['hmdb_id'].notna()][['biomarker_id', 'name', 'hmdb_id']]
hmdb_map.to_csv('exposome_hmdb_map.tsv', sep='\t', index=False)

print(f"HMDB mappings: {len(hmdb_map)}")
EOF

# Map to PubChem
python3 << 'EOF'
import pandas as pd

biomarkers = pd.read_csv('exposome_data/biomarkers.tsv', sep='\t')

# Extract PubChem mappings
pubchem_map = biomarkers[biomarkers['pubchem_cid'].notna()][['biomarker_id', 'name', 'pubchem_cid']]
pubchem_map.to_csv('exposome_pubchem_map.tsv', sep='\t', index=False)

print(f"PubChem mappings: {len(pubchem_map)}")
EOF
```

## Exposure Category Analysis

```bash
python3 << 'EOF'
import pandas as pd

exposures = pd.read_csv('exposome_data/exposures.tsv', sep='\t')
associations = pd.read_csv('exposome_data/associations.tsv', sep='\t')

# Count associations per category
merged = associations.merge(exposures, on='exposure_id')
category_counts = merged.groupby('category').size().sort_values(ascending=False)

print("Associations by exposure category:")
for cat, count in category_counts.items():
    print(f"  {cat}: {count}")
EOF
```

## Related Resources

- [Schema Documentation](./schema.md)
- [HMDB](../hmdb/) - Human Metabolome Database
- [FooDB](../../6.1.food.composition/foodb/) - Food compound database
