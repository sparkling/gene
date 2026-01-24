---
id: download-drugbank
title: "DrugBank Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# DrugBank Download Instructions

## Quick Start

```bash
# After registration, download from releases page
curl -L -o drugbank_all_full_database.xml.zip \
  "https://go.drugbank.com/releases/latest/downloads/all-full-database"
```

## Prerequisites

- **Free academic registration** at DrugBank (required)
- **unzip** for extracting archives
- **XML parser** (xmlstarlet, lxml, or similar) for processing
- Approximately 5GB disk space for full database

## Registration Process

### Step 1: Create Account

1. Navigate to https://go.drugbank.com/public_users/sign_up
2. Select "Academic/Non-profit" account type
3. Provide institutional email address
4. Verify email and complete profile

### Step 2: Request Download Access

1. Log in to DrugBank
2. Navigate to Releases: https://go.drugbank.com/releases
3. Select latest release
4. Accept CC BY-NC 4.0 license terms
5. Download link becomes available

## Download Methods

### Method 1: Web Interface (Recommended)

1. Log in at https://go.drugbank.com
2. Navigate to https://go.drugbank.com/releases/latest
3. Download desired files directly

### Method 2: Command Line (with authentication)

```bash
# Set credentials (after registration)
DRUGBANK_EMAIL="your@email.com"
DRUGBANK_TOKEN="your_download_token"

# Download full database XML
curl -L -u "${DRUGBANK_EMAIL}:${DRUGBANK_TOKEN}" \
  -o drugbank_all_full_database.xml.zip \
  "https://go.drugbank.com/releases/5-1-12/downloads/all-full-database"

# Download structures (SDF)
curl -L -u "${DRUGBANK_EMAIL}:${DRUGBANK_TOKEN}" \
  -o drugbank_all_structures.sdf.zip \
  "https://go.drugbank.com/releases/5-1-12/downloads/all-structures"
```

### Method 3: Specific Data Types

```bash
# Approved drugs only
curl -L -o approved_drugs.xml.zip \
  "https://go.drugbank.com/releases/latest/downloads/approved-full-database"

# Drug-drug interactions
curl -L -o drug_interactions.csv.zip \
  "https://go.drugbank.com/releases/latest/downloads/drug-drug-interactions"

# Drug targets
curl -L -o drug_targets.csv.zip \
  "https://go.drugbank.com/releases/latest/downloads/target-all-polypeptide-ids"
```

## File Inventory

### Core Database Files

| File | Size | Description |
|------|------|-------------|
| drugbank_all_full_database.xml.zip | ~300 MB | Complete database (XML) |
| drugbank_all_structures.sdf.zip | ~150 MB | Chemical structures (SDF) |
| drugbank_all_drug_links.csv.zip | ~5 MB | External identifiers |
| drugbank_vocabulary.csv.zip | ~10 MB | Drug synonyms |

### Subset Files

| File | Size | Description |
|------|------|-------------|
| drugbank_approved*.xml.zip | ~100 MB | Approved drugs only |
| drugbank_experimental*.xml.zip | ~150 MB | Experimental drugs |
| drugbank_nutraceuticals*.xml.zip | ~20 MB | Nutraceuticals |

### Interaction Data

| File | Size | Description |
|------|------|-------------|
| drug-drug-interactions.csv | ~50 MB | DDI data |
| drug-food-interactions.csv | ~2 MB | Food interactions |

### Target Data

| File | Size | Description |
|------|------|-------------|
| target-all-polypeptide-ids.csv | ~5 MB | All targets with IDs |
| enzyme-all-polypeptide-ids.csv | ~2 MB | Drug-metabolizing enzymes |
| transporter-all-polypeptide-ids.csv | ~1 MB | Drug transporters |
| carrier-all-polypeptide-ids.csv | ~500 KB | Drug carriers |

## Post-Download Processing

```bash
# Extract XML database
unzip drugbank_all_full_database.xml.zip
# Creates: full database.xml (~1.5 GB uncompressed)

# Rename for easier handling
mv "full database.xml" drugbank_full.xml

# Extract specific fields with xmlstarlet
xmlstarlet sel -t -m "//drug" \
  -v "drugbank-id[@primary='true']" -o $'\t' \
  -v "name" -o $'\t' \
  -v "description" -n \
  drugbank_full.xml > drugs_basic.tsv

# Parse with Python
python3 << 'EOF'
import xml.etree.ElementTree as ET

tree = ET.parse('drugbank_full.xml')
root = tree.getroot()

ns = {'db': 'http://www.drugbank.ca'}

for drug in root.findall('db:drug', ns):
    drug_id = drug.find('db:drugbank-id[@primary="true"]', ns).text
    name = drug.find('db:name', ns).text
    print(f"{drug_id}\t{name}")
EOF
```

## Verification

```bash
# Verify zip integrity
unzip -t drugbank_all_full_database.xml.zip

# Check XML structure
head -100 drugbank_full.xml

# Count drug entries
grep -c "<drug type=" drugbank_full.xml

# Validate XML
xmllint --noout drugbank_full.xml
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major releases | Quarterly |
| Minor updates | Monthly |
| Interaction updates | Weekly |

## Common Issues

- **Access denied**: Ensure academic registration is complete and license accepted
- **Large XML parsing**: Use streaming parsers (SAX, iterparse) for memory efficiency
- **Namespace issues**: DrugBank XML uses namespace `http://www.drugbank.ca`
- **Character encoding**: Files are UTF-8 encoded; ensure proper handling
- **Version compatibility**: Structure may change between versions; check release notes

## API Access (Commercial)

For programmatic access without download:

```bash
# DrugBank API (requires commercial subscription)
curl -X GET "https://api.drugbank.com/v1/drugs/DB00945" \
  -H "Authorization: Bearer YOUR_API_KEY"
```

## License Restrictions

| Use Case | Permitted |
|----------|-----------|
| Academic research | Yes (CC BY-NC 4.0) |
| Non-commercial applications | Yes |
| Commercial use | Requires license |
| Redistribution | Attribution required |

## Related Downloads

- **DrugBank Structures**: SDF format for cheminformatics
- **DrugBank External Links**: Mapping to PubChem, ChEMBL, etc.
- **DrugBank Sequences**: Target protein sequences (FASTA)
