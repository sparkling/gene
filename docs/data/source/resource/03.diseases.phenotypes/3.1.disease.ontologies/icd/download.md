---
id: download-icd
title: "International Classification of Diseases (ICD) Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# International Classification of Diseases (ICD) Download Instructions

## Quick Start

```bash
# Download ICD-10-CM (US Clinical Modification) - No registration
wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2024/icd10cm_order_2024.txt

# ICD-11 requires WHO API access (registration required)
```

## Prerequisites

- **wget** or **curl** for downloads
- **WHO ICD API token** for ICD-11 programmatic access
- **Registration** at id.who.int for ICD-11 API
- Approximately 500MB disk space for full ICD-10-CM

## Registration Requirements

| Version | Registration |
|---------|--------------|
| ICD-10-CM (CDC) | None required |
| ICD-10 (WHO) | WHO terms acceptance |
| ICD-11 (WHO) | API registration required |

## Download Methods

### Method 1: ICD-10-CM (CDC - US Version)

```bash
# Current year release
wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2024/icd10cm_order_2024.txt

# Tabular files
wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2024/icd10cm_tabular_2024.xml

# Code descriptions (long and short)
wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2024/icd10cm-CodesDescriptions-2024.zip

# Drug/Chemical index
wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2024/icd10cm_drug_2024.txt

# Guidelines
wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/ICD10CM/2024/icd10cm_guidelines_2024.pdf
```

### Method 2: ICD-11 WHO API

```bash
# 1. Register at https://icd.who.int/icdapi
# 2. Get client credentials

# Set credentials
CLIENT_ID="your_client_id"
CLIENT_SECRET="your_client_secret"

# Get OAuth token
TOKEN=$(curl -s -X POST "https://icdaccessmanagement.who.int/connect/token" \
  -d "client_id=${CLIENT_ID}" \
  -d "client_secret=${CLIENT_SECRET}" \
  -d "scope=icdapi_access" \
  -d "grant_type=client_credentials" | jq -r '.access_token')

# Search for diseases
curl "https://id.who.int/icd/entity/search?q=diabetes" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Accept: application/json" \
  -H "Accept-Language: en" \
  -o icd11_diabetes_search.json

# Get entity details
curl "https://id.who.int/icd/entity/1954798927" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Accept: application/json" \
  -o icd11_entity.json

# Get linearization (MMS)
curl "https://id.who.int/icd/release/11/2023-01/mms" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Accept: application/json" \
  -o icd11_mms_root.json
```

### Method 3: ICD-10 WHO Downloads

```bash
# WHO ICD-10 browser export (requires acceptance of terms)
# Visit: https://icd.who.int/browse10/2019/en

# ClaML format (XML)
# Available after registration at WHO classifications site
```

### Method 4: UMLS Integration

```bash
# ICD codes are available through UMLS (requires UMLS license)
# Download from: https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html

# After UMLS download, extract ICD:
unzip umls-2024AA-metathesaurus.zip
cd META
grep "ICD10CM" MRCONSO.RRF > icd10cm_umls.txt
grep "ICD11" MRCONSO.RRF > icd11_umls.txt
```

### Method 5: Alternative Sources

```bash
# CMS (Centers for Medicare & Medicaid)
wget https://www.cms.gov/files/zip/2024-icd-10-cm-codes-file.zip

# AAPC (American Academy of Professional Coders)
# Commercial resource - https://www.aapc.com/icd-10/codes/

# ICD.Codes (community resource)
# https://icd.codes/icd10cm - Web browsing only
```

## File Inventory

### ICD-10-CM (CDC)

| File | Size | Description |
|------|------|-------------|
| icd10cm_order_2024.txt | ~15 MB | Codes in tabular order |
| icd10cm_tabular_2024.xml | ~50 MB | Full hierarchy (XML) |
| icd10cm_codes_2024.txt | ~10 MB | All valid codes |
| icd10cm_drug_2024.txt | ~5 MB | Drug/chemical poisoning index |
| icd10cm_guidelines_2024.pdf | ~3 MB | Official coding guidelines |

### ICD-11 (WHO API)

| Endpoint | Description |
|----------|-------------|
| /icd/entity/ | Foundation entities |
| /icd/release/11/ | Linearization releases |
| /icd/entity/search | Full-text search |
| /icd/release/11/mms | MMS linearization |

## Post-Download Processing

```bash
# Parse ICD-10-CM order file
python3 << 'EOF'
import pandas as pd

# Read fixed-width format
# Columns: order_number, code, header_flag, short_desc, long_desc
df = pd.read_fwf('icd10cm_order_2024.txt',
                  colspecs=[(0,5), (6,13), (14,15), (16,76), (77,400)],
                  names=['order', 'code', 'header', 'short_desc', 'long_desc'])

# Remove header rows (header=0 means valid code)
valid_codes = df[df['header'] == '0']
print(f"Valid codes: {len(valid_codes)}")

# Save as TSV
valid_codes.to_csv('icd10cm_codes.tsv', sep='\t', index=False)

# Get diabetes codes (E08-E13)
diabetes = valid_codes[valid_codes['code'].str.startswith(('E08', 'E09', 'E10', 'E11', 'E12', 'E13'))]
print(f"Diabetes codes: {len(diabetes)}")
diabetes.to_csv('icd10cm_diabetes.tsv', sep='\t', index=False)
EOF

# Parse ICD-11 API response
python3 << 'EOF'
import json
import pandas as pd

with open('icd11_diabetes_search.json') as f:
    data = json.load(f)

results = []
for item in data.get('destinationEntities', []):
    results.append({
        'id': item.get('id'),
        'title': item.get('title'),
        'code': item.get('theCode'),
        'chapter': item.get('chapter')
    })

df = pd.DataFrame(results)
df.to_csv('icd11_search_results.tsv', sep='\t', index=False)
print(df)
EOF

# Build ICD hierarchy from XML
python3 << 'EOF'
import xml.etree.ElementTree as ET

tree = ET.parse('icd10cm_tabular_2024.xml')
root = tree.getroot()

# Extract all codes with parents
codes = []
for diag in root.iter('diag'):
    code = diag.find('name')
    desc = diag.find('desc')
    if code is not None and desc is not None:
        codes.append({
            'code': code.text,
            'description': desc.text
        })

print(f"Total codes: {len(codes)}")
EOF
```

## Verification

```bash
# Check ICD-10-CM file
head -20 icd10cm_order_2024.txt
wc -l icd10cm_order_2024.txt

# Verify code format
grep "^....E11" icd10cm_order_2024.txt | head -10

# Check XML structure
head -50 icd10cm_tabular_2024.xml

# Test ICD-11 API response
cat icd11_diabetes_search.json | jq '.destinationEntities | length'
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| ICD-11 v2025-01 | 2025-01-01 | API-based | Current |
| ICD-10-CM 2025 | 2024-10-01 | ~100 MB | Current |
| ICD-10-CM 2024 | 2023-10-01 | ~95 MB | Archived |

### Version Notes

ICD-11 (current WHO version):
- 17,000+ diagnostic categories
- Foundation layer with 100,000+ entities
- Multiple linearizations (MMS, Primary Care, etc.)
- REST API for programmatic access

ICD-10-CM 2025 (US Clinical Modification):
- 72,750+ valid diagnosis codes
- FY2025 effective October 1, 2024
- Annual updates from CDC/CMS

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://id.who.int/icd` (ICD-11) |
| Rate Limit | Token-based access |
| Auth Required | Yes (OAuth2 for ICD-11 API) |
| Documentation | https://icd.who.int/icdapi |

## Update Schedule

| Version | Frequency |
|---------|-----------|
| ICD-10-CM | Annual (October 1) |
| ICD-11 | Continuous updates |
| WHO ICD-10 | Periodic revisions |

## Common Issues

- **API rate limits**: ICD-11 API has rate limits; implement backoff
- **Code format**: ICD-10-CM codes have decimals (E11.9), ICD-11 uses different format
- **Header vs valid**: In order file, header=1 means category header, not billable
- **Language**: ICD-11 API supports multiple languages via Accept-Language header
- **Token expiry**: OAuth tokens expire; refresh periodically

## ICD-10 Chapter Structure

| Chapter | Range | Description |
|---------|-------|-------------|
| I | A00-B99 | Infectious diseases |
| II | C00-D49 | Neoplasms |
| III | D50-D89 | Blood diseases |
| IV | E00-E89 | Endocrine/metabolic |
| V | F01-F99 | Mental disorders |
| VI | G00-G99 | Nervous system |
| IX | I00-I99 | Circulatory system |
| X | J00-J99 | Respiratory system |

## Cross-Reference Mappings

```bash
# UMLS provides ICD mappings to other terminologies
# Extract ICD-SNOMED mappings from UMLS
grep "ICD10CM" MRMAP.RRF | grep "SNOMEDCT" > icd_snomed_map.txt

# Map ICD to MONDO
python3 << 'EOF'
import pronto

# Requires MONDO ontology
mondo = pronto.Ontology("mondo.obo")

# Find ICD cross-references
for term in mondo.terms():
    for xref in term.xrefs:
        if xref.id.startswith("ICD"):
            print(f"{term.id}\t{term.name}\t{xref.id}")
EOF
```

## Related Resources

- [MeSH](../mesh/) - Medical Subject Headings
- [MONDO](../mondo/) - Disease ontology with ICD mappings
- [SNOMED CT](../../../../05.standards.ontologies/5.1.core.ontologies/snomed.ct/) - Clinical terminology
