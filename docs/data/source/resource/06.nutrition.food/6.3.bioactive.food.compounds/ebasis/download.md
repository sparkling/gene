---
id: download-ebasis
title: "eBASIS Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
access_type: restricted
---

# eBASIS Download Instructions

## Access Classification

| Aspect | Status |
|--------|--------|
| Access Type | **Restricted - Registration + Data Use Agreement** |
| Registration | Required (Free) |
| Bulk Downloads | Requires EuroFIR Agreement |
| API Access | Not Public |
| Data Export | Limited (Individual queries) |
| Academic Use | Free with agreement |
| Commercial Use | License Required |

---

## Important Notice

eBASIS requires registration and a data use agreement for access. Bulk downloads are restricted and require contacting EuroFIR directly.

## Access Requirements

| Requirement | Details |
|-------------|---------|
| Registration | Required (free) |
| Data Agreement | EuroFIR data use agreement |
| Bulk Access | Contact EuroFIR |

## Registration

1. Visit: https://ebasis.eurofir.org/
2. Click "Register"
3. Complete registration form
4. Agree to terms of use
5. Confirm email

## Download Methods

### Method 1: Web Interface (Individual Records)

```
1. Log in to eBASIS
2. Navigate to Browse > Foods or Components
3. Select items of interest
4. View composition data
5. Export individual records (if permitted)
```

### Method 2: Query Builder

```
1. Log in to eBASIS
2. Use Query function
3. Build search criteria:
   - Food group: Fruits
   - Component group: Flavonoids
   - Quality score: >= 3
4. View/export results
```

### Method 3: Contact EuroFIR for Bulk Data

```
For research requiring bulk data:

Email: info@eurofir.org
Subject: eBASIS Data Access Request

Include:
- Institution name
- Research purpose
- Specific data requirements
- Intended use
- Publications planned
```

### Method 4: API Access (Limited)

eBASIS does not provide public API access. For programmatic access, contact EuroFIR about partnership options.

## Data Available Through Web Interface

| Data Type | Access Method |
|-----------|---------------|
| Food composition | Browse/Query |
| Component information | Browse |
| References | View per record |
| Quality scores | Included |

## Example Manual Data Collection

```bash
# For personal research notes (comply with ToS)

# 1. Create research notes structure
mkdir -p ebasis_notes/{foods,components,composition}

# 2. Document findings
cat > ebasis_notes/composition/apple_flavonoids.md << 'EOF'
# Apple Flavonoid Composition
Source: eBASIS (EuroFIR)
Access Date: 2026-01-23
Registration Required: Yes

## Food Information
- Food: Apple, raw, with skin
- Scientific name: Malus domestica
- LanguaL code: A0148

## Flavonoid Content (mg/100g edible portion)
| Component | Mean | Min | Max | Quality |
|-----------|------|-----|-----|---------|
| Quercetin | 4.01 | 2.10 | 6.45 | 4 |
| Catechin | 8.50 | 5.20 | 12.30 | 3 |
| Epicatechin | 7.20 | 4.10 | 10.50 | 4 |

## Data Quality
- Quality scores based on EuroFIR criteria
- Scores range from 1 (poor) to 5 (excellent)

## Citation
eBASIS - EuroFIR. https://ebasis.eurofir.org/
EOF
```

## Alternative Data Sources

For openly accessible bioactive compound data, consider:

| Source | Access | Coverage |
|--------|--------|----------|
| [FooDB](../../6.1.food.composition/foodb/) | Free | 70,000+ compounds |
| [Phenol-Explorer](../phenol.explorer/) | Free | Polyphenols |
| [PhytoHub](../phytohub/) | Free | Phytochemicals |
| [USDA FoodData](../../6.1.food.composition/usda.fooddata/) | Free | Flavonoids subset |

## Phenol-Explorer Data (Open Access Alternative)

```bash
# Phenol-Explorer provides similar data openly
# Download polyphenol composition data

wget http://phenol-explorer.eu/downloads/composition_data.csv

# Contains:
# - Food items
# - Polyphenol compounds
# - Concentrations (mg/100g)
# - References
```

## Data Use Agreement Terms

Typical EuroFIR agreement terms include:
- Data for stated research purpose only
- Attribution required in publications
- No redistribution without permission
- Report publications to EuroFIR
- Acknowledge funding sources

## EuroFIR Contact Information

| Type | Contact |
|------|---------|
| General inquiries | info@eurofir.org |
| Data access | data@eurofir.org |
| Website | https://www.eurofir.org/ |

## Integration with Other Databases

```bash
# Cross-reference eBASIS compounds with PubChem
# Using CAS numbers from manual records

python3 << 'EOF'
import requests

# Example: Quercetin (CAS: 117-39-5)
cas = "117-39-5"

# Get PubChem compound
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/quercetin/JSON"
response = requests.get(url)

if response.ok:
    data = response.json()
    cid = data['PC_Compounds'][0]['id']['id']['cid']
    print(f"PubChem CID for Quercetin: {cid}")
EOF
```

## Update Schedule

| Content | Frequency |
|---------|-----------|
| New composition data | As literature published |
| Quality re-evaluation | Periodic |
| Database version | Major updates annually |

## Common Questions

### Is bulk download available?
Not publicly. Contact EuroFIR for research data access agreements.

### What is the data quality score?
EuroFIR quality rating from 1-5 based on documentation, sampling, and analytical methods.

### Can I use data in publications?
Yes, with proper citation and adherence to data use agreement terms.

### How does eBASIS differ from Phenol-Explorer?
eBASIS covers broader bioactive categories with EuroFIR quality evaluation; Phenol-Explorer focuses on polyphenols with open access.

---

## Data Structure Overview (For Registered Users)

### Available Data Categories

| Category | Content | Export Options |
|----------|---------|----------------|
| Component Data | 300+ bioactive compounds | Browse/Query |
| Food Composition | 500+ foods, 25,000+ data points | Browse/Query |
| Quality Scores | EuroFIR quality evaluation (1-5) | Included with data |
| References | Primary literature citations | View per record |

### Data Fields Available Through Interface

**Component Fields:**
- Component name (common and systematic)
- Chemical identifiers (CAS, PubChem CID, InChIKey)
- Component group/subgroup classification
- Molecular properties (formula, weight)
- Related compounds

**Composition Data Fields:**
- Food name and scientific name
- Component name and group
- Concentration value with units (mg/100g)
- Value type (mean, median, single)
- Min/max range and standard deviation
- Number of samples
- Data quality score (1-5)
- Analytical method used
- Primary reference

**Quality Score Breakdown:**
- Food description quality
- Sampling procedure quality
- Analytical method quality
- Analytical QC quality
- Number of samples adequacy
- Component identification quality

---

## Programmatic Alternatives

For researchers needing programmatic access to bioactive compound data, these open alternatives provide similar information:

### Phenol-Explorer (Recommended Alternative)
```bash
# Free, open access polyphenol database
# Direct download available

# Download composition data
wget http://phenol-explorer.eu/downloads/composition_data.csv

# Download component data
wget http://phenol-explorer.eu/downloads/compounds.csv

# Download food data
wget http://phenol-explorer.eu/downloads/foods.csv

# Data includes:
# - 500+ polyphenols
# - 400+ foods
# - Chromatographic and after hydrolysis values
# - Full references
```

### FooDB (Comprehensive Alternative)
```bash
# Free food compound database with API

# Search compounds in foods
curl "https://foodb.ca/compounds.json?q=quercetin"

# Search foods containing compound
curl "https://foodb.ca/foods.json?q=apple"

# Full database download available
wget https://foodb.ca/public/system/downloads/foodb_2020_04_07_csv.tar.gz

# Includes:
# - 70,000+ compounds
# - 800+ foods
# - Concentrations and references
```

### USDA FoodData Central (Flavonoids)
```bash
# Free USDA flavonoid data

# API access
curl "https://api.nal.usda.gov/fdc/v1/foods/search?query=apple&dataType=SR%20Legacy&api_key=DEMO_KEY"

# Bulk download
wget https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_sr_legacy_food_csv_2021-10-28.zip

# Includes limited flavonoid data for common foods
```

### PhytoHub
```bash
# Free phytochemical database
# Web interface: http://phytohub.eu/

# Provides:
# - Dietary phytochemical data
# - Metabolite information
# - Food sources
# No bulk download, but browsable interface
```

---

## Comparison of eBASIS vs Open Alternatives

| Feature | eBASIS | Phenol-Explorer | FooDB | USDA |
|---------|--------|-----------------|-------|------|
| Bioactive compounds | 300+ | 500+ (polyphenols only) | 70,000+ | Limited |
| Foods covered | 500+ | 400+ | 800+ | 8,000+ |
| Quality scores | Yes (EuroFIR) | No | No | No |
| Bulk download | Agreement required | Yes (free) | Yes (free) | Yes (free) |
| API access | No | No | Yes | Yes |
| Data licensing | EuroFIR agreement | CC BY | CC BY | Public domain |

---

## EuroFIR Data Request Process

For bulk access to eBASIS data:

1. **Prepare Request Documentation:**
   - Institution name and type
   - Principal investigator details
   - Research project description
   - Specific data requirements
   - Intended use and publications
   - Data security measures

2. **Contact EuroFIR:**
   ```
   Email: data@eurofir.org
   Subject: eBASIS Data Access Request - [Institution]

   Include:
   - Completed data request form
   - Research protocol summary
   - IRB/Ethics approval (if applicable)
   ```

3. **Review Process:**
   - EuroFIR reviews request (2-4 weeks)
   - Data use agreement sent if approved
   - Sign and return agreement
   - Receive data access credentials or files

4. **Agreement Terms:**
   - Data for stated purpose only
   - No redistribution
   - Attribution in publications
   - Report publications to EuroFIR
   - Data destruction upon project completion

---

## Related Resources

- [Schema Documentation](./schema.md)
- [FooDB](../../6.1.food.composition/foodb/) - Open food compound database
- [Phenol-Explorer](../phenol.explorer/) - Open polyphenol database
