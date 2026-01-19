# Data Integration Guide: Traditional Medicine → Genes → Pathways

**Created:** January 2026
**Purpose:** Practical guide for downloading and integrating traditional medicine data with genetic/pathway information

---

## Executive Summary

This guide consolidates the research from our comprehensive database analysis to provide step-by-step instructions for:
1. **Downloading** traditional medicine data (TCM, Ayurveda, Kampo, Western Herbal)
2. **Linking** compounds and formulations to genes/proteins
3. **Mapping** to biochemical pathways

**Key Finding:** All the data needed exists in publicly accessible databases. The integration challenge is harmonizing identifiers across systems.

---

## Quick Reference: Best Databases by Use Case

| Need | Database | Why |
|------|----------|-----|
| TCM compound → target | **BATMAN-TCM 2.0** | 2.3M predicted interactions, REST API |
| Ayurveda compound → target | **IMPPAT 2.0** | 27,365 predicted interactions, CC BY 4.0 |
| Kampo formula → target | **KampoDB** | 62,906 proteins, CC BY-SA 4.0 |
| Western herbal → pathway | **Dr. Duke's + KEGG** | CC0 license, bulk CSV |
| All pathways | **Reactome API** | Free, comprehensive, REST API |
| Gene/protein lookup | **UniProt** | Gold standard identifiers |

---

## Part 1: Downloading Data

### 1.1 TCM Databases

#### BATMAN-TCM 2.0 (Recommended - Has API)

**URL:** http://bionet.ncpsb.org.cn/batman-tcm/

**Download Methods:**

```bash
# REST API - Get ingredients for a formula
curl "http://bionet.ncpsb.org.cn/batman-tcm/api/ingredient?formula=mahuang"

# REST API - Get predicted targets for an ingredient
curl "http://bionet.ncpsb.org.cn/batman-tcm/api/target?ingredient=ephedrine"

# Bulk Download (tab-delimited)
# Navigate to: http://bionet.ncpsb.org.cn/batman-tcm/download.php
# Files available:
#   - formula_ingredient.txt (54,832 formulas → 39,171 ingredients)
#   - ingredient_target_known.txt (17,068 known interactions)
#   - ingredient_target_predicted.txt (2,319,272 predicted interactions, AUC=0.97)
```

**Data Schema:**
```
ingredient_target_predicted.txt:
Column 1: Ingredient ID
Column 2: Ingredient Name
Column 3: Target Gene Symbol
Column 4: Target UniProt ID
Column 5: Prediction Score (0-1)
```

**License:** CC BY-NC 4.0 (non-commercial; contact for commercial)

---

#### TCMSP (Spider Tool Available)

**URL:** https://tcmsp-e.com/tcmsp.php

**Download Method:**
```bash
# Clone the spider tool
git clone https://github.com/shujuecn/TCMSP-Spider

# Install dependencies
pip install requests pandas openpyxl

# Run spider to extract data
python tcmsp_spider.py --herb "Ephedra sinica"
```

**Output:** Excel file with:
- Herb → Ingredients (29,384 total)
- Ingredients → Targets (3,311 targets)
- Targets → Diseases (837 diseases)
- ADME properties (OB%, DL, Caco-2, BBB, HL)

**License:** ODbL 1.0 (commercial OK with attribution)

---

#### TCMBank (Bulk Download)

**URL:** https://tcmbank.cn/

**Download:**
- Navigate to Download section
- Get: Herb-Ingredient-Target-Disease mapping files
- 9,192 herbs, 61,966 ingredients, 15,179 targets

**License:** CC BY 4.0

---

### 1.2 Ayurveda Databases

#### IMPPAT 2.0 (Best for Ayurveda)

**URL:** https://cb.imsc.res.in/imppat/

**Download Methods:**
```bash
# Web export - Tab-separated files
# Navigate to any search result → Export option

# Structure files (per compound)
# Available formats: SDF, MOL, MOL2, PDB, PDBQT

# GitHub - Analysis code
git clone https://github.com/asamallab/IMPPAT2
```

**Key Data:**
- 4,010 medicinal plants
- 17,967 phytochemicals
- 27,365 predicted compound-target interactions
- 5,042 human target proteins (from STITCH, score ≥700)

**Target Data Format:**
```
Phytochemical ID | Phytochemical Name | Target UniProt ID | Target Gene | STITCH Score
IMPPAT001234    | Curcumin           | P04637           | TP53        | 850
```

**License:** CC BY 4.0 (commercial OK)

---

#### NPACT (Cancer Targets)

**URL:** https://webs.iiitd.edu.in/raghava/npact/

**Coverage:** 1,574 plant compounds with validated cancer targets

**Download:** Web export, curated compound-target pairs

---

### 1.3 Kampo Databases

#### KampoDB (Primary Kampo Source)

**URL:** https://wakanmoview.inm.u-toyama.ac.jp/kampo/

**Data Available:**
- 298 Kampo formulas
- 180 crude drugs
- 3,002 natural compounds
- 62,906 proteins (docking predictions)
- 460 known + 1,369 predicted targets

**Access:** Web queries; limited bulk download
**License:** CC BY-SA 4.0

**Compound IDs:** Map to KNApSAcK IDs, which link to PubChem

---

#### STORK (Official Japanese Reference)

**URL:** http://mpdb.nibiohn.go.jp/stork/

**Coverage:** All 148 approved Kampo formulas in Japan
- Standardized compositions
- Crude drug ratios
- Package insert data (English)

---

### 1.4 Western Herbal Databases

#### DSLD - NIH Dietary Supplement Label Database (Best API)

**URL:** https://dsld.od.nih.gov
**API:** https://api.ods.od.nih.gov/dsld/v9/

```bash
# Search products by ingredient
curl "https://api.ods.od.nih.gov/dsld/v9/browse-products?ingredient=curcumin"

# Get full label data
curl "https://api.ods.od.nih.gov/dsld/v9/label/12345"

# Get ingredient groups
curl "https://api.ods.od.nih.gov/dsld/v9/ingredient-groups"
```

**Coverage:** 200,000+ supplement labels
**License:** CC0 (public domain)

---

#### Dr. Duke's Phytochemical Database (Bulk CSV)

**URL:** https://phytochem.nal.usda.gov
**Bulk Download:** https://data.nal.usda.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases

```bash
# Download bulk data
wget https://data.nal.usda.gov/system/files/Duke-Source-CSV.zip
unzip Duke-Source-CSV.zip

# Key tables:
# - plants.csv (plant names, families)
# - chemicals.csv (compounds, CAS numbers)
# - activities.csv (biological effects)
# - plant_chemicals.csv (plant → compound links)
# - chemical_activities.csv (compound → activity links)
```

**License:** CC0 (public domain)

---

#### Health Canada LNHPD (REST API)

**URL:** https://health-products.canada.ca/lnhpd-bdpsnh/
**API:** https://health-products.canada.ca/api/natural-licences/

```bash
# Get all medicinal ingredients
curl "https://health-products.canada.ca/api/natural-licences/medicinalingredient/?lang=en&type=json"

# Get product licenses
curl "https://health-products.canada.ca/api/natural-licences/productlicence/?lang=en&type=json"
```

**License:** Open Government License - Canada

---

## Part 2: Linking Compounds to Genes/Proteins

### 2.1 Identifier Mapping Strategy

The key challenge is mapping between different ID systems:

```
Compound IDs:        PubChem CID, CAS, InChIKey, SMILES
Protein IDs:         UniProt ID, Gene Symbol, Entrez Gene ID
Pathway IDs:         KEGG pathway, Reactome stable ID, GO terms
```

**Recommended Mapping Flow:**
```
Traditional Medicine DB → PubChem CID → UniProt Target → Reactome Pathway
```

### 2.2 PubChem as Central Hub

Use PubChem to normalize compound identifiers:

```bash
# Get compound info by name
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/curcumin/JSON"

# Get compound targets (BioAssay data)
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/969516/assaysummary/JSON"

# Get compound by InChIKey
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/VFLDPWHFBUODDF-FCXRPNKRSA-N/JSON"
```

### 2.3 UniProt for Protein/Gene Mapping

```bash
# Get protein info by UniProt ID
curl "https://rest.uniprot.org/uniprotkb/P04637.json"

# Map gene symbol to UniProt
curl "https://rest.uniprot.org/uniprotkb/search?query=gene:TP53+AND+organism_id:9606&format=json"

# Get pathway annotations
curl "https://rest.uniprot.org/uniprotkb/P04637.json?fields=xref_reactome,xref_kegg"
```

---

## Part 3: Mapping to Biochemical Pathways

### 3.1 Reactome API (Recommended)

**URL:** https://reactome.org/ContentService/

```bash
# Get pathways for a gene
curl "https://reactome.org/ContentService/data/pathways/low/UniProt/P04637"

# Get pathway details
curl "https://reactome.org/ContentService/data/query/R-HSA-69620"

# Get pathway diagram data
curl "https://reactome.org/ContentService/exporter/diagram/R-HSA-69620.json"

# Search pathways by keyword
curl "https://reactome.org/ContentService/search/query?query=methylation&species=Homo+sapiens"
```

**Key Endpoints:**
- `/data/pathways/low/UniProt/{id}` - Get pathways for a protein
- `/data/participant/{id}/pathways` - Get pathways containing entity
- `/interactors/static/protein/{id}` - Get protein interactions
- `/analysis/identifiers/projection` - Analyze gene/protein list

### 3.2 KEGG API

```bash
# Get pathway list
curl "https://rest.kegg.jp/list/pathway/hsa"

# Get genes in a pathway
curl "https://rest.kegg.jp/link/hsa/hsa00010"  # Glycolysis pathway

# Get pathway image
curl "https://rest.kegg.jp/get/hsa00010/image"

# Convert gene IDs
curl "https://rest.kegg.jp/conv/uniprot/hsa:7157"  # TP53 gene
```

### 3.3 WikiPathways (Community Pathways)

**URL:** https://www.wikipathways.org/

```bash
# Search pathways
curl "https://webservice.wikipathways.org/findPathwaysByText?query=methylation&format=json"

# Get pathway info
curl "https://webservice.wikipathways.org/getPathway?pwId=WP4220&format=json"

# Download pathway as GPML
curl "https://webservice.wikipathways.org/getPathwayAs?pwId=WP4220&fileType=gpml"
```

---

## Part 4: Complete Integration Pipeline

### 4.1 Example: TCM Herb → Pathway

```python
"""
Pipeline: Ephedra (Ma Huang) → Ephedrine → Targets → Pathways
"""

import requests
import json

# Step 1: Get ingredients from BATMAN-TCM
# (Using mock URL - actual endpoint may differ)
batman_url = "http://bionet.ncpsb.org.cn/batman-tcm/api"
herb_data = requests.get(f"{batman_url}/ingredient?formula=mahuang").json()

# Step 2: For each ingredient, get predicted targets
for ingredient in herb_data['ingredients']:
    targets = requests.get(f"{batman_url}/target?ingredient={ingredient['id']}").json()

    # Step 3: Map targets to UniProt
    for target in targets['targets']:
        uniprot_id = target.get('uniprot_id')
        if not uniprot_id:
            # Map gene symbol to UniProt
            gene = target['gene_symbol']
            uniprot_resp = requests.get(
                f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:9606&format=json"
            ).json()
            uniprot_id = uniprot_resp['results'][0]['primaryAccession']

        # Step 4: Get Reactome pathways
        pathways = requests.get(
            f"https://reactome.org/ContentService/data/pathways/low/UniProt/{uniprot_id}"
        ).json()

        print(f"Ingredient: {ingredient['name']}")
        print(f"  Target: {target['gene_symbol']} ({uniprot_id})")
        print(f"  Pathways: {[p['displayName'] for p in pathways[:3]]}")
```

### 4.2 Example: Ayurveda Compound → Pathway

```python
"""
Pipeline: Curcumin (from Turmeric) → Targets → Pathways
"""

# Step 1: Get curcumin data from IMPPAT (or PubChem)
pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
curcumin = requests.get(f"{pubchem_url}/compound/name/curcumin/JSON").json()
cid = curcumin['PC_Compounds'][0]['id']['id']['cid']

# Step 2: Get bioassay targets from PubChem
targets = requests.get(f"{pubchem_url}/compound/cid/{cid}/assaysummary/JSON").json()

# Step 3: Get known targets from IMPPAT (using STITCH scores)
# IMPPAT includes: 27,365 predicted interactions for 1,294 phytochemicals
# Filter by STITCH score >= 700 for high confidence

# Step 4: Map to pathways via Reactome
# (same as above)
```

---

## Part 5: Data Model for Integration

### Recommended Schema

```sql
-- Compounds from all sources
CREATE TABLE compounds (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255),
    pubchem_cid INTEGER,
    inchikey VARCHAR(27),
    smiles TEXT,
    source VARCHAR(50),  -- 'tcm', 'ayurveda', 'kampo', 'western'
    source_db VARCHAR(100),  -- 'batman-tcm', 'imppat', 'kampodb', etc.
    source_id VARCHAR(100)
);

-- Compound-Target relationships
CREATE TABLE compound_targets (
    id SERIAL PRIMARY KEY,
    compound_id INTEGER REFERENCES compounds(id),
    target_uniprot VARCHAR(20),
    target_gene_symbol VARCHAR(50),
    interaction_type VARCHAR(50),  -- 'known', 'predicted'
    confidence_score FLOAT,
    source_db VARCHAR(100),
    evidence_pmid INTEGER[]
);

-- Target-Pathway relationships
CREATE TABLE target_pathways (
    id SERIAL PRIMARY KEY,
    target_uniprot VARCHAR(20),
    pathway_id VARCHAR(50),  -- Reactome stable ID
    pathway_name VARCHAR(255),
    pathway_source VARCHAR(20)  -- 'reactome', 'kegg', 'wikipathways'
);

-- Formulations (TCM formulas, Ayurveda formulations, Kampo)
CREATE TABLE formulations (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255),
    tradition VARCHAR(50),  -- 'tcm', 'ayurveda', 'kampo'
    source_db VARCHAR(100),
    source_id VARCHAR(100)
);

-- Formulation-Compound relationships
CREATE TABLE formulation_compounds (
    id SERIAL PRIMARY KEY,
    formulation_id INTEGER REFERENCES formulations(id),
    compound_id INTEGER REFERENCES compounds(id),
    plant_source VARCHAR(255),
    plant_part VARCHAR(100)
);
```

---

## Part 6: Licensing Summary

| Database | License | Commercial Use |
|----------|---------|----------------|
| **BATMAN-TCM 2.0** | CC BY-NC 4.0 | Contact required |
| **TCMSP** | ODbL 1.0 | Yes, with attribution |
| **TCMBank** | CC BY 4.0 | Yes |
| **IMPPAT 2.0** | CC BY 4.0 | Yes |
| **KampoDB** | CC BY-SA 4.0 | Yes, share-alike |
| **DSLD** | CC0 | Yes (public domain) |
| **Dr. Duke's** | CC0 | Yes (public domain) |
| **Health Canada LNHPD** | Open Gov License | Yes |
| **Reactome** | CC BY 4.0 | Yes |
| **UniProt** | CC BY 4.0 | Yes |
| **PubChem** | Public Domain | Yes |

---

## Part 7: Recommended Integration Priority

### Phase 1: Core Pipeline (MVP)

1. **BATMAN-TCM 2.0** - TCM compound → target (API + bulk)
2. **IMPPAT 2.0** - Ayurveda compound → target (CC BY 4.0)
3. **Reactome API** - Target → pathway mapping
4. **UniProt API** - Protein/gene normalization

### Phase 2: Expanded Coverage

5. **KampoDB** - Kampo formula → target
6. **TCMSP** - TCM with ADME properties
7. **Dr. Duke's** - Western herbal bulk data
8. **DSLD API** - Supplement label data

### Phase 3: Deep Integration

9. **HERB 2.0** - Gene expression data
10. **WikiPathways** - Community pathways
11. **KEGG** - Metabolic pathways
12. **PharmGKB** - Pharmacogenomics

---

## Appendix: API Quick Reference

```bash
# BATMAN-TCM 2.0
http://bionet.ncpsb.org.cn/batman-tcm/api/

# PubChem
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{domain}/{search_type}/{search_term}/{output}

# UniProt
https://rest.uniprot.org/uniprotkb/search?query={query}&format=json

# Reactome
https://reactome.org/ContentService/data/pathways/low/UniProt/{uniprot_id}

# KEGG
https://rest.kegg.jp/{operation}/{argument}

# WikiPathways
https://webservice.wikipathways.org/{operation}?{params}&format=json

# DSLD (NIH Supplements)
https://api.ods.od.nih.gov/dsld/v9/{endpoint}

# Health Canada LNHPD
https://health-products.canada.ca/api/natural-licences/{endpoint}?lang=en&type=json
```

---

*Guide compiled from comprehensive database research, January 2026*
