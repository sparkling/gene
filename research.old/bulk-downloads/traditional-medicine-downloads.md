# Traditional Medicine Database Bulk Downloads

This document provides comprehensive information on bulk download methods for Traditional Chinese Medicine (TCM), Ayurveda, and related natural products databases.

---

## Table of Contents

1. [BATMAN-TCM 2.0](#1-batman-tcm-20)
2. [TCMSP](#2-tcmsp)
3. [TCMBank](#3-tcmbank)
4. [HERB 2.0](#4-herb-20)
5. [IMPPAT 2.0](#5-imppat-20)
6. [KampoDB](#6-kampodb)
7. [Dr. Duke's Phytochemical Database](#7-dr-dukes-phytochemical-database)
8. [LOTUS](#8-lotus)
9. [KNApSAcK](#9-knapsack)
10. [Processing Recommendations](#10-processing-recommendations)

---

## 1. BATMAN-TCM 2.0

**Bioinformatics Analysis Tool for Molecular mechANism of Traditional Chinese Medicine**

### Database URL
- **Main Site**: http://bionet.ncpsb.org.cn/batman-tcm/
- **Download Page**: Navigate to "Download" section from main site

### Database Contents
| Data Type | Count |
|-----------|-------|
| Known TTIs (TCM ingredient-Target protein Interactions) | 17,068 |
| Predicted TTIs | ~2,319,272 |
| Total TTIs | ~2.3 million |

### Bulk Download Files

| File | Format | Description |
|------|--------|-------------|
| All TTI Data | Tab-delimited text (.txt) | All known and predicted TCM ingredient-target protein interactions |

### File Format Details
- **Format**: Tab-delimited text files
- **Encoding**: UTF-8
- **Schema**: Ingredient ID, Target ID, Interaction Score, Evidence Type

### API Access
BATMAN-TCM 2.0 provides a REST API for programmatic access:
```
http://bionet.ncpsb.org.cn/batman-tcm/api/[parameters]
```

**Supported Outputs**:
- JSON format (for programming)
- HTML (for browser viewing)

### License/Terms
- Academic use permitted
- Citation required for publications

### Processing Recommendations
1. Download tab-delimited files directly from Download page
2. Parse using pandas with `sep='\t'`
3. Filter by confidence score for high-quality predictions
4. Cross-reference with UniProt for target validation

---

## 2. TCMSP

**Traditional Chinese Medicine Systems Pharmacology Database and Analysis Platform**

### Database URL
- **Main Site**: https://tcmsp-e.com/

### Important Notice
**TCMSP does NOT provide official bulk download functionality.** Data must be accessed through web interface queries or spider tools.

### Spider Tool Approach

#### TCMSP-Spider (Recommended)
- **Repository**: https://github.com/shujuecn/TCMSP-Spider
- **License**: MIT

#### Installation
```bash
git clone https://github.com/shujuecn/TCMSP-Spider.git
cd TCMSP-Spider
pip3 install -r requirements.txt
```

#### Usage
1. Configure `herb_list.txt` with drug names (Chinese, Pinyin, or Latin)
2. Run `search_save_herbs.py` for targeted queries
3. Run `get_all_data.py` for comprehensive downloads

#### Output
- Format: Excel files (.xlsx)
- Location: `data/spider_data/` or `data/sample_data/`
- Contents: Drugs, ingredients, targets, diseases

### Rate Limiting Considerations
- Implement delays between requests (recommended: 1-2 seconds)
- Use rotating user agents
- Handle token authentication automatically (tool manages this)
- Monitor for IP blocking; use proxies if necessary

### Alternative: TCMBank
Due to TCMSP's lack of official bulk download, consider using **TCMBank** as a more accessible alternative (see Section 3).

---

## 3. TCMBank

**The Largest TCM Database**

### Database URL
- **Main Site**: https://tcmbank.cn/
- **Download Page**: https://tcmbank.cn/Download

### Database Statistics
| Category | Count |
|----------|-------|
| Herbs | 9,192 |
| Ingredients (unduplicated) | 61,966 |
| Targets | 15,179 |
| Diseases | 32,529 |

### Available Downloads

| Dataset | Format | Direct URL |
|---------|--------|------------|
| All-Herb | XLSX | https://tcmbank.cn/file/TCM_database/herb_all.xlsx |
| All-Ingredient | XLSX | https://tcmbank.cn/file/TCM_database/ingredient_all.xlsx |
| All-Target | XLSX | https://tcmbank.cn/file/TCM_database/gene_all.xlsx |
| All-Disease | XLSX | https://tcmbank.cn/file/TCM_database/disease_all.xlsx |
| All-mol2 | ZIP (mol2 files) | https://tcmbank.cn/file/TCM_database/all_mol2.zip |

### File Schema

#### Herbs (herb_all.xlsx)
| Column | Description |
|--------|-------------|
| TCMBank_ID | Unique identifier |
| TCM_name | Chinese name |
| TCM_name_en | English name |
| Herb_pinyin_name | Pinyin romanization |
| Herb_latin_name | Latin botanical name |
| Function | Traditional functions |
| Indication | Traditional indications |

#### Ingredients (ingredient_all.xlsx)
| Column | Description |
|--------|-------------|
| TCMBank_ID | Unique identifier |
| Ingredient_Name | Compound name |
| Ingredient_Alias | Alternative names |
| Smiles | SMILES structure |
| Molecular_Formula | Chemical formula |
| Molecular_Weight | MW in Daltons |
| mol2_path | Path to 3D structure file |
| OB_score | Oral bioavailability score |

#### Targets (gene_all.xlsx)
| Column | Description |
|--------|-------------|
| TCMBank_ID | Unique identifier |
| Gene_Name | Official gene symbol |
| Gene_alias | Alternative names |
| Description | Gene function description |
| Chromosome | Chromosomal location |

### Cross-References
TCMBank provides links to: CAS, DrugBank, PubChem, MeSH, OMIM, DO, ETCM, HERB

### License/Terms
- Free for academic/non-commercial use
- Citation required

### Processing Recommendations
```python
import pandas as pd

# Load all datasets
herbs = pd.read_excel('herb_all.xlsx')
ingredients = pd.read_excel('ingredient_all.xlsx')
targets = pd.read_excel('gene_all.xlsx')
diseases = pd.read_excel('disease_all.xlsx')

# For 3D structures, extract mol2 files
import zipfile
with zipfile.ZipFile('all_mol2.zip', 'r') as zip_ref:
    zip_ref.extractall('mol2_structures/')
```

---

## 4. HERB 2.0

**High-throughput Experiment- and Reference-guided Database of TCM**

### Database URL
- **Main Site**: http://herb.ac.cn/
- **No registration required**

### Database Contents
| Category | Count |
|----------|-------|
| Herbs | 6,892 |
| Ingredients | 44,595 |
| Formulae (new) | 6,743 |
| Gene Targets | 15,515 |
| Diseases | 30,170 |
| Clinical Trials | 8,558 |
| Meta-analyses | 8,032 |
| High-throughput Experiments | 2,231 |
| Curated References | 6,644 |

### Gene Expression Data Access

#### Re-analyzed Datasets
- 6,164 gene expression profiles from 1,037 GEO experiments
- Mapped to CMap (Connectivity Map) for drug-TCM connections

#### Connectivity Mapping Interface
1. Upload gene expression signature (upregulated/downregulated gene lists)
2. Select reference perturbation dataset
3. Available datasets: 22 herbs, 198 ingredients, 4 formulae, 308 diseases

### Clinical Trial Data Access

#### Data Sources
- ClinicalTrials.gov (8,558 trials)
- PROSPERO database (8,032 meta-analyses)

#### Clinical Conclusions
- 1,941 clinical trials with extracted conclusions
- 593 meta-analyses with clear conclusions

### Knowledge Graph
HERB 2.0 provides knowledge graph with 9 entity types:
- Herbs, Ingredients, Formulae, Targets, Diseases
- Clinical trials, Meta-analyses
- High-throughput experiments, Curated references

### Download Options
**Note**: HERB 2.0 does not appear to offer bulk database downloads. Data access is through:
1. Web interface browsing/searching
2. Gene expression analysis interface
3. API queries (if available)

**Recommendation**: Contact database maintainers for bulk data access requests.

### License/Terms
- Free access, no registration required
- Built with Python-Flask, MySQL backend

---

## 5. IMPPAT 2.0

**Indian Medicinal Plants, Phytochemistry And Therapeutics**

### Database URL
- **Main Site**: https://cb.imsc.res.in/imppat
- **Version**: 2.0 (Released June 17, 2022)
- **GitHub**: https://github.com/asamallab/IMPPAT2

### Database Contents
| Category | Count |
|----------|-------|
| Indian Medicinal Plants | 4,010 |
| Phytochemicals | 17,967 |
| Therapeutic Uses | 1,095 |
| 3D Structures Available | 17,910 |

### Export Options from Web Interface

#### Per-Compound Downloads
From each phytochemical's dedicated page, download:

| Format | Type | Description |
|--------|------|-------------|
| SDF | 2D/3D | Structure Data File |
| MOL | 2D/3D | MDL Molfile |
| MOL2 | 2D/3D | Tripos Mol2 |
| PDB | 3D only | Protein Data Bank format |
| PDBQT | 3D only | AutoDock format |

#### Text Formats
- SMILES
- InChI
- InChIKey

### Bulk Download Considerations

**No official bulk download is available.** Options for obtaining bulk data:

1. **Web Scraping** (with appropriate rate limiting)
2. **Contact Authors**: Request data directly from database maintainers
3. **GitHub Scripts**: Use provided Python scripts for analysis

### GitHub Repository Contents

**URL**: https://github.com/asamallab/IMPPAT2

| Script | Purpose |
|--------|---------|
| ChemicalSimilarityNetwork.py | Calculate Tanimoto coefficients |
| ChemicalStructureImages.py | Generate SVG/PNG images |
| DruglikenessProperties.py | Evaluate drug-likeness |
| MolecularProperties.py | Compute physicochemical properties |
| MolecularScaffolds.py | Analyze molecular scaffolds |
| MurckoScaffold.py | Murcko scaffold decomposition |

### Computed Properties Available
- Physicochemical properties
- Drug-likeness scores (multiple schemes)
- Predicted ADMET properties

### License/Terms
- Code: MIT License
- Database: Academic use (cite publication)

### Processing Recommendations
1. For small datasets: Use web interface downloads
2. For bulk data: Contact authors at IMSc, Chennai
3. For structure analysis: Clone GitHub repo and use provided scripts

---

## 6. KampoDB

**Database of Japanese Kampo Medicines**

### Database URL
- **Main Site**: http://wakanmoview.inm.u-toyama.ac.jp/kampo/
- **Maintainer**: Institute of Natural Medicine, University of Toyama

### Database Contents
| Category | Count |
|----------|-------|
| Kampo Medicines | 42 |
| Crude Drugs | 54 |
| Constituent Compounds | 1,230 |
| Known Target Proteins | 460 |
| Predicted Target Proteins | 1,369 |

### Data Structure (4 Layers)
1. Kampo medicine query
2. Crude drugs forming the Kampo medicine
3. Constituent compounds in crude drugs
4. Target proteins (known + predicted)

### Cross-Reference to KNApSAcK
- Chemical structures sourced from KNApSAcK and PubChem
- Represented by KCF-S (KEGG Chemical Function and Substructures) descriptors

### Data Availability
**Note**: KampoDB does not provide direct bulk download functionality.

**Access Methods**:
1. Web interface browsing
2. Individual compound/formula queries
3. Contact maintainers for bulk data requests

### Functional Annotations
- Biological pathway annotations
- Molecular function annotations
- Docking simulation predictions
- Machine learning-based target predictions

### License/Terms
- **License**: CC BY-SA 4.0
- Original data freely available under Creative Commons

### Processing Recommendations
1. Use web interface for targeted queries
2. For comprehensive access, contact Institute of Natural Medicine
3. Cross-reference compound data with KNApSAcK for additional information

---

## 7. Dr. Duke's Phytochemical Database

**USDA Phytochemical and Ethnobotanical Databases**

### Database URLs
- **Interactive Search**: https://phytochem.nal.usda.gov/
- **Data Catalog**: https://catalog.data.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases-0849e
- **Ag Data Commons**: https://agdatacommons.nal.usda.gov/articles/dataset/Dr_Duke_s_Phytochemical_and_Ethnobotanical_Databases/24660351

### Bulk CSV Download

#### Primary File
| File | URL | Format |
|------|-----|--------|
| Duke-Source-CSV.zip | https://ndownloader.figshare.com/files/43363335 | ZIP (CSV files) |

#### Data Dictionary
| File | Description |
|------|-------------|
| DrDukesDatabaseDataDictionary-prelim.csv | Column descriptions for all tables |

### File Contents
The ZIP archive contains multiple CSV files covering:
- Species information
- Phytochemical compounds
- Biological activities
- Ethnobotanical uses
- References

### Database Scope
- Pharmaceutical applications
- Nutritional research
- Biomedical research
- Alternative therapies
- Herbal products

### License/Terms
- **License**: Creative Commons CCZero (CC0)
- **Public Domain**: Free for any use without restriction

### Processing Recommendations
```python
import pandas as pd
import zipfile
import os

# Extract ZIP file
with zipfile.ZipFile('Duke-Source-CSV.zip', 'r') as zip_ref:
    zip_ref.extractall('duke_data/')

# Load CSV files
for file in os.listdir('duke_data/'):
    if file.endswith('.csv'):
        df = pd.read_csv(f'duke_data/{file}', encoding='utf-8')
        print(f"{file}: {len(df)} rows")
```

---

## 8. LOTUS

**Liberating Occurrence of NaTUral productS**

### Database URLs
- **Main Site**: https://lotus.naturalproducts.net/
- **Download Page**: https://lotus.naturalproducts.net/download
- **GitHub**: https://github.com/lotusnprod/lotus-web
- **Wikidata Entry**: https://www.wikidata.org/wiki/Q104225190

### Database Statistics
- 750,000+ referenced structure-organism pairs
- 290,000+ unique chemical structures
- 40,000+ distinct organisms
- 75,000+ references

### MongoDB Dump Download

#### Download Commands
```bash
# Download MongoDB dump
curl https://lotus.naturalproducts.net/download/mongo -o mongodata/LOTUSlatest.zip

# Extract
unzip mongodata/LOTUSlatest.zip -d mongodata/
cd mongodata/NPOC2021/NPOC2021/

# Restore to MongoDB
mongorestore --port 27019 --db=NPOC2021 --noIndexRestore .
```

### Zenodo Downloads

#### Primary Dataset (v10)
- **DOI**: 10.5281/zenodo.7534071
- **Published**: January 6, 2023

| File | Size | Format | MD5 |
|------|------|--------|-----|
| 230106_frozen.csv.gz | 15.2 MB | Gzip CSV | 4862ab16a7fe916f8ee6930f9bb95f15 |
| 230106_frozen_metadata.csv.gz | 93.0 MB | Gzip CSV | 324579ae4a4785bf1b50b532b0a00eec |
| **Total** | **108.2 MB** | | |

#### Frozen Archive
- **DOI**: 10.5281/zenodo.5794106
- **Timestamp**: 2021-12-20

### Wikidata Integration
LOTUS data is hosted in parallel on Wikidata:
- Community curation enabled
- Novel data additions supported
- SPARQL queries: https://w.wiki/3HLy

### Related Database: COCONUT
COCONUT (COlleCtion of Open NatUral producTs) also available:
- MongoDB dump: coconut-dump-01-2026.sql (31.91 GB)
- CSV lite: coconut_csv_lite-01-2026.zip (191 MB)
- Full CSV: coconut_csv-01-2026.zip (207.9 MB)

### License/Terms
- **License**: Creative Commons Attribution 4.0 International (CC BY 4.0)
- Open source project

### Processing Recommendations
```python
import pandas as pd
import gzip

# Load compressed CSV
with gzip.open('230106_frozen.csv.gz', 'rt') as f:
    lotus_data = pd.read_csv(f)

with gzip.open('230106_frozen_metadata.csv.gz', 'rt') as f:
    lotus_metadata = pd.read_csv(f)
```

---

## 9. KNApSAcK

**Metabolite-Species Relationship Database**

### Database URLs
- **Main Site**: http://www.knapsackfamily.com/KNApSAcK/
- **Alternative**: http://kanaya.naist.jp/KNApSAcK_Family/
- **Manual**: http://www.knapsackfamily.com/KNApSAcK/dl/Manual/KNApSAcKManual.html

### Database Statistics (as of 2024/12/24)
| Category | Count |
|----------|-------|
| Metabolites | 63,715 entries |
| Metabolite-Species Pairs | 159,095 entries |

### KNApSAcK Family Components
| Database | Content |
|----------|---------|
| KNApSAcK Core | 101,500 species-metabolite relationships |
| KNApSAcK WorldMap | 41,548 GZ-plant pairs (222 zones, 15,240 plants) |
| KAMPO DB | 336 formulae, 278 medicinal plants |
| JAMU DB | 5,310 formulae, 550 medicinal plants |
| Biological Activity DB | 2,418 activities, 33,706 relationships |

### Download Options

#### Downloadable Versions
| Format | File | Version |
|--------|------|---------|
| ZIP | KNApSAcK-v1.200.03.zip | 1.200.03 (2008/12/19) |
| LZH | KNApSAcK-v1.200.03.lzh | 1.200.03 (2008/12/19) |
| TAR.GZ | KNApSAcK-v1.200.03.tar.gz | 1.200.03 (2008/12/19) |
| JAR | KNApSAcK.jar | PC application |

**Note**: The downloadable versions are older (2008). For current data, use the web interface or contact maintainers.

### Web Access
- Keyword search interface available
- Mass spectra search for metabolomics
- Search by: accurate mass, molecular formula, metabolite name, mass spectra

### License/Terms
- Academic use permitted
- Maintained by NAIST (Nara Institute of Science and Technology)

### Processing Recommendations
1. For current data: Use web interface queries
2. For bulk access: Contact maintainers at NAIST
3. Cross-reference with KampoDB for Kampo-specific data

---

## 10. Processing Recommendations

### General Workflow

```python
#!/usr/bin/env python3
"""
Traditional Medicine Database Processing Pipeline
"""

import pandas as pd
import requests
import zipfile
import gzip
from pathlib import Path

# Create data directory
data_dir = Path('traditional_medicine_data')
data_dir.mkdir(exist_ok=True)

# 1. TCMBank Downloads
tcmbank_files = {
    'herbs': 'https://tcmbank.cn/file/TCM_database/herb_all.xlsx',
    'ingredients': 'https://tcmbank.cn/file/TCM_database/ingredient_all.xlsx',
    'targets': 'https://tcmbank.cn/file/TCM_database/gene_all.xlsx',
    'diseases': 'https://tcmbank.cn/file/TCM_database/disease_all.xlsx',
    'mol2': 'https://tcmbank.cn/file/TCM_database/all_mol2.zip'
}

for name, url in tcmbank_files.items():
    print(f"Downloading TCMBank {name}...")
    response = requests.get(url)
    suffix = '.xlsx' if not name == 'mol2' else '.zip'
    with open(data_dir / f'tcmbank_{name}{suffix}', 'wb') as f:
        f.write(response.content)

# 2. Dr. Duke's Database
print("Downloading Dr. Duke's Database...")
duke_url = 'https://ndownloader.figshare.com/files/43363335'
response = requests.get(duke_url)
with open(data_dir / 'Duke-Source-CSV.zip', 'wb') as f:
    f.write(response.content)

# 3. LOTUS from Zenodo
print("Downloading LOTUS data...")
# Note: Adjust Zenodo URLs based on latest version
lotus_base = 'https://zenodo.org/records/7534071/files/'
for file in ['230106_frozen.csv.gz', '230106_frozen_metadata.csv.gz']:
    response = requests.get(f'{lotus_base}{file}')
    with open(data_dir / file, 'wb') as f:
        f.write(response.content)

print("Downloads complete!")
```

### Data Integration Strategy

1. **Standardize Identifiers**
   - Map to PubChem CIDs where possible
   - Use InChIKey for structure matching
   - Map gene targets to UniProt IDs

2. **Structure Validation**
   - Validate SMILES with RDKit
   - Standardize stereochemistry
   - Remove salts and solvents

3. **Cross-Database Linking**
   - TCMBank <-> TCMSP via compound names/structures
   - KampoDB <-> KNApSAcK via compound IDs
   - LOTUS <-> Wikidata via QIDs

### Storage Recommendations

| Database | Recommended Storage |
|----------|---------------------|
| Small (<1GB) | SQLite, CSV |
| Medium (1-10GB) | PostgreSQL, MySQL |
| Large (>10GB) | MongoDB, PostgreSQL with partitioning |
| Structures | Separate mol2/SDF directories |

### Citation Requirements

When using these databases, cite the original publications:

- **BATMAN-TCM 2.0**: Liu et al., Nucleic Acids Research 2024
- **TCMSP**: Ru et al., Journal of Cheminformatics 2014
- **TCMBank**: Xu et al., Chemical Science 2023
- **HERB 2.0**: Gao et al., Nucleic Acids Research 2024
- **IMPPAT 2.0**: Vivek-Ananth et al., ACS Omega 2023
- **KampoDB**: Sawada et al., Scientific Reports 2018
- **Dr. Duke's**: James A. Duke, USDA
- **LOTUS**: Rutz et al., eLife 2022
- **KNApSAcK**: Shinbo et al., Plant Cell Physiology 2006

---

## Summary Table

| Database | Bulk Download | Format | Size | License |
|----------|--------------|--------|------|---------|
| BATMAN-TCM 2.0 | Yes (Download page) | Tab-delimited | ~100MB+ | Academic |
| TCMSP | No (use spider) | - | - | - |
| TCMBank | Yes (direct URLs) | XLSX, mol2 | Variable | Free/Academic |
| HERB 2.0 | Limited | Web interface | - | Free |
| IMPPAT 2.0 | Per-compound | SDF, MOL, MOL2 | ~17K compounds | Academic |
| KampoDB | No | Web interface | ~1K compounds | CC BY-SA 4.0 |
| Dr. Duke's | Yes (CSV ZIP) | CSV | ~50MB | CC0 |
| LOTUS | Yes (MongoDB, CSV) | CSV.gz, MongoDB | 108 MB | CC BY 4.0 |
| KNApSAcK | Yes (older version) | ZIP, JAR | Variable | Academic |

---

*Document created: 2026-01-18*
*Last updated: 2026-01-18*
