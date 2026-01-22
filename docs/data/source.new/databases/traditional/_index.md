---
title: "Traditional Medicine Databases"
parent: ../_index.md
world: 2
last_updated: 2026-01-22
status: draft
---

# Traditional Medicine Databases (World 2)

Databases for traditional medicine systems including TCM, Ayurveda, Kampo, Western Herbal, and global ethnobotanical knowledge.

## Database Catalog

| Database | System | Tier | Coverage | Access Method | Size |
|----------|--------|------|----------|---------------|------|
| **BATMAN-TCM** | TCM | 1 | Formulas, herbs, compounds | Web, API (unofficial) | 8,000+ compounds |
| **HERB** | TCM | 1 | Herb-gene-disease | Web scraping | 7,000+ herbs |
| **TCMSID** | TCM | 1 | TCM syndrome-gene | Web, Download | 3,500+ syndromes |
| **IMPPAT** | Ayurveda | 1 | Phytochemicals | Download | 9,500+ compounds |
| **NPASS** | Multi-system | 2 | Natural product activity | Web, Download | 35,000+ compounds |
| **KampoDB** | Kampo | 2 | Japanese herbal formulas | Literature-based | 200+ formulas |
| **TM-MC** | Multi-system | 2 | Traditional medicine compounds | Web | 30,000+ compounds |
| **ETCM** | TCM | 2 | Network pharmacology | Web | 7,000+ herbs |
| **SymMap** | TCM | 2 | Symptom mapping | Web | 1,717 TCM symptoms |
| **YaTCM** | TCM | 3 | Yet another TCM database | Web | Comprehensive |
| **TCMIP** | TCM | 3 | TCM integrated platform | Web | Multi-omics |
| **Indigenous Databases** | Global | 3 | Regional ethnobotany | Varies | Scattered |

## Traditional Medicine Systems

### Traditional Chinese Medicine (TCM)

#### BATMAN-TCM (Tier 1)
- **Purpose**: Bioinformatic Analysis Tool for Molecular mechANism of TCM
- **Content**: 8,000+ compounds, 50+ herbs, target predictions
- **Use**: Formula analysis, target prediction, network pharmacology
- **Access**: http://bionet.ncpsb.org/batman-tcm/
- **Integration**: Links to KEGG, DrugBank, PubChem

#### HERB (Tier 1)
- **Purpose**: High-throughput Experiment and Reference-guided database of traditional Chinese medicine
- **Content**: 7,000+ herbs, 50,000+ ingredients, 500+ diseases
- **Use**: Herb-gene-disease relationships
- **Access**: http://herb.ac.cn/ (web scraping required)
- **Integration**: Gene targets, pathways, diseases

#### TCMSID (Tier 1)
- **Purpose**: TCM Syndrome-Ingredient-Disease
- **Content**: 3,500+ syndromes, 1,000+ herbs, disease mappings
- **Use**: Syndrome differentiation, personalized TCM
- **Access**: http://www.tcmsid.com/
- **Integration**: ICD codes, symptom ontologies

### Ayurvedic Medicine

#### IMPPAT (Tier 1)
- **Purpose**: Indian Medicinal Plants, Phytochemistry And Therapeutics
- **Content**: 9,500+ phytochemicals, 1,700+ plants, 2,700+ therapeutic uses
- **Use**: Ayurvedic compound-disease relationships
- **Access**: https://cb.imsc.res.in/imppat/ (bulk download available)
- **Integration**: PubChem, ChEBI, therapeutic categories

### Kampo (Japanese Herbal Medicine)

#### KampoDB (Tier 2)
- **Purpose**: Kampo formula database
- **Content**: 200+ standardized formulas, ingredient lists
- **Use**: Kampo prescription analysis
- **Access**: Literature-based extraction required
- **Integration**: Similar to TCM herbs, KEGG pathways

### Western Herbal Medicine

#### NAPRALERT (Tier 3)
- **Purpose**: Natural Products Alert
- **Content**: 200,000+ scientific articles on natural products
- **Use**: Western herbal research
- **Access**: Subscription required
- **Integration**: Chemical structures, biological activities

## Cross-System Databases

### NPASS (Tier 2)
- **Purpose**: Natural Product Activity and Species Source
- **Content**: 35,000+ natural products, 1,500+ species, activity data
- **Use**: Cross-system natural product analysis
- **Access**: http://bidd.group/NPASS/index.php
- **Integration**: PubChem, ChEBI, UniProt

### TM-MC (Tier 2)
- **Purpose**: Traditional Medicine Molecular Compounds
- **Content**: 30,000+ compounds from multiple traditions
- **Use**: Multi-system compound analysis
- **Access**: Web interface
- **Integration**: Chemical structures, biological targets

## Data Integration Workflow

```
Traditional Formula/Herb
         ↓
BATMAN-TCM / HERB / IMPPAT
         ↓
    Compounds
         ↓
┌────────┴────────┐
│                 │
NPASS          COCONUT
│                 │
Activity       Structure
│                 │
└────────┬────────┘
         ↓
   Target Genes
         ↓
    Pathways
```

## Access Strategies

### Tier 1 Databases

#### BATMAN-TCM
```bash
# Web interface
# http://bionet.ncpsb.org/batman-tcm/

# Unofficial API (if available)
curl "http://bionet.ncpsb.org/batman-tcm/api/formula?name=SiJunZiTang"
```

#### HERB
```python
# Web scraping approach
import requests
from bs4 import BeautifulSoup

url = "http://herb.ac.cn/herb/detail?id=HERB000001"
response = requests.get(url)
soup = BeautifulSoup(response.content, 'html.parser')
# Parse herb data
```

#### IMPPAT
```bash
# Bulk download available
wget https://cb.imsc.res.in/imppat/downloads/all_compounds.csv
wget https://cb.imsc.res.in/imppat/downloads/plant_therapeutic_use.csv
```

### Data Extraction Challenges

| Database | Challenge | Solution |
|----------|-----------|----------|
| BATMAN-TCM | No official API | Web scraping with rate limiting |
| HERB | Dynamic content | Selenium/Playwright |
| TCMSID | Limited bulk access | Incremental scraping |
| IMPPAT | Large files | Streaming download |
| NPASS | Mixed formats | Custom parsers |

## Traditional System Mapping

### TCM Syndrome → ICD Mapping
```
TCM Syndrome (TCMSID)
         ↓
Symptom Ontology
         ↓
    ICD-10/11
         ↓
Western Disease
```

### Herb → Compound → Gene Pipeline
```
TCM Herb (HERB)
         ↓
Compounds (BATMAN-TCM)
         ↓
Target Genes (HERB, NPASS)
         ↓
Pathways (KEGG, Reactome)
         ↓
Western Disease (DisGeNET)
```

## Update Strategy

### Tier 1 Updates
- **BATMAN-TCM**: Check quarterly for new formulas
- **HERB**: Semi-annual full refresh
- **TCMSID**: Annual updates
- **IMPPAT**: Check for new releases annually

### Tier 2 Updates
- **NPASS**: Annual refresh
- **TM-MC**: Check biannually
- **ETCM**: As needed

## Storage Estimates

| Database | Storage Required | Format |
|----------|------------------|--------|
| BATMAN-TCM | 500 MB | HTML, JSON |
| HERB | 1-2 GB | HTML, CSV |
| TCMSID | 200-500 MB | HTML, TSV |
| IMPPAT | 100-200 MB | CSV, JSON |
| NPASS | 500 MB - 1 GB | CSV, SDF |

## Data Quality Considerations

### Standardization Issues
- Herb name variations (Chinese, pinyin, Latin)
- Compound structure inconsistencies
- Syndrome description variations
- Dosage unit differences

### Integration Challenges
- Mapping TCM syndromes to ICD codes
- Linking traditional herbs to modern chemical databases
- Reconciling different therapeutic classification systems
- Handling regional variations in formulas

## Navigation

- **Parent**: [Database Sources](../_index.md)
- **Related**: [Compounds](../compounds/_index.md), [Nutrition](../nutrition/_index.md)
