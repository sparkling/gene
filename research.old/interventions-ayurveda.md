# Ayurveda, Siddha, and Unani Medicine Data Sources

## Executive Summary

This document catalogs available data sources for traditional Indian medicine systems (Ayurveda, Siddha, Unani, and Sowa Rigpa) for integration into a genetics/health knowledge base platform. The Indian traditional medicine landscape offers several well-curated databases with varying levels of accessibility, from fully open-access resources to restricted government repositories.

**Key Findings:**
- **Best Overall Resource**: IMPPAT 2.0 - Most comprehensive, open access (CC BY 4.0), includes predicted protein targets
- **Largest Plant Coverage**: OSADHI - 6,959 unique plants, 27,440 phytochemicals
- **Best for Formulations**: GRAYU - 1,039 formulations mapped to plants, compounds, and diseases
- **Most Comprehensive Traditional Knowledge**: TKDL - 454,000+ formulations but restricted access
- **Best for Anti-cancer Research**: NPACT - 1,574 compounds with validated cancer targets

---

## 1. IMPPAT 2.0 (Indian Medicinal Plants, Phytochemistry And Therapeutics)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://cb.imsc.res.in/imppat/ |
| **Maintainer** | The Institute of Mathematical Sciences (IMSc), Chennai, India |
| **Latest Version** | 2.0 (Released June 17, 2022) |
| **Publication** | ACS Omega, 8(9):8827-8845, 2023 |

### Content Coverage
| Category | Count |
|----------|-------|
| Medicinal Plants | 4,010 |
| Phytochemicals | 17,967 |
| Therapeutic Uses | 1,095 |
| Plant-Phytochemical Associations | 5-fold increase over v1.0 |
| Plant-Therapeutic Associations | 5-fold increase over v1.0 |
| Predicted Phytochemical-Target Interactions | 27,365 |
| Human Target Proteins | 5,042 |
| Drug-like Phytochemicals | 1,335 (filtered subset) |

### Data Associations
Two primary relationship types:
1. **Plant -> Plant Part -> Phytochemical**
2. **Plant -> Plant Part -> Therapeutic Use**

### Access Method
- **Web Interface**: Browse, basic search, advanced search with filters
- **Export**: Tab-separated lists via web interface export option
- **Structure Downloads**: SDF, MOL, MOL2, PDB, PDBQT formats for individual phytochemicals
- **API**: No public API documented

### Data Format & Schema
For each phytochemical, six data tabs are available:
1. **Summary**: IMPPAT ID, IUPAC name, synonyms, chemical class, InChI, InChI key, SMILES, CAS number
2. **Physicochemical**: Molecular weight, logP, H-bond donors/acceptors, rotatable bonds
3. **Drug-likeness**: Lipinski, Ghose, Veber, Egan scores
4. **ADMET**: Predicted absorption, distribution, metabolism, excretion, toxicity properties
5. **Descriptors**: Chemical descriptors
6. **Predicted Human Targets**: Protein interaction predictions (score >= 700)

### Gene/Protein Targets
**YES** - Includes 27,365 predicted interactions between 1,294 phytochemicals and 5,042 human target proteins (high confidence interactions only, combined score >= 700).

### Licensing
**Creative Commons Attribution 4.0 International (CC BY 4.0)** - Permits use, sharing, adaptation, distribution and reproduction in any medium or format with appropriate credit.

### Data Sources
- 100+ books on traditional Indian medicine
- 7,000+ published research articles
- Existing phytochemical databases

### GitHub Repository
https://github.com/asamallab/IMPPAT2 - Contains codes for cheminformatics analysis

---

## 2. OSADHI (Online Structural and Analytics Database for Herbs of India)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://neist.res.in/osadhi/ |
| **Maintainer** | CSIR-North East Institute of Science and Technology (NEIST), Jorhat, India |
| **Publication** | Computational Biology and Chemistry (2022) |

### Content Coverage
| Category | Count |
|----------|-------|
| Total Medicinal Plants Identified | 21,238 |
| Unique Medicinal Plants | 6,959 |
| Plant Families | 348 |
| Phytochemicals | 27,440 (for unique plants) |
| Detailed Phytochemical Analysis | 22,314 |
| Therapeutic Uses | 2,477 |
| States/Union Territories Covered | 28 states + 8 UTs |

### Key Features (4Ds Concept)
1. **Documentation** - Traditional knowledge
2. **Digitization** - Vernacular names, plant parts, therapeutic use, taxonomy
3. **Deposition** - Geographical distribution
4. **Data Science** - Chemoinformatics analysis

### Data Structure
Four major sections:
1. **Traditional Knowledge**: Vernacular names, plant parts used, therapeutic use, taxonomy
2. **Geographical Classification**: Distribution across Indian states and Union Territories
3. **Phytochemicals**: SMILES, InChIKey, IUPAC, 2D & 3D structures
4. **Chemoinformatics**: Physicochemical properties, ADMET, NPClassifier classification, antiviral potency predictions (Random Forest & XGBoost)

### Access Method
- **Web Interface**: Open access
- **Structure Downloads**: 2D and 3D structures available
- **API**: Not documented

### Gene/Protein Targets
**NO** - Focus is on plant-compound relationships and properties, not target predictions.

### Licensing
Open access database (specific license terms not explicitly stated on website)

---

## 3. GRAYU (Graph-based Database integrating Ayurvedic formulations)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://caps.ncbs.res.in/GRAYU/ |
| **Maintainer** | National Centre for Biological Sciences (NCBS), Bangalore, India |
| **Publication** | Frontiers in Pharmacology (2025) |

### Content Coverage
| Category | Count |
|----------|-------|
| Ayurvedic Formulations | 1,039 |
| Indigenous Plants | 12,743-12,949 |
| Phytochemicals | 129,542 |
| Diseases | 13,480 |
| Plant-Phytochemical Associations | 1,382,362 |
| Plant-Disease Associations | 116,824 |
| Plant-Formulation Associations | 2,405 |
| Formulation-Disease Associations | 4,087 |

### Key Features
- **Graph-based Architecture**: Specifically designed for Ayurveda
- **Disease Ontology Mapping**: Links Sanskrit nosology with MeSH and DOID disease ontologies
- **Network Pharmacology Ready**: Designed for computational drug discovery

### Data Sources Integration
Aggregated from: CMAUP, FooDB, HMDB, IMPPAT, NPASS, PubChem

### Access Method
- **Web Interface**: Category-specific searches, Interactive Knowledge Graph
- **Advanced Search**: Available
- **Download**: Not explicitly documented
- **API**: Not documented

### Gene/Protein Targets
**INDIRECT** - Through disease ontology mappings and phytochemical associations; no direct target data in database itself.

### Licensing
**Research & Education Only** - No clinical use permitted. Citation required when publishing results.

---

## 4. DHARA (Digital Helpline for Ayurveda Research Articles)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://dharaonline.org/ |
| **Maintainer** | AVP Research Foundation (formerly AVTAR), The Ayurvedic Trust, Coimbatore, Tamil Nadu |
| **Funding History** | CCRAS/Ministry of AYUSH (2010-2011), AVP Research Foundation (2012-present) |

### Content Coverage
| Category | Count |
|----------|-------|
| Total Indexed Articles | 52,000+ (earlier data) / 10,583 (current website) |
| Journals Covered | 4,100+ worldwide |
| Journals (Ayurveda-specific) | 17 |
| Journals (CAM-related) | 18 |
| Journals (Mainstream medical) | 846 |
| Full Text Available | 4,064 (3,185 free, 780 paid) |
| Abstract-Only Articles | 5,368 |
| Title-Only Articles | 1,151 |
| Unique Authors | 24,129 |

### Access Method
- **Web Interface**: Keyword search, advanced search with Boolean operators, field tags
- **Free Access**: No registration required for basic search
- **Download**: Links to full text based on journal policies
- **API**: Not available

### Data Format
Bibliographic index format - titles, abstracts, links to full text

### Gene/Protein Targets
**NO** - This is a literature index, not a compound/target database.

### Licensing
Free public service for academic research

### Importance
Less than 10% of articles indexed in DHARA can be retrieved using keyword "Ayurveda" in international databases like PubMed, making DHARA essential for comprehensive Ayurveda literature searches.

---

## 5. TKDL (Traditional Knowledge Digital Library)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://www.tkdl.res.in/ |
| **Maintainer** | Council of Scientific and Industrial Research (CSIR) & Ministry of AYUSH, Government of India |
| **Established** | 2001 |

### Content Coverage (as of March 2022)
| System | Formulations Transcribed |
|--------|-------------------------|
| **Ayurveda** | 119,269 |
| **Unani** | 236,399 |
| **Siddha** | 54,689 |
| **Yoga** | 4,151 |
| **Sowa Rigpa** | 4,377 |
| **Total** | 454,885+ |

### Classification System (TKRC)
Traditional Knowledge Resource Classification (TKRC):
- ~5,000 subgroups created for Ayurveda, Unani, Siddha, and Yoga
- Correlates traditional terminology with modern terminology
- Recognized by WIPO, incorporated into IPC under A61K 36/00

### Languages Available
Content transcribed into 5 international languages: English, Japanese, French, German, Spanish

### Access Method
**RESTRICTED ACCESS**
- **Current Access**: Patent Offices only under TKDL Access Agreement
  - 16 Patent Offices have access (EPO, USPTO, JPO, UKIPO, CIPO, DPMA, IP Australia, Indian Patent Office, Chile, Malaysia, etc.)
- **Future Access**: Paid subscription model being implemented for:
  - Businesses (herbal healthcare, AYUSH, pharmaceuticals, nutraceuticals, personal care, FMCG)
  - Research institutions (public and private)
  - Educational institutions
  - ISM practitioners, knowledge holders, patentees

### Search Features
- Full text search and retrieval
- Boolean expression search (AND, OR, AND NOT, NEAR)
- Proximity search
- Field search
- Phrase search
- IPC and TKRC code searches

### Gene/Protein Targets
**NO** - Focus is on traditional formulations for patent prior art, not molecular targets.

### Licensing
**Highly Restricted** - Access agreements with built-in safeguards on non-disclosure

### Contact
Head, CSIR-Traditional Knowledge Digital Library Unit: headtkdl@csir.res.in

---

## 6. AYUSH Research Portal

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://ayushportal.nic.in/ or https://arp.ayush.gov.in/ |
| **Maintainer** | Ministry of AYUSH, Government of India / CCRAS-NIIMH, Hyderabad |

### Content Coverage
| Category | Count |
|----------|-------|
| Article Views | 2,111,343+ |
| Downloads | 292,395+ |
| Systems Covered | Ayurveda, Yoga, Naturopathy, Unani, Siddha, Sowa-Rigpa, Homeopathy |

### Research Categories
- Clinical Research
- Pre-clinical Research
- Drug Research
- Fundamental Research

### Access Method
- **Web Interface**: Search AYUSH terminology, research articles, journals
- **Registration**: Sign-up available for different user categories
- **Download**: Individual article downloads
- **API**: Not documented

### Gene/Protein Targets
**NO** - Focus on research articles and publications, not compound-target data.

### Licensing
Academic use only. Ministry of AYUSH not responsible for findings/claims in uploaded content.

---

## 7. NPACT (Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | http://crdd.osdd.net/raghava/npact/ |
| **Maintainer** | CSIR-Institute of Microbial Technology (IMTECH), Chandigarh, India |
| **Publication** | Nucleic Acids Research, 41, D1124-D1129, 2013 |

### Content Coverage
| Category | Count |
|----------|-------|
| Bioactive Compounds | 1,574 |
| Cancer Cell Lines | 353 |
| Cancer Types | 27 |
| Compound Classes | 19 |
| Compound-Cell Line Interactions | ~5,214 (in vitro) |
| Compound-Target Interactions | ~1,980 (experimentally validated) |
| Commercial Suppliers | 50 |

### Data Fields Per Entry
- Compound name, NPACT ID, IUPAC name, synonyms
- Chemical class, InChI, InChI key, SMILES, SMART, CAS number
- Physical, elemental, and topological properties
- Cancer type, cell lines, inhibitory values (IC50, ED50, EC50, GI50)
- Molecular targets
- Commercial suppliers
- Drug-likeness scores

### Access Method
- **Web Interface**: Browse by alphabetical order, cancer type, cell line, class name, supplier
- **Structure Download**: MOL format for 3D structures
- **Bulk Download**: Available through ZINC database (ZINC Catalog NPACT)
- **API**: Not documented

### Gene/Protein Targets
**YES** - Contains ~1,980 experimentally validated compound-target interactions for anti-cancer compounds.

### Licensing
Open access (specific license not stated)

---

## 8. CMAUP (Collective Molecular Activities of Useful Plants)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://bidd.group/CMAUP/ |
| **Maintainer** | Bioinformatics and Drug Design Group, National University of Singapore |
| **Publication** | Nucleic Acids Research, 47(D1), D1118-D1127, 2019 |

### Content Coverage
| Category | Count |
|----------|-------|
| Plants | 5,645 (2,567 medicinal, 170 food, 1,567 edible, 3 agricultural, 119 garden) |
| Countries/Regions Represented | 153 |
| Phytochemicals | ~48,000 |
| Active Plant Ingredients | 47,645 |
| Targets | 646 |
| KEGG Pathways | 234 |
| Gene Ontologies | 2,473 |
| Diseases | 656 |

### Key Features
- Collective landscapes of multiple targets (ChEMBL target classes)
- Activity levels in 2D target-ingredient heatmaps
- Regulated gene ontologies, biological pathways, processes, diseases

### Access Method
- **Web Interface**: Search by keywords, plant usage classes, species families, targets, KEGG pathways, gene ontologies, diseases (ICD code), geographical locations
- **Interactive World Map**: Browse plants by location
- **Download**: All data freely downloadable
- **API**: Not documented

### Gene/Protein Targets
**YES** - Contains 646 targets with activity data from ChEMBL target classes.

### Licensing
Freely accessible (specific license not stated)

---

## 9. NPASS (Natural Product Activity and Species Source)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | http://bidd.group/NPASS |
| **Maintainer** | Bioinformatics and Drug Design Group, National University of Singapore |
| **Latest Version** | 3rd version (2023 update) |
| **Publication** | Nucleic Acids Research (2023) |

### Content Coverage (2023 Update)
| Category | Count |
|----------|-------|
| Unique Natural Products | 35,032+ (original); significantly expanded |
| Species Sources | 25,041+ |
| Targets | ~7,700 |
| NP-Species Pairs | 288,002+ |
| NP-Target Pairs | 222,092+ |
| NP Composition/Concentration Records | ~95,000 |
| NPs with Estimated Activity Profiles | ~66,600 |

### Species Sources Breakdown
- Viridiplantae (plants): 67.8%
- Metazoan (animals): 9.4%
- Fungi: 7.9%
- Bacteria: 6.7%

### Data Fields
- NP names, activity values (IC50, Ki, EC50, GI50, MIC)
- Target types (proteins, protein families, complexes, nucleic acids, cell lines, tissue targets)
- Species source taxonomy
- Drug-likeness properties
- ADMET properties
- Chemical Checker estimated activity profiles

### Access Method
- **Web Interface**: Browse by NP names, targets, target types, species, molecular weights
- **Download**: All data freely downloadable from entry pages
- **Download Summary Table**: Select specific data sections
- **API**: Not documented

### Gene/Protein Targets
**YES** - Contains ~222,092 NP-target pairs with experimentally determined activity values.

### Licensing
Freely accessible (specific license not stated)

---

## 10. Dr. Duke's Phytochemical and Ethnobotanical Databases

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://phytochem.nal.usda.gov/ |
| **Maintainer** | USDA National Agricultural Library |
| **Original Developer** | James A. Duke, USDA Economic Botany Laboratory |

### Content Coverage
Comprehensive repository of ethnobotanical data including:
- Species information
- Phytochemicals
- Biological activities
- Ethnobotanical uses (including Ayurvedic applications)

### Access Method
- **Web Interface**: Plant, chemical, bioactivity, ethnobotany searches
- **Download**: Search results downloadable in spreadsheet form
- **Bulk Download**: Raw CSV database tables available
- **Data Dictionary**: DrDukesDatabaseDataDictionary-prelim.csv available
- **API**: Not available

### Download Sources
- **Ag Data Commons**: https://data.nal.usda.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases
- **Figshare**: https://agdatacommons.nal.usda.gov/articles/dataset/Dr_Duke_s_Phytochemical_and_Ethnobotanical_Databases/24660351
- **Data.gov**: https://catalog.data.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases-0849e

### Gene/Protein Targets
**LIMITED** - Contains biological activity data but not systematic target predictions.

### Licensing
**Creative Commons CC0 Public Domain** - Completely open for any use.

---

## 11. AromaDb (CIMAP Aroma Molecules Database)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://aromadb.cimapbioinfo.in/ or https://bioinfo.cimap.res.in/aromadb/ |
| **Maintainer** | CSIR-Central Institute of Medicinal and Aromatic Plants (CIMAP), Lucknow, India |
| **Publication** | Frontiers in Plant Science, 2018 |

### Content Coverage
| Category | Count |
|----------|-------|
| Aroma Chemical Structures | 1,321 |
| Fragrance Types | 357 |
| Commercially Used Plants | 166 |
| High-Yielding Varieties/Chemotypes | 148 |

### Data Fields
- IUPAC name, synonyms, CAS registry number
- Chemical classification, functional groups
- Molecular weight
- 2D and 3D structures (volatile <300 MW, medium <500 MW)
- Pharmacokinetics (ADMET) data
- Bioactivities of essential oil/aroma compounds
- Biological pathways information
- Cross-references

### Access Method
- **Web Interface**: Simple and advanced search
- **Structure Download**: 2D & 3D structures freely downloadable
- **Similarity Search**: Based on structural similarity
- **Property Search**: Based on physicochemical and ADMET properties
- **API**: Not documented

### Gene/Protein Targets
**LIMITED** - Contains bioactivity data but not systematic target mapping.

### Licensing
Open access (specific license not stated)

---

## 12. sCentInDB (Essential Oil Chemical Profiles Database)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://cb.imsc.res.in/scentindb/ |
| **Maintainer** | The Institute of Mathematical Sciences (IMSc), Chennai, India |
| **Publication** | Molecular Diversity, 2025 |

### Content Coverage
| Category | Count |
|----------|-------|
| Indian Medicinal Plants | 554 |
| Essential Oil Profiles | 2,170 |
| Chemicals | 3,420 |
| Plant-Part-Therapeutic Use Associations | 471 |
| Plant-Part-Odor Associations | 120 |
| Plant-Part-Color Associations | 218 |
| Source Research Articles | 778 |

### Key Features
- Data compiled at plant-part level
- Sensory attributes (odor, color)
- Pharmacokinetic data included
- FAIR-compliant design

### Access Method
- **Web Interface**: Open access, no account required
- **Download**: CSV format (machine-readable)
- **HTML Display**: Human-readable format
- **API**: Not documented

### Gene/Protein Targets
**NO** - Focus on essential oil chemical profiles, not molecular targets.

### Licensing
Academic research use permitted freely

---

## 13. Indian Medicinal Plants Database (FRLHT/NMPB)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://www.medicinalplants.in/ |
| **Maintainer** | Foundation for Revitalisation of Local Health Traditions (FRLHT), Bangalore |
| **Funding** | National Medicinal Plants Board (NMPB) |

### Content Coverage
| Category | Count |
|----------|-------|
| Botanical Names | 7,263 |
| Vernacular Names | 150,000+ (in 10 Indian languages) |
| Authentic Plant Images | 5,000+ |
| Indian Systems of Medicine Covered | 6 (Ayurveda, Unani, Siddha, Homeopathy, etc.) |

### Data Fields
- Botanical names and synonyms
- Vernacular names in regional languages
- Medicinal uses by system of medicine
- Plant images linked to botanical entities

### Access Method
- **Web Interface**: Registration required (user ID and password)
- **Search**: By botanical name or vernacular name within each system of medicine
- **API**: Not available

### Gene/Protein Targets
**NO** - Focus on botanical and traditional use information.

### Licensing
**Research and Development Only** - Commercial use strictly prohibited.

---

## 14. IMPDB (Indian Medicinal Phytochemical Database)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://impdb.org.in/ |
| **Publication** | Journal of Computational Biology and Chemistry, 2022 |
| **Technical Stack** | PHP (front-end), MySQL (back-end) |

### Content Coverage
- Phytochemicals from "Database on Medicinal Plants Used in Ayurveda" (CCRAS publication)
- Physicochemical properties calculated using Discovery Studio v4.5
- Drug-likeness properties

### Key Finding
99.96-99.99% of phytochemicals in IMPDB are novel compared to 2,068 FDA-approved drugs.

### Access Method
- **Web Interface**: Standard database interface
- **Download**: Not explicitly documented
- **API**: Not documented

### Gene/Protein Targets
**NO** - Focus on druggability analysis, not target predictions.

### Licensing
Not explicitly stated

---

## 15. ABIM (Annotated Bibliography of Indian Medicine)

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://indianmedicine.eldoc.ub.rug.nl/ |
| **Maintainer** | University of Groningen, Netherlands |
| **Developer** | Jan Meulenbeld |
| **Platform** | EPrints 3, University of Southampton |

### Content Coverage
| Category | Count |
|----------|-------|
| Total Citations | 50,000+ |

### Search Features
- Simple search (no Boolean operators)
- Displays diacritical marks in Sanskrit terms

### Access Method
- **Web Interface**: Search and browse
- **API**: Not available

### Gene/Protein Targets
**NO** - Bibliographic database only.

### Licensing
Academic/research use

---

## 16. e-Charak Portal

### Database Information
| Field | Details |
|-------|---------|
| **URL** | https://echarak.ayush.gov.in/ |
| **Maintainer** | Ministry of AYUSH, Government of India / CCRAS |

### Content Coverage
- 22,000+ references (reprints/abstracts) for 16 selected medicinal plants
- Classical literature from Ayurveda, Unani, Siddha, Homeopathy
- Modern literature from various books and journals

### Categories
1. Botany
2. Chemistry
3. Pharmacology
4. Pharmacy
5. Miscellaneous
6. Reprints in other languages

### Access Method
- **Web Interface**: Registration required
- **Search**: By plant name
- **API**: Not available

### Gene/Protein Targets
**NO** - Literature reference database.

### Licensing
Research and development only. Commercial use prohibited.

---

## Siddha Medicine Databases

### Expert System for Siddha (eSS)

| Field | Details |
|-------|---------|
| **Status** | Pilot project |
| **Content** | 110 medicinal preparations from 2 Siddhars (Agathiyar and Theran) |
| **Purpose** | Extract, annotate, and visualize patterns in Siddha medical prescriptions |

### TKDL Siddha
| Category | Formulations |
|----------|--------------|
| Siddha Formulations | 54,689 |

**Note**: TKDL remains the largest repository of Siddha formulations but has restricted access.

### Challenges
- Entire Siddha literature is in Tamil language (poems/padal format)
- Language barrier limits international accessibility
- Need for structured extraction systems

---

## Unani Medicine Databases

### TKDL Unani
| Category | Formulations |
|----------|--------------|
| Unani Formulations | 236,399 (largest in TKDL) |

### National Formulary of Unani Medicine (NFUM)
- 6 volumes published containing 1,228 standardized formulations

### Unani Pharmacopoeia of India (UPI)
- 298 single drugs (6 volumes, Part I)
- 100 compound drugs (2 volumes, Part II)

### CCRUM (Central Council for Research in Unani Medicine)
| Field | Details |
|-------|---------|
| **URL** | https://ccrum.res.in/ |
| **Research Units** | 31 institutes/units across India |

---

## Sowa Rigpa (Tibetan Medicine) Databases

### Dataset of Materia Medica in Sowa Rigpa

| Field | Details |
|-------|---------|
| **Publication** | ScienceDirect, Data in Brief (2020) |
| **Content** | ~700 commonly used materia medica |
| **Format** | Unicode Tibetan font, Wylie transliteration, phonetic transcription |

### TKDL Sowa Rigpa
| Category | Formulations |
|----------|--------------|
| Sowa Rigpa Formulations | 4,377 |

### Recognition
Government of India officially recognized Sowa Rigpa in 2010 as a "science of healing."

---

## Comparison Matrix

| Database | Plants | Compounds | Formulations | Targets | Open Access | Bulk Download | API |
|----------|--------|-----------|--------------|---------|-------------|---------------|-----|
| **IMPPAT 2.0** | 4,010 | 17,967 | - | Yes (27,365) | CC BY 4.0 | Export | No |
| **OSADHI** | 6,959 | 27,440 | - | No | Yes | Partial | No |
| **GRAYU** | 12,743 | 129,542 | 1,039 | Indirect | Research only | No | No |
| **TKDL** | - | - | 454,885+ | No | Restricted | No | No |
| **DHARA** | - | - | - | No | Yes | No | No |
| **NPACT** | - | 1,574 | - | Yes (1,980) | Yes | MOL files | No |
| **CMAUP** | 5,645 | 48,000 | - | Yes (646) | Yes | Yes | No |
| **NPASS** | - | 35,000+ | - | Yes (222,092) | Yes | Yes | No |
| **Dr. Duke's** | Many | Many | - | Limited | CC0 | CSV | No |
| **AromaDb** | 166 | 1,321 | - | Limited | Yes | 2D/3D | No |
| **sCentInDB** | 554 | 3,420 | - | No | Yes | CSV | No |
| **FRLHT/NMPB** | 7,263 | - | - | No | Registration | No | No |

---

## Recommendations for Knowledge Base Integration

### Priority 1: Primary Data Sources (Highest Value)
1. **IMPPAT 2.0**
   - Best combination of compound coverage, target data, and open licensing
   - Export to tab-separated format via web interface
   - Structure files in multiple formats

2. **NPASS**
   - Excellent for experimentally validated activity data
   - Free bulk download
   - Quantitative activity values (IC50, Ki, etc.)

3. **CMAUP**
   - Good global plant coverage including Indian species
   - Target and pathway data
   - Free download

### Priority 2: Supplementary Sources
4. **OSADHI** - Geographic and traditional knowledge context
5. **GRAYU** - Formulation-level data (unique value)
6. **NPACT** - Anti-cancer specific targets
7. **Dr. Duke's** - CC0 license, ethnobotanical context

### Priority 3: Reference Sources
8. **DHARA** - Literature discovery
9. **AYUSH Portal** - Research articles
10. **sCentInDB/AromaDb** - Essential oil-specific data

### Data Integration Strategy

1. **Core Compound-Target Data**
   - Start with IMPPAT 2.0 predicted targets
   - Supplement with NPACT validated cancer targets
   - Add NPASS quantitative activity data
   - Cross-reference with CMAUP pathway data

2. **Traditional Knowledge Layer**
   - Plant-compound associations from IMPPAT + OSADHI
   - Formulation data from GRAYU
   - Therapeutic use mappings

3. **Disease Ontology Mapping**
   - Leverage GRAYU's MeSH/DOID mappings
   - Connect traditional indications to modern disease terms

4. **Literature Support**
   - Index via DHARA and AYUSH Portal
   - Link compounds to supporting research

### Technical Considerations

- **No APIs available**: All databases require web scraping or manual export
- **Data harmonization needed**: Different identifier systems across databases
- **PubChem integration**: Most databases reference PubChem for compound standardization
- **STRING/ChEMBL linkage**: For target prediction validation

---

## References

1. Vivek-Ananth RP, et al. (2023). IMPPAT 2.0: An Enhanced and Expanded Phytochemical Atlas of Indian Medicinal Plants. ACS Omega, 8(9):8827-8845.
2. Mangal M, et al. (2013). NPACT: Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target database. Nucleic Acids Research, 41, D1124-D1129.
3. Zeng X, et al. (2019). CMAUP: a database of collective molecular activities of useful plants. Nucleic Acids Research, 47(D1), D1118-D1127.
4. Kumar Y, et al. (2018). AromaDb: A Database of Medicinal and Aromatic Plant's Aroma Molecules. Frontiers in Plant Science.
5. Baskaran SP, et al. (2025). sCentInDB: a database of essential oil chemical profiles of Indian medicinal plants. Molecular Diversity.

---

*Document generated: January 2026*
*Last updated: January 17, 2026*
