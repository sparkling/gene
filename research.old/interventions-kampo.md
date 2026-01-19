# Kampo Medicine Data Sources Research

**Research Date:** January 17, 2026
**Purpose:** Comprehensive inventory of databases and data sources for Japanese traditional medicine (Kampo) for a genetics/health knowledge base platform.

---

## Executive Summary

This document catalogs all available Kampo and related traditional medicine databases, including coverage statistics, access methods, data formats, and licensing information. Kampo medicine is a Japanese adaptation of traditional Chinese medicine (TCM), and thus many TCM databases contain relevant Kampo data.

### Key Findings

| Database | Records | Target Data | Download | API | License |
|----------|---------|-------------|----------|-----|---------|
| **KampoDB** | 298 formulas, 3,002 compounds | Yes (62,906 proteins) | Limited | No | CC BY-SA 4.0 |
| **KNApSAcK KAMPO** | 336 formulas, 278 plants | No | Partial | URL patterns | Non-commercial |
| **TradMPD** | 80 Kampo extracts, 80 crude drugs | No | No | No | Academic |
| **STORK** | 148 Kampo formulas | No | No | No | Open reference |
| **TM-MC 2.0** | 192 JP herbs, 34,107 compounds | Yes (13,992 targets) | Yes | No | Open |
| **BATMAN-TCM 2.0** | 54,832 formulas, 39,171 ingredients | Yes (2.3M TTIs) | Yes | Yes | Academic |
| **TCMSP** | 499 herbs, 29,384 compounds | Yes (3,311 targets) | Limited | No | ODbL |
| **HERB 2.0** | 7,263 herbs, 49,258 ingredients | Yes (12,933 targets) | Yes | No | Open |

---

## 1. Kampo-Specific Databases

### 1.1 KampoDB

**The primary Kampo-specific database with molecular target predictions.**

| Attribute | Details |
|-----------|---------|
| **URL** | https://wakanmoview.inm.u-toyama.ac.jp/kampo/ |
| **Maintainer** | Institute of Natural Medicine, University of Toyama |
| **Last Updated** | Active (current version) |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Kampo formulas | 298 |
| Crude drugs | 180 |
| Natural compounds | 3,002 |
| Proteins/genes | 62,906 |
| Docking simulation results | 3,063,505 |
| Known target proteins | 460 |
| Predicted target proteins | 1,369 |

#### Data Structure

KampoDB consists of three main components:

1. **Natural Medicines List**: Query by Kampo medicine name (e.g., "kakkonto")
2. **Functional Analysis**: Biological processes and pathways of target proteins
3. **Target Prediction**: Docking simulations and machine learning predictions

**Hierarchical Structure:**
- Layer 1: Kampo medicine
- Layer 2: Crude drugs (constituent herbs)
- Layer 3: Constituent compounds
- Layer 4: Target proteins

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search and browse functionality |
| Bulk download | Limited | Not explicitly documented |
| API | No | No programmatic access documented |

#### Data Format
- Web-based queries return structured data
- Compound IDs correspond to KNApSAcK IDs
- Chemical structures linked to KNApSAcK and PubChem

#### Licensing
**CC BY-SA 4.0** - All original data are available under Creative Commons Attribution-ShareAlike 4.0 license, permitting sharing and adaptation with appropriate attribution.

#### Gene/Protein Targets
**Yes** - Includes:
- 460 known target proteins (experimental evidence)
- 1,369 potential target proteins (predicted via docking/ML)
- Functional annotations for biological pathways
- Molecular function annotations

#### Key Publication
Sawada R, et al. (2018) "KampoDB, database of predicted targets and functional annotations of natural medicines." *Scientific Reports* 8:11216.

---

### 1.2 WAKANYAKU Database Portal (University of Toyama)

**Umbrella portal hosting multiple traditional medicine databases.**

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.inm.u-toyama.ac.jp/en/database/ |
| **Maintainer** | Institute of Natural Medicine, University of Toyama |

#### Component Databases

| Database | URL | Content |
|----------|-----|---------|
| TradMPD | https://dentomed.toyama-wakan.net/index_en | Genetic, chemical, biological data on crude drugs & Kampo |
| WAKANYAKU Wiki | https://www.inm.u-toyama.ac.jp/wiki/ | Natural medicine reviews, LC-MS profiles, Edo-period texts |
| ETHMEDmmm | https://ethmed.toyama-wakan.net/SearchEn/ | Ethnomedicine database, Museum specimens |
| KampoDB | https://wakanmoview.inm.u-toyama.ac.jp/kampo/ | Compound-protein interactions |
| Shouruihonzou DB | https://ethmed.toyama-wakan.net/honzou | Historical Chinese materia medica text |

---

### 1.3 TradMPD (Traditional Medical & Pharmaceutical Database)

**First repository of experimental data for Kampo natural medicines.**

| Attribute | Details |
|-----------|---------|
| **URL** | https://dentomed.toyama-wakan.net/index_en |
| **Maintainer** | Section of Pharmacognosy, Institute of Natural Medicine, University of Toyama |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Compounds (Wakanyaku Library) | ~80 |
| Crude drug hot water extracts | 80 |
| Kampo extracts | 80 |

#### Data Structure

**General Information:**
- Scientific and herbological information on crude drugs
- Genetic information on medicinal plants
- Medical information on Kampo formulas
- Literature searches on Kampo formulas
- Application of Kampo formulas for diseases
- Chemical compound data

**Analytical Results:**
- LC-MS profiling data
- Biological activity data for:
  - Neurodegenerative diseases
  - Allergy & inflammation
  - Cancer
  - Lifestyle diseases

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search functionality |
| Bulk download | No | Not available |
| API | No | Not documented |

#### Licensing
Academic/research use. Supported by Japanese Ministry of Education grants (2008-2012).

#### Gene/Protein Targets
**No** - Focuses on genetic authentication of crude drugs and chemical profiling, not target prediction.

---

### 1.4 STORK (Standards of Reporting Kampo Products)

**Official reference database for all 148 approved Kampo formulas in Japan.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://mpdb.nibiohn.go.jp/stork/ |
| **Maintainer** | National Institutes of Biomedical Innovation, Health and Nutrition (NIBIOHN) + NIHS |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Ethical Kampo drugs (prescription) | 148 |
| Total approved crude drugs | 241 |
| Crude drug preparations | 5 |

#### Data Structure

Provides standardized information on each of the 148 ethical Kampo drugs marketed in Japan:
- Formula composition
- Crude drug components with ratios
- Manufacturer information
- Package insert equivalent data (English)
- Glycyrrhiza dose information

#### Purpose
STORK corresponds to Item 5 (Intervention) of the CONSORT Statement 2010, providing uniform citations for Kampo products in clinical research articles.

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Browse/search all 148 formulas |
| Bulk download | No | Reference use only |
| API | No | Not available |

#### Licensing
**Open reference** - Free to use for citation in research articles.

#### Gene/Protein Targets
**No** - This is a product standardization database, not a molecular mechanisms database.

---

### 1.5 NIBIOHN Comprehensive Medicinal Plant Database (MPDB)

| Attribute | Details |
|-----------|---------|
| **URL** | http://mpdb.nibiohn.go.jp/ |
| **Maintainer** | Research Center for Medicinal Plant Resources, NIBIOHN |

#### Content
- Common names
- Crude drug information
- Plant species data
- Compound data
- Sample information
- List of Kampo resources

#### Access Method
Web-based interface. Note: Server may be unreliable (connection issues observed during research).

---

### 1.6 EKAT (Evidence Reports of Kampo Treatment)

**Clinical evidence database for Kampo medicines.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.jsom.or.jp/medical/ebm/ere/index.html |
| **Maintainer** | Japan Society for Oriental Medicine (JSOM) |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| RCT structured abstracts | 500+ |
| EKAT 2010 version | 345 RCTs |

#### Data Structure
Structured abstracts of randomized controlled trials (RCTs), each comprising 12 components:
- Available in Japanese, English, and Korean
- Linked to Cochrane Library (CENTRAL) as a Specialized Register

#### Access Method
Web-based access. Annual updates published on JSOM website.

#### Licensing
Free access for physicians and researchers.

#### Gene/Protein Targets
**No** - Clinical evidence database, not molecular mechanisms.

---

### 1.7 Metabolomics.jp Wiki - Kampo Section

**Wiki-based repository for crude drugs and Kampo medicine.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://metabolomics.jp/wiki/ |
| **Maintainer** | Research community (MediaWiki-based) |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Crude drugs | 158 |
| Kampo prescriptions | 348 |

#### Data Structure
- Taxonomic information
- Chemical information
- MediaWiki structure with SQL backend
- Japanese Pharmacopoeia entries on "Index:JPR" page

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Wiki browsing |
| MediaWiki API | Partial | Can download full page contents via MediaWiki API |
| SQL export | No | Backend is SQL but not directly accessible |

#### Licensing
Wiki-based, open access for viewing.

---

### 1.8 ETHMEDmmm (Data Base of Ethno-medicines in the world)

| Attribute | Details |
|-----------|---------|
| **URL** | https://ethmed.toyama-wakan.net/SearchEn/ |
| **Maintainer** | Museum of Materia Medica, Institute of Natural Medicine, University of Toyama |

#### Content
- ~30,000 herbal medicine specimens
- Crude drug specimens stored in the Museum
- Academic information on traditional medicine systems

#### Access Method
- Free consultation for reference
- **Downloading restricted** - Contents and photographs cannot be diverted without authorization

#### Licensing
**Restricted** - Must follow prescribed procedures for imagery data use.

---

## 2. KNApSAcK Family Databases

**Integrated metabolite-plant species databases including KAMPO section.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.knapsackfamily.com/knapsack_core/top.php |
| **Maintainer** | Nara Institute of Science and Technology (NAIST) |

### 2.1 KNApSAcK Core Database

| Data Type | Count |
|-----------|-------|
| Metabolite entries | 63,715 |
| Metabolite-species pairs | 159,095 |
| Species entries | 24,749 |

### 2.2 KNApSAcK KAMPO Database

| Data Type | Count |
|-----------|-------|
| Kampo formulae | 336 |
| Medicinal plants | 278 |

### 2.3 KNApSAcK WorldMap Database

| Data Type | Count |
|-----------|-------|
| GZ-plant pair entries | 41,548 |
| Geographic zones | 222 |
| Medicinal/edible plants | 15,240 |

### 2.4 Biological Activity Database

| Data Type | Count |
|-----------|-------|
| Biological activities | 2,418 |
| Plant-activity relationships | 33,706 |

#### Search Capabilities
- Organism names
- Metabolite identifiers
- Molecular formulas
- C_ID, CAS_ID, INCHI-KEY, INCHI-CODE, SMILES

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search functionality |
| URL pattern access | Yes | `info.php?sname=[type]&word=[query]` |
| Bulk download | Partial | Via Metabolomics Search Engine section |
| API | No | No formal REST API |

#### Licensing
**Non-commercial** - "Any content included in KNApSAcK database cannot be re-distributed or used for commercial purposes" without contacting administrators (skanaya@gtc.naist.jp).

#### Gene/Protein Targets
**No** - Focuses on metabolite-species relationships, not protein targets.

---

## 3. Traditional Chinese Medicine Databases (Kampo-Relevant)

Since Kampo is derived from Traditional Chinese Medicine (TCM), these databases contain relevant compound and target data.

### 3.1 BATMAN-TCM 2.0

**Most comprehensive TCM database with API access.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://bionet.ncpsb.org/batman-tcm |
| **Maintainer** | National Center for Protein Sciences Beijing |

#### Content Coverage (v2.0)

| Data Type | Count |
|-----------|-------|
| Formulae | 54,832 |
| Herbs | 8,404 |
| Ingredients | 39,171 |
| Known TTIs | 17,068 |
| Predicted TTIs | 2,319,272 |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Query by formula, herb, compound |
| Bulk download | Yes | All known/predicted TTIs downloadable |
| API | **Yes** | JSON format, URL-based parameters |

#### API Details
- Build URL with selected parameters
- Request type, output format, input item
- Returns JSON (computer-readable) or hypertext
- Response time: ~3 seconds per herb query

#### Data Format
- JSON (API)
- Downloadable datasets

#### Licensing
Academic/research use.

#### Gene/Protein Targets
**Yes** - Extensive target predictions:
- 17,068 known TCM ingredient-Target Interactions
- 2,319,272 high-confidence predicted TTIs
- KEGG/GO/disease enrichment analysis
- Network visualization

---

### 3.2 TCMSP (Traditional Chinese Medicine Systems Pharmacology)

| Attribute | Details |
|-----------|---------|
| **URL** | https://tcmsp-e.com/tcmsp.php |
| **Maintainer** | Northwest A&F University |
| **Version** | 2.3 |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Chinese herbs | 499 |
| Ingredients | 29,384 |
| Targets | 3,311 |
| Associated diseases | 837 |

#### ADME Properties Available
- Human oral bioavailability
- Half-life
- Drug-likeness
- Caco-2 permeability
- Blood-brain barrier
- Lipinski's rule of five
- Aqueous solubility

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Search by herb, chemical, InChIKey, CAS, target, disease |
| Bulk download | Limited | Not explicitly documented |
| API | No | Web interface only |

#### Licensing
**ODbL 1.0** (Open Database License) - Database structure
**DBCL 1.0** (Database Contents License) - Individual contents

#### Gene/Protein Targets
**Yes** - 3,311 targets with drug-target networks and drug-target-disease networks.

---

### 3.3 TCMBank

**Largest downloadable non-commercial TCM database.**

| Attribute | Details |
|-----------|---------|
| **URL** | https://TCMBank.cn/ |
| **Maintainer** | Artificial Intelligence Medical Center, Sun Yat-sen University |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Herbs | 9,192 |
| Ingredients | 61,966 |
| Targets | 15,179 |
| Diseases | 32,529 |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Search by herb, ingredient, target, disease |
| Bulk download | **Yes** | Largest downloadable TCM database |
| API endpoint | Partial | `/api` endpoint exists |

#### Data Format
- 3D structures in mol2 format
- Cross-reference links to CAS, DrugBank, PubChem, MeSH, OMIM, DO, ETCM, HERB

#### Licensing
**Non-commercial** - Free for research use.

#### Gene/Protein Targets
**Yes** - 15,179 targets with pathway and disease associations.

---

### 3.4 HERB 2.0

**High-throughput experiment-guided TCM database.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://herb.ac.cn/ |
| **Maintainer** | Research community |

#### Content Coverage (v2.0)

| Data Type | Count |
|-----------|-------|
| Herbs | 7,263 |
| Ingredients | 49,258 |
| Targets | 12,933 |
| Diseases | 28,212 |
| High-throughput experiments | 2,231 |
| Curated references | 6,644 |
| Clinical trials | 8,558 |
| Meta-analyses | 8,032 |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search and browse |
| Download | Yes | Custom scripts downloadable (probe2gene) |
| API | No | Not documented |

#### Licensing
Open access.

#### Gene/Protein Targets
**Yes** - 12,933 targets from manual curation and database integration.

---

### 3.5 SymMap 2.0

**TCM database with symptom mapping.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.symmap.org/ |
| **Maintainer** | Research community |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| TCM symptoms | 1,717 |
| Herbs | 499 |
| Modern medicine symptoms | 961 |
| Diseases | 5,235 |
| Herbal ingredients | 19,595 |
| Target genes | 4,302 |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search functionality |
| Download | Yes | Button at upper right on database interface |
| API | No | Not documented |

#### Licensing
Open access.

#### Gene/Protein Targets
**Yes** - 4,302 target genes from HIT, TCMSP, HPO, DrugBank, NCBI gene.

---

### 3.6 HIT 2.0 (Herbal Ingredients' Targets)

| Attribute | Details |
|-----------|---------|
| **URL** | http://hit2.badd-cao.net |
| **Maintainer** | Research community |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Compound-target activity pairs | 10,031 |
| Targets | 2,208 |
| Ingredients | 1,237 |
| Reputable herbs | 1,250+ |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Search by herb, ingredient, target |
| Download | Partial | "My-target" feature for personal curation |
| API | No | Not documented |

#### Crosslinks
TTD, DrugBank, KEGG, PDB, UniProt, Pfam, NCBI, TCM-ID

#### Gene/Protein Targets
**Yes** - Primary focus on herbal ingredient-target relationships.

---

### 3.7 TM-MC 2.0

**Northeast Asian traditional medicine database covering Korean, Chinese, and Japanese pharmacopoeias.**

| Attribute | Details |
|-----------|---------|
| **URL** | https://tm-mc.kr |
| **Maintainer** | Korean research community |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Medicinal materials | 635 |
| Chemical compounds | 34,107 (21,306 de-duplicated) |
| Targets | 13,992 |
| Diseases | 27,997 |
| Prescriptions | 5,075 (2,393 de-duplicated) |

#### Japanese Pharmacopoeia Coverage
- **192 medicinal materials** from Japanese Pharmacopoeia included
- Cross-referenced with Korean (454) and Chinese (556) pharmacopoeias

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search functionality |
| Download | **Yes** | Downloadable files available |
| API | No | Not documented |

#### Data Format
- Compounds deposited in PubChem (22,643)
- Identifiers and pharmacokinetic properties provided

#### Licensing
Open access.

#### Gene/Protein Targets
**Yes** - 13,992 targets with modern disease associations.

---

### 3.8 TCMID 2.0 (Traditional Chinese Medicine Integrated Database)

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.megabionet.org/tcmid/ |
| **Download** | http://47.100.169.139/tcmid/download/ |
| **Maintainer** | Research community |

#### Content Coverage

| Data Type | Count (Original) | Count (v2.0 additions) |
|-----------|------------------|------------------------|
| Prescriptions | ~47,000 | +15 |
| Herbs | 8,159 | - |
| Compounds | 25,210 | +18,203 |
| Drugs | 6,828 | +1,356 |
| Diseases | 3,791 | +842 |
| Related targets | 17,521 | +82 |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search |
| Bulk download | Yes | Download page available |
| API | No | Not documented |

#### Licensing
Free for research use. Citation requested for significant downloads.

#### Gene/Protein Targets
**Yes** - 17,521+ related targets with Drugbank/OMIM/PubChem links.

---

### 3.9 YaTCM

| Attribute | Details |
|-----------|---------|
| **URL** | http://cadd.pharmacy.nankai.edu.cn/yatcm/home |
| **Maintainer** | Nankai University |

#### Content
- Prescriptions, herbs, ingredients, diseases, targets, pathways
- Multi-voting similarity ensemble approach (MV-SEA) for target prediction

#### Tools
- Similarity search
- Substructure search
- Network analysis
- Pathway analysis

#### Gene/Protein Targets
**Yes** - Target prediction via MV-SEA methodology.

---

### 3.10 SuperTCM

| Attribute | Details |
|-----------|---------|
| **URL** | http://tcm.charite.de/supertcm |
| **Maintainer** | Charite - Universitatsmedizin Berlin |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| TCM drugs (herbs) | 6,516 |
| Botanical species | 5,372 |
| Active ingredients | 55,772 |
| Targets | 543 |
| KEGG pathways | 254 |
| Diseases | 8,634 |

#### Features
- KEGG Global Maps visualization
- Pathway projection in iPath3.0

#### Licensing
Open access.

#### Gene/Protein Targets
**Yes** - 543 targets in 254 KEGG pathways.

---

## 4. General Natural Products & Chemical Databases

### 4.1 STITCH

**Chemical-protein interaction database.**

| Attribute | Details |
|-----------|---------|
| **URL** | http://stitch.embl.de/ (Note: v5.0 from 2015, unsupported) |
| **Maintainer** | EMBL |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Organisms | 2,031 |
| Chemicals | 500,000 |
| Proteins | 9.6 million |
| Interactions | 1.6 billion |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Interactive network view |
| Download | Yes | Available under Creative Commons |
| API | Yes | Documented in Access section |

#### Licensing
**Creative Commons** - Separate commercial licensing for subset.

#### Gene/Protein Targets
**Yes** - Integrates data from ChEMBL, PDSP Ki Database, BindingDB, PDB, DrugBank, GLIDA, Matador, TTD, CTD.

---

### 4.2 COCONUT 2.0

**Largest open natural products database.**

| Attribute | Details |
|-----------|---------|
| **URL** | https://coconut.naturalproducts.net |
| **Download** | https://coconut.naturalproducts.net/download |
| **GitHub** | https://github.com/Steinbeck-Lab/coconut |
| **Maintainer** | Steinbeck Lab |

#### Content Coverage
- Integrates 63+ open NP resources
- Natural product structures with stereochemical forms
- Literature references
- Producing organisms
- Geographic distribution
- Precomputed molecular properties

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Browse and search |
| Bulk download | **Yes** | CSV, SDF, SQL dump formats |
| API | **Yes** | REST API (OpenAPI compliant) |

#### Download Formats
- CSV (general, mass spectrometry, molecular descriptors)
- SDF
- SQL dump
- Monthly releases on Zenodo

#### Licensing
**CC0** - Public domain, no restrictions.

#### Gene/Protein Targets
**Limited** - Focuses on chemical structures and properties, not targets.

---

### 4.3 PubChem

| Attribute | Details |
|-----------|---------|
| **URL** | https://pubchem.ncbi.nlm.nih.gov/ |
| **Downloads** | https://pubchem.ncbi.nlm.nih.gov/docs/downloads |
| **Maintainer** | NCBI/NIH |

#### Content
- World's largest collection of freely accessible chemical information
- Includes ~24,000 natural products (as of ChEMBL_02)

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search functionality |
| FTP bulk download | **Yes** | Millions of structures downloadable |
| API | **Yes** | Programmatic access available |

#### Data Format
- SDF
- JSON
- XML
- Various chemical notation formats

#### Licensing
**Public domain** - Free for any use.

---

### 4.4 ChEMBL

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/chembl/ |
| **API Docs** | https://chembl.gitbook.io/chembl-interface-documentation/ |
| **Maintainer** | EMBL-EBI |

#### Content Coverage

| Data Type | Count |
|-----------|-------|
| Compounds | 1.7 million |
| Activities | 14 million |
| Natural products | 24,000+ |

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Full search |
| FTP download | **Yes** | RDF, flatfiles available |
| REST API | **Yes** | Full programmatic access |
| Python client | Yes | With pagination and caching |

#### API Features
- Substructure search
- Similarity search
- SMILES, InChI Key, ChEMBL_ID support
- Pagination (up to 1000 per page)

#### Licensing
Open access for academic use.

---

### 4.5 DrugBank

| Attribute | Details |
|-----------|---------|
| **URL** | https://go.drugbank.com/ |
| **API Portal** | https://dev.drugbank.com/ |
| **Maintainer** | DrugBank Inc. |

#### Content
- 500,000+ drugs and drug products
- Herbs and Natural Products category (DBCAT003727)

#### Access Method

| Method | Available | Details |
|--------|-----------|---------|
| Web interface | Yes | Account required |
| Download | Partial | Academic license available |
| API | **Yes** | Clinical API, Discovery API |

#### Licensing
- Academic licenses for non-commercial research
- Commercial licenses available

---

## 5. Comparison Matrix

### Coverage Comparison

| Database | Formulas | Herbs | Compounds | Targets | Kampo-Specific |
|----------|----------|-------|-----------|---------|----------------|
| KampoDB | 298 | 180 | 3,002 | 62,906 | Yes |
| KNApSAcK KAMPO | 336 | 278 | - | - | Yes |
| STORK | 148 | - | - | - | Yes |
| TM-MC 2.0 | 5,075 | 635 | 34,107 | 13,992 | Partial (JP) |
| BATMAN-TCM 2.0 | 54,832 | 8,404 | 39,171 | 2.3M TTIs | No (TCM) |
| TCMSP | - | 499 | 29,384 | 3,311 | No (TCM) |
| TCMBank | - | 9,192 | 61,966 | 15,179 | No (TCM) |
| HERB 2.0 | - | 7,263 | 49,258 | 12,933 | No (TCM) |

### Access Method Comparison

| Database | Web | Download | API | Commercial OK |
|----------|-----|----------|-----|---------------|
| KampoDB | Yes | Limited | No | Yes (CC BY-SA) |
| KNApSAcK | Yes | Partial | URL patterns | No |
| TradMPD | Yes | No | No | Academic only |
| STORK | Yes | No | No | Yes |
| BATMAN-TCM 2.0 | Yes | Yes | **Yes** | Academic |
| TCMSP | Yes | Limited | No | Yes (ODbL) |
| TCMBank | Yes | **Yes** | Partial | Non-commercial |
| COCONUT | Yes | **Yes** | **Yes** | Yes (CC0) |
| ChEMBL | Yes | **Yes** | **Yes** | Academic |
| PubChem | Yes | **Yes** | **Yes** | Yes |

---

## 6. Recommendations for Integration

### Primary Kampo Data Sources

1. **KampoDB** - Best for Kampo-specific target predictions
   - Use for: Formula-compound-target relationships
   - Limitation: 298 formulas (not all 148 approved + others)
   - Integration: Web scraping or contact maintainers

2. **STORK** - Reference for all 148 approved Kampo formulas
   - Use for: Official formula compositions
   - Limitation: No molecular data
   - Integration: Web scraping

3. **TM-MC 2.0** - Japanese Pharmacopoeia herbs with targets
   - Use for: 192 JP herbs with compound/target data
   - Integration: Direct download available

### Supplementary TCM Data Sources

4. **BATMAN-TCM 2.0** - Best API access for target predictions
   - Use for: Large-scale ingredient-target interactions
   - Integration: REST API (JSON)

5. **TCMBank** - Largest downloadable database
   - Use for: Comprehensive herb/compound/target mappings
   - Integration: Bulk download

6. **COCONUT** - Open natural products
   - Use for: Chemical structure data with CC0 license
   - Integration: REST API + bulk downloads

### Data Integration Strategy

```
Layer 1: Kampo Formulas (STORK + KampoDB)
    |
Layer 2: Crude Drugs/Herbs (TM-MC 2.0 + KNApSAcK)
    |
Layer 3: Compounds (TCMBank + COCONUT + PubChem)
    |
Layer 4: Targets (BATMAN-TCM 2.0 + KampoDB + ChEMBL)
    |
Layer 5: Pathways/Diseases (KEGG + DisGeNET)
```

---

## 7. Technical Notes

### ID Mapping

| Source | Primary ID | Cross-references |
|--------|------------|------------------|
| KampoDB | Internal | KNApSAcK, PubChem |
| KNApSAcK | C_ID | CAS, InChI, SMILES |
| PubChem | CID | CAS, InChI, SMILES |
| ChEMBL | ChEMBL_ID | PubChem, InChI |
| UniProt | UniProt ID | NCBI Gene, Ensembl |

### Recommended Data Formats

- **Compounds**: SDF, SMILES, InChI (from PubChem/COCONUT)
- **Targets**: UniProt IDs (for gene integration)
- **Interactions**: JSON (from BATMAN-TCM API)
- **Pathways**: KEGG IDs

---

## 8. References

### Key Publications

1. Sawada R, et al. (2018) "KampoDB, database of predicted targets and functional annotations of natural medicines." *Scientific Reports* 8:11216.

2. Liu Z, et al. (2024) "BATMAN-TCM 2.0: an enhanced integrative database for known and predicted interactions between traditional Chinese medicine ingredients and target proteins." *Nucleic Acids Research* 52(D1):D1110-D1116.

3. Kim SK, et al. (2024) "TM-MC 2.0: an enhanced chemical database of medicinal materials in Northeast Asian traditional medicine." *BMC Complementary Medicine and Therapies* 24:37.

4. Fang S, et al. (2021) "HERB: a high-throughput experiment- and reference-guided database of traditional Chinese medicine." *Nucleic Acids Research* 49(D1):D1197-D1206.

5. Ru J, et al. (2014) "TCMSP: a database of systems pharmacology for drug discovery from herbal medicines." *Journal of Cheminformatics* 6:13.

### Database URLs Summary

| Database | URL |
|----------|-----|
| KampoDB | https://wakanmoview.inm.u-toyama.ac.jp/kampo/ |
| WAKANYAKU Portal | https://www.inm.u-toyama.ac.jp/en/database/ |
| TradMPD | https://dentomed.toyama-wakan.net/index_en |
| STORK | http://mpdb.nibiohn.go.jp/stork/ |
| KNApSAcK | http://www.knapsackfamily.com/knapsack_core/top.php |
| BATMAN-TCM 2.0 | http://bionet.ncpsb.org/batman-tcm |
| TCMSP | https://tcmsp-e.com/tcmsp.php |
| TCMBank | https://TCMBank.cn/ |
| HERB | http://herb.ac.cn/ |
| SymMap | http://www.symmap.org/ |
| TM-MC 2.0 | https://tm-mc.kr |
| COCONUT | https://coconut.naturalproducts.net |
| PubChem | https://pubchem.ncbi.nlm.nih.gov/ |
| ChEMBL | https://www.ebi.ac.uk/chembl/ |
| STITCH | http://stitch.embl.de/ |
| EKAT | http://www.jsom.or.jp/medical/ebm/ere/index.html |

---

*Document compiled from web research conducted January 2026. URLs and statistics should be verified before implementation as databases are regularly updated.*
