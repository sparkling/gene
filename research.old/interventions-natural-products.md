# Cross-System and Integrative Natural Product Databases

## Research Summary

This document catalogues databases that span multiple traditional medicine systems OR focus on natural products generally. These resources are essential for building a comprehensive genetics/health knowledge base that includes natural product interventions.

---

## Table of Contents

1. [Major Aggregated Natural Product Databases](#1-major-aggregated-natural-product-databases)
2. [Traditional Medicine System Databases](#2-traditional-medicine-system-databases)
3. [Food and Dietary Compound Databases](#3-food-and-dietary-compound-databases)
4. [Specialized Natural Product Databases](#4-specialized-natural-product-databases)
5. [Target Prediction and Network Pharmacology Tools](#5-target-prediction-and-network-pharmacology-tools)
6. [Bioactivity and Interaction Databases](#6-bioactivity-and-interaction-databases)
7. [ADMET and Toxicity Prediction Tools](#7-admet-and-toxicity-prediction-tools)
8. [Integration Recommendations](#8-integration-recommendations)

---

## 1. Major Aggregated Natural Product Databases

### 1.1 LOTUS (Natural Products Online)

| Attribute | Details |
|-----------|---------|
| **URL** | https://lotus.naturalproducts.net/ and https://www.wikidata.org/ |
| **Maintainer** | Wikidata community / Originally developed by Pierre-Marie Allard et al. |
| **Content** | 750,000+ referenced structure-organism pairs from comprehensive literature curation |
| **Systems Integrated** | Cross-system; captures natural products from all biological sources |
| **Access Method** | Wikidata SPARQL queries, web interface, bulk download |
| **Data Format** | Wikidata format (SPARQL), SDF, flat tables |
| **Licensing** | **CC0 (Creative Commons Zero)** - completely open |
| **Molecular Data** | SMILES, InChI, InChIKey, 2D structures |
| **Target Predictions** | No direct target predictions; links to other databases |

**Key Features:**
- First major natural products database hosted on Wikidata
- Community-curated with transparent provenance
- Referenced structure-organism pairs (not just structures)
- Dual hosting: Wikidata (community curation) + lotus.naturalproducts.net (user-friendly search)
- Note: The naturalproducts.net version is being phased out in favor of Wikidata

**Schema/Structure:**
- Structure (chemical entity) linked to Organism (species)
- Each pair linked to literature reference
- Taxonomic classification for organisms
- Chemical classification for structures

**Data Access Example (SPARQL):**
```sparql
SELECT ?compound ?compoundLabel ?organism ?organismLabel
WHERE {
  ?compound wdt:P235 ?inchikey .
  ?compound p:P703 ?statement .
  ?statement ps:P703 ?organism .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 100
```

---

### 1.2 COCONUT (Collection of Open Natural Products)

| Attribute | Details |
|-----------|---------|
| **URL** | https://coconut.naturalproducts.net/ |
| **Maintainer** | Steinbeck Lab, Friedrich Schiller University Jena, Germany |
| **Content** | 695,133 unique natural product structures (COCONUT 2.0, September 2024) |
| **Systems Integrated** | Aggregates 50+ open NP databases including UNPD, ZINC NP, TCM databases |
| **Access Method** | REST API, bulk download, web interface |
| **Data Format** | SDF (2D/3D), CSV, SQL dump (31.91 GB full database) |
| **Licensing** | **CC0 (Creative Commons Zero)** - no attribution required |
| **Molecular Data** | SMILES, InChI, InChIKey, 2D structures, 3D structures (RDKit generated) |
| **Target Predictions** | No direct targets; provides NP-likeness, synthetic accessibility, QED scores |

**Key Features:**
- Largest open-access aggregated NP database
- FAIR principles compliant
- Includes stereochemistry information (82,220 without stereocenters, 539,350 with preserved stereochemistry)
- REST API with OpenAPI specification
- Regular updates with versioning

**Download Options:**
| File Type | Size | Description |
|-----------|------|-------------|
| SDF 2D Lite | 287.6 MB | Basic 2D structures |
| SDF 2D Full | 691.7 MB | Full 2D with metadata |
| SDF 3D | 305.3 MB | RDKit-generated 3D coordinates |
| CSV Lite | 191 MB | Tabular format |
| CSV Full | 207.9 MB | Full metadata |
| SQL Dump | 31.91 GB | Complete database |

**API Documentation:** https://coconut.naturalproducts.net/api-documentation

**Schema Fields:**
- coconut_id, SMILES, InChI, InChIKey
- molecular_formula, molecular_weight
- np_likeness_score, synthetic_accessibility_score, qed_score
- data_sources (array of contributing databases)
- taxonomic information when available

---

### 1.3 NPASS (Natural Product Activity & Species Source)

| Attribute | Details |
|-----------|---------|
| **URL** | https://bidd.group/NPASS/ |
| **Maintainer** | BIDD Group (Bioinformatics and Drug Design), National University of Singapore |
| **Content** | 204,023 NPs, 8,764 targets, 48,940 species, 1,048,756 activity records |
| **Systems Integrated** | Cross-system; plants, microorganisms, marine organisms |
| **Access Method** | Web interface, download page |
| **Data Format** | Downloadable datasets (format unspecified in documentation) |
| **Licensing** | Free for academic use |
| **Molecular Data** | Chemical structures, calculated properties |
| **Target Predictions** | **Yes** - estimated activity profiles for ~66,600 NPs, plus experimental data |

**Key Features:**
- Focus on **quantitative experimental activity data** (IC50, Ki, EC50, GI50, MIC)
- Species source with collection location/date when available
- ADMET properties computed for all NPs
- Drug-likeness properties
- Composition/concentration data for NPs in species

**2023 Update Additions:**
- ~95,000 records of composition/concentration values
- Extended activity values for ~43,200 NPs against ~7,700 targets
- ~66,600 NPs with estimated activity profiles
- New species types: co-cultured microbes, engineered microbes

**Target Types:**
- Proteins and protein families
- Protein complexes
- Nucleic acids
- Subcellular targets
- Microbial/pathogenic organisms
- Cell lines
- Tissue targets

---

### 1.4 NPAtlas (Natural Products Atlas)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.npatlas.org/ |
| **Maintainer** | Linington Research Group, Simon Fraser University, Vancouver |
| **Content** | 36,545 compounds (as of latest release, database version 2024_09) |
| **Systems Integrated** | **Microbial only** - bacteria, fungi, cyanobacteria |
| **Access Method** | RESTful API, bulk download, web interface |
| **Data Format** | PostgreSQL (with RDKit cartridge), JSON API responses |
| **Licensing** | **CC BY 4.0** (Creative Commons Attribution 4.0) |
| **Molecular Data** | SMILES, InChI, InChIKey, structures |
| **Target Predictions** | Biosynthetic classification, chemical ontology |

**Key Features:**
- Comprehensive coverage of microbially-derived NPs from peer-reviewed literature
- Includes lichens and mushrooms (higher fungi)
- Molecular weight cutoff: 3000 Da
- Structure-based queries via RDKit PostgreSQL extension
- SNAP-MS tool for MS/MS-based compound family annotation
- Cross-links to other NP resources

**Exclusions:**
- Plant-derived compounds (unless also found in microbes)
- Marine invertebrates
- Diatoms
- Primary metabolites
- Patents and conference proceedings

**Data Updates:**
- 590 structural corrections incorporated
- 1,347 papers curated in current release
- Includes reassignments and revisions

---

### 1.5 SuperNatural 3.0

| Attribute | Details |
|-----------|---------|
| **URL** | http://bioinf-applied.charite.de/supernatural_3/ |
| **Maintainer** | Charite - Universitatsmedizin Berlin |
| **Content** | 449,058 natural compounds |
| **Systems Integrated** | Cross-system with vendor availability information |
| **Access Method** | Web interface (no login required) |
| **Data Format** | Web-based queries; download options unspecified |
| **Licensing** | Free for academic use |
| **Molecular Data** | Structural and physicochemical information |
| **Target Predictions** | **Yes** - mechanism of action, pathways, disease indications |

**Key Features:**
- Seven main function modules:
  1. Compound Search
  2. Mechanism of Actions
  3. Pathways
  4. Target Library
  5. Disease Indication
  6. COVID-19 virtual screening
  7. Organoleptic properties (taste prediction)

- Confidence scoring system:
  - Score 1.0: Has taxonomy info + linked to 3+ NP databases
  - Score 0.5: No taxonomy but linked to 1+ NP database

- Disease-specific drug-like space prediction:
  - Antiviral, antibacterial, antimalarial, anticancer
  - CNS-targeted compounds

---

### 1.6 UNPD (Universal Natural Products Database)

| Attribute | Details |
|-----------|---------|
| **URL** | Original site (pkuxxj.pku.edu.cn/UNPD) - **NO LONGER ACCESSIBLE** |
| **Maintainer** | Originally Peking University |
| **Content** | 229,358 natural product structures |
| **Status** | **Discontinued** - data preserved in other databases |
| **Access Method** | Via COCONUT or ISDB |

**Where to Access UNPD Data:**
1. **COCONUT**: Integrated 156,865 UNPD entries
2. **ISDB (In-Silico MS/MS Database)**: https://oolonek.github.io/ISDB/
   - Contains 170,602 NPs from UNPD with predicted MS/MS spectra
   - SDF version downloadable

---

## 2. Traditional Medicine System Databases

### 2.1 TCMSP (Traditional Chinese Medicine Systems Pharmacology)

| Attribute | Details |
|-----------|---------|
| **URL** | https://tcmsp-e.com/ |
| **Maintainer** | Lab of Systems Pharmacology, Northwest University, China |
| **Content** | 499 herbs, 29,384 ingredients, 3,311 targets, 837 diseases |
| **System** | Traditional Chinese Medicine (registered in Chinese Pharmacopoeia) |
| **Access Method** | Web interface, limited download |
| **Data Format** | Web queries; structured data export |
| **Licensing** | Open Database License |
| **Molecular Data** | SMILES, structures, ADME properties |
| **Target Predictions** | **Yes** - compound-target and target-disease networks |

**ADME Properties Provided (12 properties):**
- Human oral bioavailability (OB)
- Half-life
- Drug-likeness
- Caco-2 permeability
- Blood-brain barrier penetration
- Lipinski's Rule of Five compliance

**Network Features:**
- Herb-Compound-Target-Disease (H-C-T-D) networks
- Automatic network visualization
- Mechanism of action exploration

---

### 2.2 HERB (High-throughput Experiment- and Reference-guided Database of TCM)

| Attribute | Details |
|-----------|---------|
| **URL** | http://herb.ac.cn/ |
| **Maintainer** | Academic consortium (Chinese institutions) |
| **Content** | 7,263 herbs, 49,258 ingredients, 12,933 targets, 28,212 diseases |
| **System** | Traditional Chinese Medicine with modern medicine mapping |
| **Access Method** | Web interface |
| **Data Format** | Web-based queries |
| **Licensing** | Academic use |
| **Molecular Data** | Chemical structures from multiple TCM databases |
| **Target Predictions** | **Yes** - via CMap-style transcriptomic mapping |

**Unique Features:**
- Re-analyzed 6,164 gene expression profiles from 1,037 high-throughput experiments
- Maps TCM herbs/ingredients to 2,837 modern drugs via pharmacotranscriptomics
- 1,966 manually curated references
- Links TCM to CMap (Connectivity Map) drug database

**HERB 2.0 Updates:**
- 8,558 clinical trials and 8,032 meta-analyses
- 2,231 high-throughput experiments
- 6,644 curated references
- Knowledge graph representations
- Auto-analysis interface for user gene expression profiles

---

### 2.3 SymMap

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.symmap.org/ |
| **Maintainer** | Chinese academic institutions |
| **Content** | 499 herbs, 1,717 TCM symptoms, 961 MM symptoms, 5,235 diseases, 19,595 ingredients, 4,302 targets |
| **System** | TCM with Modern Medicine symptom mapping |
| **Access Method** | Web interface |
| **Data Format** | Web-based queries |
| **Licensing** | Academic use |
| **Molecular Data** | Via linked databases |
| **Target Predictions** | Indirect via network analysis |

**Unique Features:**
- Maps TCM symptoms to UMLS (Unified Medical Language System)
- Expert consensus from 17 leading TCM practitioners
- Phenotypic AND molecular level integration
- Heterogeneous network of all components

---

### 2.4 BATMAN-TCM 2.0

| Attribute | Details |
|-----------|---------|
| **URL** | http://bionet.ncpsb.org/batman-tcm |
| **Maintainer** | Beijing Proteome Research Center |
| **Content** | 17,068 known TTIs + ~2.3 million predicted TTIs |
| **System** | Traditional Chinese Medicine |
| **Access Method** | Web interface with analysis tools |
| **Data Format** | Analysis results, network visualizations |
| **Licensing** | Academic use |
| **Molecular Data** | Ingredient structures |
| **Target Predictions** | **Yes** - similarity-based target prediction |

**Analysis Functions:**
1. TCM ingredients' target prediction
2. Pathway enrichment (KEGG)
3. Gene Ontology functional term enrichment
4. Disease enrichment
5. Network visualization
6. Multi-TCM comparison

**Performance:**
- Query time reduced from 80 seconds to 3 seconds (v2.0)
- Parallelization with R: doParallel
- In-memory databases for fast access

---

### 2.5 IMPPAT 2.0 (Indian Medicinal Plants, Phytochemistry And Therapeutics)

| Attribute | Details |
|-----------|---------|
| **URL** | https://cb.imsc.res.in/imppat/ |
| **Maintainer** | Institute of Mathematical Sciences (IMSc), Chennai, India |
| **Content** | 4,010 plants, 17,967 phytochemicals, 1,095 therapeutic uses |
| **System** | Indian Traditional Medicine (Ayurveda, Siddha, Unani) |
| **Access Method** | Web interface, download |
| **Data Format** | 2D/3D structures downloadable |
| **Licensing** | Academic use |
| **Molecular Data** | 2D/3D structures, SMILES, physicochemical properties |
| **Target Predictions** | Drug-likeness scores, predicted ADMET properties |

**Key Features:**
- Largest digital database of Indian medicinal plants
- Phytochemicals provided at plant part level
- Multiple drug-likeness scoring schemes
- ADMET predictions computed
- Used for COVID-19 drug discovery studies

**Source Material:**
- 100+ books on traditional Indian medicine
- 7,000+ published research articles
- Manual curation and digitization

---

### 2.6 CMAUP (Collective Molecular Activities of Useful Plants)

| Attribute | Details |
|-----------|---------|
| **URL** | https://bidd.group/CMAUP/ |
| **Maintainer** | BIDD Group, National University of Singapore |
| **Content** | 5,765 plants, 47,645 ingredients, 646 targets, 656 diseases |
| **Systems Integrated** | Cross-system - 153 countries/regions |
| **Access Method** | Web interface, searchable by geography |
| **Data Format** | Web-based queries |
| **Licensing** | Free for academic use |
| **Molecular Data** | Ingredient structures |
| **Target Predictions** | **Yes** - target heatmaps, pathway associations |

**Plant Categories:**
- 2,567 medicinal plants
- 170 food plants
- 1,567 edible plants
- 3 agricultural plants
- 119 garden plants

**2024 Update Features:**
- Human transcriptomic changes overlapping with 1,152 targets
- 74 diseases from 20,027 patient samples
- Clinical information for 185 plants in 691 clinical trials
- Drug development info for 4,694 drug-producing plants

**Visualization:**
- Interactive world map for plant distribution
- 2D target-ingredient heatmaps
- KEGG pathway highlighting

---

### 2.7 HIT 2.0 (Herbal Ingredients' Targets)

| Attribute | Details |
|-----------|---------|
| **URL** | http://hit2.badd-cao.net |
| **Maintainer** | Chinese academic institutions |
| **Content** | 1,237 ingredients, 2,208 targets, 10,031 compound-target activity pairs |
| **System** | Cross-system herbal medicine |
| **Access Method** | Web interface with auto-mining |
| **Data Format** | Web-based; manual curation tools |
| **Licensing** | Academic use |
| **Molecular Data** | Ingredient structures |
| **Target Predictions** | Curated from literature |

**Unique Features:**
- Auto-mining from daily PubMed releases
- "My-target" system for personal curation
- Cross-links to TTD, DrugBank, KEGG, PDB, UniProt

**Target Types Covered:**
- Directly activated/inhibited genes/proteins
- Protein binders
- Enzyme substrates/products
- Regulated genes under treatment

---

### 2.8 KNApSAcK Family

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.knapsackfamily.com/KNApSAcK/ |
| **Maintainer** | Nara Institute of Science and Technology (NAIST), Japan |
| **Content** | 50,048 metabolites, 20,741 species, 101,500 species-metabolite relationships |
| **Systems Integrated** | Plants, microorganisms, traditional medicines (KAMPO, JAMU) |
| **Access Method** | Web interface, database download |
| **Data Format** | Various (downloadable) |
| **Licensing** | Academic use |
| **Molecular Data** | Structures, mass spectra (multiple ionization modes) |
| **Target Predictions** | Biological activity relationships |

**Database Family Components:**
| Database | Content |
|----------|---------|
| KNApSAcK Core | 101,500 species-metabolite relationships |
| KNApSAcK WorldMap | 41,548 geographic zone-plant pairs |
| KAMPO DB | 336 Japanese herbal formulae, 278 plants |
| JAMU DB | 5,310 Indonesian formulae, 550 plants |
| Biological Activity DB | 2,418 activities, 33,706 plant-activity relationships |

**Mass Spectrometry Features:**
- Search by accurate mass
- Multiple ionization modes
- Metabolomics research support

---

## 3. Food and Dietary Compound Databases

### 3.1 FooDB

| Attribute | Details |
|-----------|---------|
| **URL** | https://foodb.ca/ |
| **Maintainer** | The Metabolomics Innovation Centre (TMIC), University of Alberta |
| **Content** | World's largest food constituent database |
| **Systems Integrated** | Food chemistry across all food sources |
| **Access Method** | REST API (beta), bulk download |
| **Data Format** | Multiple formats via downloads |
| **Licensing** | Free; citation required for publications |
| **Molecular Data** | Chemical structures, physicochemical properties |
| **Target Predictions** | No direct targets |

**Content Categories:**
- Macronutrients
- Micronutrients
- Flavor compounds
- Color compounds
- Aroma compounds
- Texture compounds

**Access Points:**
- Food Browse (by food source)
- Compound Browse (by chemical)
- API: https://foodb.ca/api_doc
- Downloads: https://foodb.ca/downloads

---

### 3.2 Phenol-Explorer

| Attribute | Details |
|-----------|---------|
| **URL** | http://phenol-explorer.eu/ |
| **Maintainer** | INRA (French National Institute for Agricultural Research) and collaborators |
| **Content** | 35,000+ content values, 500 polyphenols, 400+ foods |
| **Systems Integrated** | Food polyphenols globally |
| **Access Method** | Web interface, advanced search |
| **Data Format** | Web queries |
| **Licensing** | Free academic access |
| **Molecular Data** | Chemical structures |
| **Target Predictions** | Metabolism and pharmacokinetics data |

**Unique Features:**
- First comprehensive polyphenol food content database
- 60,000+ original values from 1,300+ publications
- Expert evaluation before inclusion

**Version History:**
| Version | Feature |
|---------|---------|
| 1.0 (2009) | Initial polyphenol content data |
| 2.0 | 380 metabolites, pharmacokinetics, 236 publications |
| 3.0 | Food processing effects on 161 polyphenols |
| 3.6 (2015) | 1,451 new lignan content values |

**Data Types:**
- Polyphenol content in raw foods
- Metabolite identification in biofluids
- Pharmacokinetic parameters
- Processing effects

---

## 4. Specialized Natural Product Databases

### 4.1 NuBBEDB (Brazilian Biodiversity)

| Attribute | Details |
|-----------|---------|
| **URL** | https://nubbe.iq.unesp.br/portal/nubbedb.html |
| **Maintainer** | NuBBE (Nuclei of Bioassays, Ecophysiology and Biosynthesis of Natural Products), UNESP, Brazil |
| **Content** | 2,147 compounds (1,688 from plants, 325 semi-synthetic, 109 from microorganisms) |
| **Systems Integrated** | Brazilian biodiversity focus |
| **Access Method** | Web interface |
| **Data Format** | Web-based, NMR spectroscopic data |
| **Licensing** | Free academic access |
| **Molecular Data** | Structures, NMR data, molecular descriptors |
| **Target Predictions** | Pharmacological properties |

**Content Breakdown:**
- 78% plant-derived
- 15% semi-synthetic
- 5% microbial
- 1.6% biotransformation products
- 0.2% marine

**Unique Data:**
- NMR spectroscopic data
- Species geographic locations
- Toxicological information

---

### 4.2 NPCARE (Natural Products for Cancer Regulation)

| Attribute | Details |
|-----------|---------|
| **URL** | http://silver.sejong.ac.kr/npcare |
| **Maintainer** | Sejong University, Korea |
| **Content** | 6,578 compounds, 2,566 fractional extracts, 1,952 species, 1,107 cell lines, 34 cancer types |
| **Systems Integrated** | Plants, marine organisms, fungi, bacteria |
| **Access Method** | Web interface |
| **Data Format** | Web queries with hyperlinks to external DBs |
| **Licensing** | Free academic access |
| **Molecular Data** | PubChem IDs, structures |
| **Target Predictions** | Cancer-related gene/protein targets |

**Anticancer Activity Criteria:**
1. Growth inhibition of cancer cell lines
2. Downregulation of oncogenes
3. Upregulation of tumor suppressor genes

**External Links:**
- HUGO Gene Nomenclature Committee
- Ensembl Genome Browser
- UniProt
- COSMIC (Catalogue of Somatic Mutations in Cancer)
- PubChem

---

### 4.3 GNPS (Global Natural Products Social Molecular Networking)

| Attribute | Details |
|-----------|---------|
| **URL** | https://gnps.ucsd.edu/ |
| **Maintainer** | Dorrestein Lab, UCSD |
| **Content** | 1,800+ public datasets, 490,000+ MS files, 1.2 billion tandem mass spectra |
| **Systems Integrated** | All natural product sources (MS/MS focus) |
| **Access Method** | Web interface, data upload/download, API |
| **Data Format** | MS/MS spectra, MGF, mzML |
| **Licensing** | Open access with sharing requirements |
| **Molecular Data** | MS/MS spectra, molecular networks |
| **Target Predictions** | Structure annotation via spectral matching |

**Core Functions:**
1. Molecular networking
2. Spectral library search
3. Community-curated reference spectra
4. MASST search (search against all public datasets)
5. DEREPLICATOR for peptidic NPs

**Integration:**
- MassIVE data repository
- GNPS spectral libraries
- Cross-platform annotation

---

### 4.4 MIBiG (Minimum Information about a Biosynthetic Gene cluster)

| Attribute | Details |
|-----------|---------|
| **URL** | https://mibig.secondarymetabolites.org/ |
| **Maintainer** | International consortium |
| **Content** | Experimentally validated biosynthetic gene clusters |
| **Systems Integrated** | Bacteria, fungi, plants |
| **Access Method** | Web interface, bulk download, API |
| **Data Format** | JSON, standardized annotations |
| **Licensing** | Open access |
| **Molecular Data** | Product structures, cross-links to NP Atlas, PubChem |
| **Target Predictions** | Biosynthetic pathway predictions |

**Version 3.0 Features:**
- 661 new entries
- New structures and bioactivities
- Enzyme domain substrate specificities
- Standardized BGC annotations

**Data Fields:**
- Genomic locus coordinates
- Producing organism taxonomy
- Biosynthetic class
- Compound names
- Literature references
- Gene functions
- Monomer identity

---

### 4.5 Dictionary of Natural Products (DNP)

| Attribute | Details |
|-----------|---------|
| **URL** | https://dnp.chemnetbase.com/ |
| **Maintainer** | CRC Press / Taylor & Francis |
| **Content** | 328,000+ compounds |
| **Systems Integrated** | All natural product sources |
| **Access Method** | Subscription-based web interface |
| **Data Format** | ASCII/Oracle files, online access |
| **Licensing** | **Commercial subscription required** |
| **Molecular Data** | Comprehensive structural data |
| **Target Predictions** | Limited |

**Content Categories:**
- Alkaloids, Antibiotics, Flavonoids
- Carbohydrates, Tannins, Terpenoids
- Steroids, Polyketides, Polypyrroles
- Lignans, Amino acids, Peptides

**Updates:** ~10,000 new entries per year

**Note:** This is a **commercial database** - not suitable for open knowledge base without licensing agreement.

---

### 4.6 ZINC Natural Products Subset

| Attribute | Details |
|-----------|---------|
| **URL** | https://zinc15.docking.org/substances/subsets/natural-products/ |
| **Maintainer** | Shoichet Lab, UCSF |
| **Content** | Millions of compounds (subset of full ZINC database) |
| **Systems Integrated** | Commercially available NPs |
| **Access Method** | Web download, scripted access |
| **Data Format** | SMILES, SDF, mol2 |
| **Licensing** | Free for research |
| **Molecular Data** | 2D/3D structures, docking-ready formats |
| **Target Predictions** | No direct predictions |

**Features:**
- Purchasability information
- Pre-filtered subsets (lead-like, drug-like)
- Ready for virtual screening
- Download via scripts for large datasets

---

## 5. Target Prediction and Network Pharmacology Tools

### 5.1 SwissTargetPrediction

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.swisstargetprediction.ch |
| **Maintainer** | SIB Swiss Institute of Bioinformatics |
| **Method** | 2D/3D molecular similarity |
| **Database** | 376,342 active compounds, 3,068 targets |
| **Species** | Human and other vertebrates |
| **Access** | Web interface (no API for batch processing) |
| **Output** | Ranked target predictions with probability scores |
| **Runtime** | 15-20 seconds per drug-like molecule |

**Input Methods:**
- SMILES string
- MarvinJS structure sketcher
- File upload

**Performance:**
- >70% accuracy for at least one correct human target in top 15

**Interoperability:**
- Links to other SIB tools (SwissADME, SwissSimilarity, etc.)

---

### 5.2 SEA (Similarity Ensemble Approach)

| Attribute | Details |
|-----------|---------|
| **URL** | https://sea.bkslab.org/ |
| **Maintainer** | Shoichet Lab, UCSF |
| **Method** | Set-wise chemical similarity among ligands |
| **Access** | Web interface |
| **Validation** | Only target prediction model with systematic experimental validation |

**Methodology:**
- Relates proteins based on ligand set similarity
- Uses various molecular fingerprints
- Best performance with Atom pair and Topological fingerprints

**Fingerprint Performance (at significance 0.05):**
| Fingerprint | Precision |
|-------------|-----------|
| Atom pair | Highest |
| Topological | 83.7% |
| Morgan | Good |
| MACCS | Moderate |
| Pharmacophore | Lowest |

**Limitations:**
- Less effective for compounds with unique/unknown structures
- Requires structural similarity to known compounds

---

### 5.3 PharmMapper

| Attribute | Details |
|-----------|---------|
| **URL** | http://lilab.ecust.edu.cn/pharmmapper/ |
| **Maintainer** | Li Lab, East China University of Science and Technology |
| **Method** | Reverse pharmacophore mapping |
| **Database** | 23,236 proteins (16,159 druggable + 51,431 ligandable models) |
| **Access** | Web interface (free, no login) |
| **Runtime** | ~1 hour for full database screen |

**Input:**
- Mol2 file with 3D coordinates (required)
- Optional parameters for speed/accuracy control

**Target Sources:**
- TargetBank
- BindingDB
- DrugBank
- Potential drug target databases

**Applications:**
- Drug repurposing
- Side effect prediction
- Target identification for natural products

---

### 5.4 Cytoscape

| Attribute | Details |
|-----------|---------|
| **URL** | https://cytoscape.org/ |
| **Type** | Network visualization and analysis platform |
| **Access** | Free download (desktop application) |
| **Licensing** | LGPL (open source) |

**Network Pharmacology Features:**
- Protein-protein interaction networks
- Drug-target interaction visualization
- Pathway enrichment visualization
- Multi-omics data integration

**Key Plugins/Apps:**
| App | Function |
|-----|----------|
| stringApp | STRING network import |
| DrugViz | Small molecule visualization |
| ClueGO | GO term analysis |
| CyTargetLinker | Drug-target network extension |

**File Formats:**
- GraphML, XGMML, SIF
- Import from igraph, Pajek, GraphViz

---

### 5.5 STRING

| Attribute | Details |
|-----------|---------|
| **URL** | https://string-db.org/ |
| **Type** | Protein-protein interaction database |
| **Content** | 24.6 million proteins, 3.1 billion interactions |
| **Access** | Web interface, API, Cytoscape app |
| **Licensing** | CC BY 4.0 |

**Applications in Network Pharmacology:**
- PPI network construction for drug targets
- Functional enrichment analysis
- Network visualization

---

## 6. Bioactivity and Interaction Databases

### 6.1 ChEMBL

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/chembl/ |
| **Maintainer** | EMBL-EBI |
| **Content** | Large-scale bioactivity database |
| **Access Method** | Web interface, REST API, bulk download |
| **Data Format** | SQL, SDF, TSV |
| **Licensing** | CC BY-SA 3.0 |
| **Molecular Data** | SMILES, InChI, structures |
| **Target Predictions** | Experimental bioactivity (Ki, Kd, IC50, EC50) |

**Natural Products Features:**
- NP flag for ~64,000 molecules (release 33)
- NP-likeness scores
- Mapping to COCONUT database
- Stereochemistry-independent NP flagging

**API Features:**
- Full programmatic access
- Multiple data types (molecules, targets, activities, documents)
- Integration-ready for external applications

**Training Course:** EMBL-EBI offers API access training

---

### 6.2 BindingDB

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.bindingdb.org/ |
| **Maintainer** | BindingDB Team (UCSD) |
| **Content** | 3.2M data points, 1.4M compounds, 11.4K targets |
| **Access Method** | Web interface, RESTful API, bulk download, KNIME workflows |
| **Data Format** | TSV, SDF, FASTA |
| **Licensing** | Free for academic use |
| **Molecular Data** | SMILES, structures |
| **Binding Data** | Experimental Ki, Kd, IC50, EC50 |

**2024 FAIR Update:**
- Publication and curation dates for each measurement
- Improved AI/ML training support

**Download Files:**
| File | Description |
|------|-------------|
| BindingDBTargetSequences.fasta | Protein sequences (7.23 MB) |
| BindingDB_CID.txt | PubChem CID mapping (22.70 MB) |
| BindingDB_UniProt.txt | UniProt mapping |
| BindingDB_DrugBankID.txt | DrugBank mapping |

**Browser Extension:**
- Chrome, Firefox, Edge, Brave
- Flags articles with BindingDB data
- Direct download to TSV

---

### 6.3 STITCH

| Attribute | Details |
|-----------|---------|
| **URL** | http://stitch.embl.de/ |
| **Maintainer** | EMBL (part of STRING family) |
| **Content** | 500K chemicals, 9.6M proteins, 1.6B interactions |
| **Access Method** | Web interface, API, bulk download |
| **Data Format** | Network files, structured data |
| **Licensing** | CC BY 4.0 |
| **Interaction Data** | Experimental, curated, text-mined, predicted |

**Data Channels:**
1. **Experiments**: ChEMBL (Ki, IC50), PDSP Ki, PDB
2. **Databases**: DrugBank, GLIDA, Matador, TTD, CTD, KEGG, Reactome
3. **Text mining**: Literature extraction
4. **Predictions**: Structure-based inference

**Features:**
- Tissue filtering
- Binding affinity visualization
- Multi-organism support (1,133 organisms)

---

### 6.4 ChEBI (Chemical Entities of Biological Interest)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/chebi/ |
| **Maintainer** | EMBL-EBI |
| **Content** | 195,000+ entries |
| **Access Method** | Web interface, web services, bulk download, libChEBI API |
| **Data Format** | Relational tables, flat files, SDF |
| **Licensing** | CC BY 4.0 (non-proprietary) |
| **Molecular Data** | Structures, InChI, SMILES |
| **Ontology** | Part of Open Biomedical Ontologies (OBO) |

**libChEBI API Features:**
- Python, Java implementations
- Automatic monthly updates
- Simplified interface to ChEBI data

**Search Options:**
- Text search
- Structure search (InChI, SMILES)
- Molecular formula
- Ontology navigation

---

### 6.5 DGIdb (Drug-Gene Interaction Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://dgidb.org/ |
| **Maintainer** | Washington University School of Medicine |
| **Content** | Drug-gene interactions from 40+ sources |
| **Access Method** | Web interface, GraphQL API, bulk download |
| **Data Format** | TSV, JSON (GraphQL) |
| **Licensing** | Open access |

**DGIdb 5.0 Features:**
- GraphQL API (replaced REST)
- Knowledge graph representations
- API playground for query testing

**Applications:**
- Druggable genome search
- Hypothesis generation
- Precision medicine support

---

### 6.6 TTD (Therapeutic Target Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://idrblab.net/ttd/ |
| **Maintainer** | IDRB, Zhejiang University |
| **Content** | 3,730 targets, 39,862 drugs |
| **Access Method** | Web interface, download |
| **Data Format** | Structured downloads |
| **Licensing** | Free (no login required) |

**Target Categories:**
- 532 successful targets
- 1,442 clinical trial targets
- 239 preclinical/patented targets
- 1,517 research targets

**Drug Categories:**
- 2,895 approved
- 11,796 clinical trial
- 5,041 preclinical/patented
- 20,130 experimental

**Natural Products:**
- "Nature-Derived Drugs" section
- Species origins and families
- Cross-links to ICD disease codes

---

## 7. ADMET and Toxicity Prediction Tools

### 7.1 ProTox 3.0

| Attribute | Details |
|-----------|---------|
| **URL** | https://tox.charite.de/ |
| **Maintainer** | Charite - Universitatsmedizin Berlin |
| **Models** | 61 toxicity prediction models |
| **Access** | Web interface (free, no login) |
| **Input** | SMILES or PubChem name |

**Toxicity Endpoints:**
| Category | Endpoints |
|----------|-----------|
| Oral Toxicity | 6 toxicity classes |
| Organ Toxicity | Hepato-, neuro-, respiratory, cardio-, nephrotoxicity |
| Toxicological Endpoints | Mutagenicity, carcinogenicity, cytotoxicity, immunotoxicity |
| Toxicological Pathways | 12 adverse outcome pathways (AOPs) |
| Toxicity Targets | 15 targets |
| MIEs | 14 molecular initiating events |
| Metabolism | 6 molecular targets |

**Output:**
- Confidence scores per endpoint
- Toxicity radar plot
- Network plot visualization

---

### 7.2 pkCSM

| Attribute | Details |
|-----------|---------|
| **URL** | https://biosig.lab.uq.edu.au/pkcsm/ |
| **Maintainer** | Biosig Lab, University of Queensland |
| **Method** | Graph-based molecular signatures |
| **Access** | Web interface (free) |
| **Privacy** | No data retention |

**Models:**
| Type | Count | Examples |
|------|-------|----------|
| Regression | 14 | Water solubility, Caco-2 permeability |
| Classification | 16 | AMES test, hERG inhibition, hepatotoxicity |

**Performance:**
- AMES test: 83.8% accuracy (vs. ToxTree: 75.8%)
- Handles datasets up to 18,000 compounds

**ADMET Categories:**
- **A**bsorption: Intestinal, Caco-2, P-gp substrate
- **D**istribution: VDss, BBB penetration, CNS permeability
- **M**etabolism: CYP450 inhibition/substrate
- **E**xcretion: Total clearance, renal OCT2 substrate
- **T**oxicity: AMES, hERG, hepatotoxicity, skin sensitization

---

### 7.3 admetSAR 3.0

| Attribute | Details |
|-----------|---------|
| **URL** | http://lmmd.ecust.edu.cn/admetsar3/ |
| **Maintainer** | East China University of Science and Technology |
| **Content** | 370,000+ experimental ADMET data, 104,652 compounds |
| **Models** | 119 prediction endpoints |
| **Access** | Web interface |

**Modules:**
1. **Search**: Experimental ADMET data + similarity search
2. **Prediction**: 119 endpoints including environmental and cosmetic risk
3. **ADMETopt**: Automatic ADMET optimization via transformation rules

**Architecture:**
- Multi-task graph neural network framework
- Read-across support via structure similarity

**Input Options:**
- SMILES string
- Structure drawing
- Batch file upload

---

### 7.4 ADMETlab 3.0

| Attribute | Details |
|-----------|---------|
| **URL** | https://admetmesh.scbdd.com/ |
| **Maintainer** | SCBDD (Computational Biology & Drug Design) |
| **Access** | Web interface, **API available** |
| **Features** | Decision support system |

**Key Features:**
- Broader coverage than previous versions
- Improved prediction performance
- API functionality for batch processing
- Decision support for compound optimization

---

## 8. Integration Recommendations

### 8.1 Priority Databases for Knowledge Base

**Tier 1 - Essential (Open, Comprehensive, Well-Maintained):**

| Database | Reason |
|----------|--------|
| **COCONUT 2.0** | Largest open aggregated NP database, CC0 license, REST API |
| **LOTUS** | Wikidata-hosted, CC0, community-curated, structure-organism pairs |
| **ChEMBL** | Gold standard bioactivity, REST API, NP flagging |
| **NPAtlas** | Comprehensive microbial NPs, CC BY 4.0, API |
| **NPASS** | Quantitative activity data with species sources |

**Tier 2 - Important (Specialized or System-Specific):**

| Database | Use Case |
|----------|----------|
| **TCMSP** | TCM herbs with ADME properties |
| **HERB 2.0** | TCM with transcriptomic validation |
| **IMPPAT 2.0** | Indian traditional medicine |
| **BindingDB** | Experimental binding affinities |
| **FooDB** | Food-derived compounds |

**Tier 3 - Supplementary:**

| Database | Use Case |
|----------|----------|
| **Phenol-Explorer** | Polyphenol content and metabolism |
| **SuperNatural 3.0** | Mechanism of action, taste profiles |
| **CMAUP** | Plant geographic distribution |
| **KNApSAcK** | Mass spectrometry support |

### 8.2 Target Prediction Pipeline

**Recommended Workflow:**

```
Input: Natural Product Structure (SMILES/InChI)
            |
            v
    +-------+-------+
    |               |
    v               v
SwissTarget    PharmMapper
(2D/3D sim)    (3D pharmacophore)
    |               |
    +-------+-------+
            |
            v
    Consensus Targets
            |
            v
    +-------+-------+-------+
    |       |       |       |
    v       v       v       v
  STRING  STITCH  ChEMBL  BindingDB
  (PPI)   (chem)  (activity) (Kd/Ki)
    |       |       |       |
    +-------+-------+-------+
            |
            v
    Network Analysis (Cytoscape)
            |
            v
    Pathway Enrichment (KEGG, GO)
```

### 8.3 Data Integration Schema

**Proposed Unified Natural Product Schema:**

```json
{
  "compound": {
    "id": "string (internal)",
    "external_ids": {
      "coconut": "string",
      "lotus_wikidata": "string",
      "pubchem_cid": "integer",
      "chembl_id": "string",
      "inchikey": "string"
    },
    "structure": {
      "smiles": "string",
      "inchi": "string",
      "inchikey": "string",
      "mol_formula": "string",
      "mol_weight": "float"
    },
    "properties": {
      "np_likeness": "float",
      "drug_likeness": "float",
      "synthetic_accessibility": "float",
      "qed_score": "float"
    }
  },
  "source_organisms": [
    {
      "taxon_id": "integer (NCBI)",
      "species_name": "string",
      "common_names": ["string"],
      "part": "string (optional)",
      "geographic_region": "string (optional)"
    }
  ],
  "traditional_medicine": {
    "systems": ["TCM", "Ayurveda", "KAMPO", "etc."],
    "therapeutic_uses": ["string"],
    "preparation_methods": ["string"]
  },
  "bioactivity": [
    {
      "target_id": "string (UniProt)",
      "target_name": "string",
      "activity_type": "IC50|Ki|EC50|etc.",
      "activity_value": "float",
      "activity_unit": "nM|uM|etc.",
      "source": "string (database)",
      "reference": "string (PMID or DOI)"
    }
  ],
  "predicted_targets": [
    {
      "target_id": "string",
      "target_name": "string",
      "prediction_method": "SwissTarget|SEA|PharmMapper",
      "probability": "float",
      "validation_status": "predicted|validated"
    }
  ],
  "admet": {
    "absorption": {},
    "distribution": {},
    "metabolism": {},
    "excretion": {},
    "toxicity": {}
  }
}
```

### 8.4 API Access Summary

| Database | API Type | Authentication | Batch Support |
|----------|----------|----------------|---------------|
| COCONUT | REST (OpenAPI) | None | Yes |
| LOTUS/Wikidata | SPARQL | None | Yes |
| ChEMBL | REST | None | Yes |
| NPAtlas | REST | None | Yes |
| BindingDB | REST | None | Yes |
| PubChem | REST (PUG) | None | Yes |
| ChEBI | Web services + libChEBI | None | Yes |
| DGIdb | GraphQL | None | Yes |
| STITCH | REST | None | Limited |
| SwissTargetPrediction | Web only | None | No |
| PharmMapper | Web only | None | No |

### 8.5 Licensing Summary

| License | Databases | Commercial Use |
|---------|-----------|----------------|
| **CC0** | COCONUT, LOTUS | Yes, unrestricted |
| **CC BY 4.0** | NPAtlas, ChEBI, STITCH, STRING | Yes, with attribution |
| **CC BY-SA 3.0** | ChEMBL | Yes, share-alike |
| **Academic Only** | TCMSP, NPASS, IMPPAT, SuperNatural | Check terms |
| **Commercial** | DNP | Subscription required |

---

## References

### Primary Publications

1. LOTUS: Rutz A, et al. (2022) "The LOTUS initiative for open knowledge management in natural products research." eLife 11:e70780.

2. COCONUT 2.0: Venkata C, et al. (2024) "COCONUT 2.0: a comprehensive overhaul and curation of the collection of open natural products database." Nucleic Acids Res. gkae1063.

3. NPASS: Zeng X, et al. (2023) "NPASS database update 2023: quantitative natural product activity and species source database for biomedical research." Nucleic Acids Res. 51(D1):D621-D628.

4. NPAtlas 3.0: van Santen JA, et al. (2024) "Natural Products Atlas 3.0: extending the database of microbially derived natural products." Nucleic Acids Res. 53(D1):D691.

5. SuperNatural 3.0: Gallo K, et al. (2023) "SuperNatural 3.0-a database of natural products and natural product-based derivatives." Nucleic Acids Res. 51(D1):D654-D659.

6. TCMSP: Ru J, et al. (2014) "TCMSP: a database of systems pharmacology for drug discovery from herbal medicines." J Cheminform. 6:13.

7. HERB 2.0: Fang S, et al. (2024) "HERB 2.0: an updated database integrating clinical and experimental evidence for traditional Chinese medicine." Nucleic Acids Res. 53(D1):D1404.

8. IMPPAT 2.0: Vivek-Ananth RP, et al. (2023) "IMPPAT 2.0: An Enhanced and Expanded Phytochemical Atlas of Indian Medicinal Plants." ACS Omega 8(9):8827-8845.

9. CMAUP 2024: Zeng X, et al. (2024) "CMAUP database update 2024: extended functional and association information of useful plants for biomedical research." Nucleic Acids Res. gkad921.

10. ChEMBL 2023: Zdrazil B, et al. (2024) "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods." Nucleic Acids Res. 52(D1):D1180-D1192.

### Review Articles

- Sorokina M, Steinbeck C (2020) "Review on natural products databases: where to find data in 2020." J Cheminform. 12:20.
- Zhang R, et al. (2024) "Network pharmacology: a crucial approach in traditional Chinese medicine research." Chin Med. 19:106.

---

*Document generated: January 2026*
*Last updated: Based on database versions available as of late 2024/early 2025*
