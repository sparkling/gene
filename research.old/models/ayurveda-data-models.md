# Ayurveda and Indian Traditional Medicine Database Data Models

This document provides detailed data models, schemas, and field descriptions for major Ayurveda and Indian traditional medicine databases.

---

## Table of Contents

1. [IMPPAT 2.0](#1-imppat-20)
2. [NPACT](#2-npact)
3. [GRAYU](#3-grayu)
4. [OSADHI](#4-osadhi)
5. [TKDL](#5-tkdl)
6. [Cross-Reference Standards](#6-cross-reference-standards)

---

## 1. IMPPAT 2.0

**URL**: https://cb.imsc.res.in/imppat/
**Full Name**: Indian Medicinal Plants, Phytochemistry And Therapeutics
**Version**: 2.0 (Released June 17, 2022)
**Reference**: [Vivek-Ananth et al., ACS Omega 2023](https://pubs.acs.org/doi/10.1021/acsomega.3c00156)

### 1.1 Database Overview

| Metric | Count |
|--------|-------|
| Indian Medicinal Plants | 4,010 |
| Phytochemicals | 17,967 |
| Therapeutic Uses | 1,095 |
| Plant-Part-Phytochemical Associations | 189,386 |
| Plant-Part-Therapeutic Use Associations | 89,733 |
| Predicted Phytochemical-Target Interactions | 27,365 |

### 1.2 Core Entity Models

#### 1.2.1 Plant Entity

| Field | Description | Data Type |
|-------|-------------|-----------|
| Scientific Name | Botanical name of the plant | String |
| Common Name | Vernacular/common name | String |
| Plant Part | Specific part used (root, leaf, bark, etc.) | String |
| IMPPAT Plant ID | Internal plant identifier | String |

#### 1.2.2 Phytochemical Entity (6-Tab Structure)

##### Tab 1: Summary

| Field | Description | Data Type |
|-------|-------------|-----------|
| IMPPAT Phytochemical ID | Unique identifier (format: IMPHY######) | String |
| Phytochemical Name | Common chemical name | String |
| Chemical Super Class | High-level chemical classification | String |
| 2D Structure | Molecular structure image | Image/SVG |
| 3D Structure | Three-dimensional structure | PDB/MOL2 |
| Molecular Scaffold | Core structural framework | String |
| SMILES | Simplified molecular-input line-entry system | String |
| InChI | International Chemical Identifier | String |
| InChIKey | Hashed InChI identifier (27 characters) | String |

**Structure File Formats Available**: SDF, MOL, MOL2, PDB, PDBQT

##### Tab 2: Physicochemical Properties

| Field | Description | Data Type | Units |
|-------|-------------|-----------|-------|
| Molecular Weight | Molecular mass | Float | g/mol |
| log P | Lipophilicity (partition coefficient) | Float | - |
| TPSA | Topological Polar Surface Area | Float | Angstrom^2 |
| HBD | Hydrogen Bond Donors | Integer | count |
| HBA | Hydrogen Bond Acceptors | Integer | count |
| Rotatable Bonds | Number of rotatable bonds | Integer | count |
| Heavy Atoms | Non-hydrogen atom count | Integer | count |
| Heteroatoms | Non-carbon heavy atoms | Integer | count |
| Ring Count | Number of rings | Integer | count |
| Aromatic Ring Count | Number of aromatic rings | Integer | count |
| Fraction sp3 | Fraction of sp3 carbons | Float | 0-1 |
| Stereocenters | Number of stereogenic centers | Integer | count |

##### Tab 3: Drug-Likeness

| Rule/Score | Description | Pass/Fail Criteria | Data Type |
|------------|-------------|-------------------|-----------|
| Lipinski's Rule of Five (RO5) | Oral bioavailability predictor | MW<=500, logP<=5, HBD<=5, HBA<=10 | Boolean + Violation Count |
| Ghose Rule | Drug-like range filter | MW:160-480, logP:-0.4-5.6, atoms:20-70, MR:40-130 | Boolean + Violation Count |
| Veber Rule | Oral bioavailability | Rotatable bonds<=10, TPSA<=140 | Status (Good/Poor) |
| Egan Rule | Passive absorption | logP<=5.88, TPSA<=131.6 | Status (Good/Poor) |
| Pfizer 3/75 Rule | Toxicity risk | logP>3 AND TPSA<75 flagged | Status (Good/Poor) |
| GSK 4/400 Rule | Oral drug-likeness | MW<=400, logP<=4 | Status (Good/Poor) |
| QEDw Score | Quantitative Estimate of Drug-likeness (weighted) | 0-1 scale | Float |

**Drug-Like Subset**: 1,335 phytochemicals pass all six rules

##### Tab 4: ADMET Properties (via SwissADME)

**Absorption:**

| Field | Description | Data Type |
|-------|-------------|-----------|
| GI Absorption | Gastrointestinal absorption prediction | Categorical (High/Low) |
| Skin Permeability (log Kp) | Transdermal permeation | Float (cm/s) |
| Water Solubility Class | ESOL model prediction | Categorical (Insoluble to Highly Soluble) |

**Distribution:**

| Field | Description | Data Type |
|-------|-------------|-----------|
| BBB Permeant | Blood-Brain Barrier permeation | Boolean (Yes/No) |
| P-gp Substrate | P-glycoprotein efflux substrate | Boolean (Yes/No) |

**Metabolism:**

| Field | Description | Data Type |
|-------|-------------|-----------|
| CYP1A2 Inhibitor | Cytochrome P450 1A2 inhibition | Boolean (Yes/No) |
| CYP2C19 Inhibitor | Cytochrome P450 2C19 inhibition | Boolean (Yes/No) |
| CYP2C9 Inhibitor | Cytochrome P450 2C9 inhibition | Boolean (Yes/No) |
| CYP2D6 Inhibitor | Cytochrome P450 2D6 inhibition | Boolean (Yes/No) |
| CYP3A4 Inhibitor | Cytochrome P450 3A4 inhibition | Boolean (Yes/No) |

**Toxicity Alerts:**

| Field | Description | Data Type |
|-------|-------------|-----------|
| PAINS Alerts | Pan-Assay Interference compounds | Integer (count) |
| Brenk Alerts | Structural toxicity alerts (105 patterns) | Integer (count) |

**Note**: 493 phytochemicals lack ADMET predictions due to SMILES length restrictions in SwissADME.

##### Tab 5: Descriptors

| Category | Description | Count |
|----------|-------------|-------|
| 2D Descriptors | Constitutional, topological, connectivity indices | 1,444 |
| 3D Descriptors | Geometrical, surface area, volume descriptors | 431 |
| **Total** | Computed using PaDEL-Descriptor | **1,875** |

**Major Descriptor Classes (PaDEL):**
- Constitutional descriptors
- Topological descriptors
- Walk and path counts
- Connectivity indices
- Atom type electrotopological state
- BCUT descriptors
- Burden eigenvalues
- McGowan volume
- Molecular linear free energy relation (MLFER)
- Ring counts
- Laggner substructure counts
- Fingerprints (10 types, 16,092 bits total)

##### Tab 6: Predicted Human Target Proteins

| Field | Description | Data Type |
|-------|-------------|-----------|
| Target Protein | Human protein target | String |
| HGNC Symbol | HUGO Gene Nomenclature Committee symbol | String |
| Confidence Score | Interaction prediction confidence | Integer (0-1000) |
| Evidence Channel | Source of prediction | Categorical |

**Confidence Score Thresholds:**
- Minimum threshold: 150 (low confidence)
- Medium confidence: >= 400
- High confidence: >= 700 (used in IMPPAT)
- Highest confidence: >= 900

**Evidence Channels:**
- Experiments (ChEMBL, PDSP Ki, PDB)
- Databases (DrugBank, KEGG, CTD)
- Text Mining (PubMed, MEDLINE)
- Predicted (structure-based)

### 1.3 Association Models

#### Plant-Phytochemical Association

| Field | Description |
|-------|-------------|
| Plant Scientific Name | Source plant |
| Plant Part | Specific plant tissue |
| IMPPAT Phytochemical ID | Associated compound |
| Phytochemical Name | Chemical name |
| Reference | Literature citation |

#### Plant-Therapeutic Use Association

| Field | Description |
|-------|-------------|
| Plant Scientific Name | Source plant |
| Plant Part | Specific plant tissue |
| Therapeutic Use | Medical application |
| Therapeutic Use ID | Standardized identifier |
| Reference | Literature citation |

### 1.4 ID System

| ID Type | Format | Example |
|---------|--------|---------|
| Phytochemical ID | IMPHY + 6 digits | IMPHY011737 |
| Therapeutic Use ID | Internal numbering | TU###### |

---

## 2. NPACT

**URL**: https://webs.iiitd.edu.in/raghava/npact/
**Full Name**: Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target Database
**Reference**: [Mangal et al., Nucleic Acids Research 2013](https://academic.oup.com/nar/article/41/D1/D1124/1052661)

### 2.1 Database Overview

| Metric | Count |
|--------|-------|
| Bioactive Compounds | 1,574 |
| Cancer Cell Lines | 353 |
| In-vitro Compound-Cell Line Interactions | 5,214 |
| Compound-Target Interactions | ~1,980 |
| Source Articles | 762 |

### 2.2 Eight-Table Schema

#### Table 1: General Information

| Field | Description | Data Type |
|-------|-------------|-----------|
| NPACT ID | Unique identifier (format: NPACT#####) | String |
| Compound Name | Chemical name | String |
| IUPAC Name | Systematic chemical name | String |
| Synonyms | Alternative names | String (comma-separated) |
| Class | Chemical class (e.g., Terpenoid, Flavonoid) | String |
| PubChem ID | PubChem compound identifier | Integer |
| CAS Number | Chemical Abstracts Service registry number | String |
| InChI | International Chemical Identifier | String |
| InChIKey | Hashed InChI | String (27 chars) |
| SMILES | Canonical SMILES notation | String |
| SMART | SMARTS pattern | String |

#### Table 2: In Vitro Activity

| Field | Description | Data Type | Units |
|-------|-------------|-----------|-------|
| Cancer Type | Type of cancer | String | - |
| Cell Line | Cancer cell line name | String | - |
| IC50 | Half maximal inhibitory concentration | Float | uM or nM |
| ED50 | Median effective dose | Float | uM or nM |
| EC50 | Half maximal effective concentration | Float | uM or nM |
| GI50 | Growth inhibition 50% | Float | uM or nM |
| Activity Type | Type of measurement used | String | - |
| Reference PMID | PubMed article ID | Integer | - |

#### Table 3: In Vivo Activity

| Field | Description | Data Type |
|-------|-------------|-----------|
| Model System | Animal model used | String |
| Protein Target | Target protein name | String |
| Observation | Experimental findings | Text |
| Reference PMID | PubMed article ID | Integer |

#### Table 4: Target Information

| Field | Description | Data Type |
|-------|-------------|-----------|
| Target Name | Protein target name | String |
| Target Type | Class of target (kinase, protease, etc.) | String |
| UniProt ID | UniProt accession | String |
| Inhibition Type | Mode of inhibition | String |

#### Table 5: Cross-Reference Links

| Database | Field | Data Type |
|----------|-------|-----------|
| PubChem | CID | Integer |
| SuperNatural Database | SN ID | String |
| HIT (Herbal Ingredients' Targets) | HIT ID | String |
| CTD (Comparative Toxicogenomics) | CTD ID | String |
| NCI-60 GI50 Data | NCI ID | String |

#### Table 6: Property Table

| Property Category | Fields |
|-------------------|--------|
| Physical | Molecular weight, melting point, boiling point |
| Elemental | Atom counts (C, H, O, N, S, etc.) |
| Topological | Wiener index, Balaban index, connectivity indices |

#### Table 7: Drug-Likeness Filters

| Filter | Description | Pass Criteria |
|--------|-------------|---------------|
| Lipinski's Rule of 5 | Oral bioavailability | <=1 violation |
| Muegge's Filter | Drug-like properties | All criteria met |
| Ghose Filter | Qualified range | Within limits |
| Veber Filter | Oral bioavailability | All criteria met |

#### Table 8: Vendor/Supplier Table

| Field | Description |
|-------|-------------|
| Supplier Name | Commercial vendor name |
| Catalog Number | Vendor catalog ID |
| Availability | Stock status |

### 2.3 Chemical Class Distribution

| Class | Compound Count |
|-------|----------------|
| Terpenoids | 513 |
| Flavonoids | 329 |
| Alkaloids | 110 |
| Polyketides | 92 |
| Steroids | Variable |
| Lignans | Variable |
| Saponins | Variable |
| Others | Variable |

### 2.4 ID System

| ID Type | Format | Example |
|---------|--------|---------|
| Compound ID | NPACT + 5 digits | NPACT00001 |

---

## 3. GRAYU

**URL**: https://caps.ncbs.res.in/GRAYU/
**Full Name**: Graph-based Database integrating Ayurvedic formulations, medicinal plants, phytochemicals and diseases
**Reference**: [bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.10.26.684703v1.full)

### 3.1 Database Overview

| Metric | Count |
|--------|-------|
| Total Nodes | 157,010 |
| Total Relationships | 1,520,687 |
| Formulations | 1,039 |
| Plants | 12,949 |
| Phytochemicals | 129,542 |
| Diseases | 13,480 |

### 3.2 Graph Database Schema (Neo4j)

#### Node Types

##### Formulation Node

| Property | Description | Data Type |
|----------|-------------|-----------|
| name | Formulation name (Sanskrit/English) | String |
| formulation_id | Internal identifier | String |
| plant_ingredient_count | Number of plant ingredients | Integer |
| disease_association_count | Number of associated diseases | Integer |

##### Plant Node

| Property | Description | Data Type |
|----------|-------------|-----------|
| scientific_name | Botanical name | String |
| plant_id | Internal identifier | String |
| common_names | Vernacular names | Array[String] |
| family | Botanical family | String |

##### Phytochemical Node

| Property | Description | Data Type |
|----------|-------------|-----------|
| pubchem_cid | PubChem Compound ID | Integer |
| name | Chemical name | String |
| molecular_formula | Chemical formula | String |
| molecular_weight | Mass | Float |
| xlogp | Lipophilicity | Float |
| tpsa | Topological polar surface area | Float |
| hbd | Hydrogen bond donors | Integer |
| hba | Hydrogen bond acceptors | Integer |
| rotatable_bonds | Rotatable bond count | Integer |
| heavy_atom_count | Non-hydrogen atoms | Integer |
| exact_mass | Exact molecular mass | Float |
| monoisotopic_mass | Monoisotopic mass | Float |
| complexity | Molecular complexity score | Float |
| inchi | InChI notation | String |
| inchikey | Hashed InChI | String |
| smiles | SMILES notation | String |
| formal_charge | Net charge | Integer |
| defined_stereocenters | Defined stereocenter count | Integer |
| total_stereocenters | Total stereocenter count | Integer |
| covalent_unit_count | Covalent units | Integer |
| isotope_atom_count | Isotope atoms | Integer |
| brenk_violation | Brenk alert flag | Boolean |
| nih_violation | NIH alert flag | Boolean |
| pains_violation | PAINS alert flag | Boolean |

##### Disease Node

| Property | Description | Data Type |
|----------|-------------|-----------|
| disease_id | Internal identifier | String |
| name | Disease name | String |
| mesh_id | MeSH descriptor ID | String |
| doid | Disease Ontology ID | String |
| ayurvedic_term | Traditional Sanskrit term | String |

### 3.3 Relationship Types

| Relationship | Source | Target | Properties |
|--------------|--------|--------|------------|
| FOUND_IN | Phytochemical | Plant | source_db (String) |
| ASSOCIATED_WITH_DISEASE | Plant | Disease | evidence_type, source |
| IS_INGREDIENT_IN | Plant | Formulation | parts_used, quantity |
| ASSOCIATED_WITH | Formulation | Disease | ayurvedic_term, indication |
| IS_SUBCLASS_OF | Disease | Disease | (none - hierarchy from MeSH/DOID) |
| HAS_SYMPTOM | Disease | Disease | (from DOID) |

### 3.4 Disease Ontology Mapping

GRAYU maps Ayurvedic nosology (Sanskrit disease terms) to standardized ontologies:

| Ontology | Source | Coverage |
|----------|--------|----------|
| MeSH | Medical Subject Headings | Tree ID: C (Diseases) |
| DOID | Disease Ontology | Via OLS4 API |

**Mapping Method**: Manual semantic similarity matching between Ayurvedic terminology and standardized ontology labels.

### 3.5 Relationship Statistics

| Relationship Type | Count |
|-------------------|-------|
| Plant-Phytochemical | 1,382,362 |
| Plant-Disease | 116,824 |
| Plant-Formulation | 2,405 |
| Formulation-Disease | 4,087 |

### 3.6 Implementation Details

| Component | Technology |
|-----------|------------|
| Database | Neo4j (native graph database) |
| Backend | Flask (Python) |
| Frontend | HTML/CSS/JavaScript |
| Visualization | Cytoscape.js |
| Data Loading | APOC `apoc.periodic.iterate` |

---

## 4. OSADHI

**URL**: https://neist.res.in/osadhi/
**Full Name**: Online Structural and Analytics based Database for Herbs of India
**Institution**: CSIR-NEIST, Jorhat, Assam
**Reference**: [Computational Biology and Chemistry 2023](https://pubmed.ncbi.nlm.nih.gov/36512929/)

### 4.1 Database Overview

| Metric | Count |
|--------|-------|
| Medicinal Plants | 6,959 |
| Plant Families | 348 |
| Unique Phytochemicals | 27,440 (22,314 with full data) |
| Therapeutic Uses | 2,477 |
| States/UTs Covered | 28 states + 8 union territories |

### 4.2 Four-Feature Schema

#### Feature 1: Traditional Knowledge

| Field | Description | Data Type |
|-------|-------------|-----------|
| Scientific Name | Botanical name | String |
| Vernacular Names | Regional language names | Array[String] |
| Family | Botanical family | String |
| Plant Parts Used | Parts with medicinal value | Array[String] |
| Therapeutic Use | Medical applications | Array[String] |
| Taxonomy | Kingdom > Phylum > Class > Order > Family > Genus > Species | Hierarchical |

#### Feature 2: Geographical Classification

| Field | Description | Data Type |
|-------|-------------|-----------|
| State | Indian state(s) where found | Array[String] |
| Union Territory | Indian UT(s) where found | Array[String] |
| Geographic Distribution | Visual map representation | GeoJSON/Image |

**Coverage**: All 28 Indian states and 8 union territories

#### Feature 3: Phytochemicals

| Field | Description | Data Type |
|-------|-------------|-----------|
| Phytochemical Name | Chemical name | String |
| IUPAC Name | Systematic name | String |
| SMILES | Canonical SMILES | String |
| InChIKey | Hashed InChI | String |
| 2D Structure | Molecular structure | Image/SDF |
| 3D Structure | Three-dimensional structure | MOL2/PDB |
| NPClassifier Class | Deep learning classification | String |
| NPClassifier Superclass | Higher-level classification | String |
| NPClassifier Pathway | Biosynthetic pathway | String |

**Classification System**: NPClassifier (deep learning framework for natural product classification)

#### Feature 4: Chemoinformatics Analysis

##### Physicochemical Properties

| Property | Description | Data Type |
|----------|-------------|-----------|
| Molecular Weight | Mass | Float (Da) |
| LogP | Lipophilicity | Float |
| TPSA | Polar surface area | Float (Angstrom^2) |
| HBD | H-bond donors | Integer |
| HBA | H-bond acceptors | Integer |
| Rotatable Bonds | Flexible bonds | Integer |

##### ADMET Properties (computed via web servers)

| Category | Fields |
|----------|--------|
| Absorption | GI absorption, skin permeability |
| Distribution | BBB permeability, P-gp substrate |
| Metabolism | CYP inhibition profiles |
| Excretion | Renal clearance predictions |
| Toxicity | Hepatotoxicity, carcinogenicity alerts |

##### Antiviral Potency Predictions

| Model | Description |
|-------|-------------|
| Random Forest | ML classification for antiviral activity |
| XGBoost | Gradient boosting for antiviral prediction |

### 4.3 Machine Learning Classifications

| Application | Model | Output |
|-------------|-------|--------|
| Therapeutic Use Classification | Indigenous ML models | Predicted therapeutic category |
| Phytochemical Classification | NPClassifier | Class, Superclass, Pathway |
| Antiviral Potency | Random Forest, XGBoost | Potency score |

---

## 5. TKDL

**URL**: https://www.tkdl.res.in/
**Full Name**: Traditional Knowledge Digital Library
**Institution**: CSIR + Ministry of Health and Family Welfare, India
**Established**: 2001

### 5.1 Database Overview

| Metric | Count |
|--------|-------|
| Formulations/Practices | ~360,000 (3.6 lakh) |
| TKRC Subgroups | ~27,000 |
| Languages Available | 5 (English, German, French, Japanese, Spanish) |
| Systems Covered | 4 (Ayurveda, Unani, Siddha, Yoga) |

### 5.2 Formulation Record Structure

| Field | Description | Data Type |
|-------|-------------|-----------|
| TKDL ID | Unique formulation identifier | String |
| Title | Formulation name | String |
| Knowledge Resource | Source text/reference | String |
| Date Since Known | Historical dating | Date/Period |
| Country | Country of origin | String (India) |
| Contact Organization | Custodian organization | String |
| Abstract on Usage | Description of use | Text |
| Keywords | Searchable terms | Array[String] |
| IPC Code | International Patent Classification | String |
| TKRC Code | Traditional Knowledge Resource Classification | String |
| System of Medicine | Ayurveda/Unani/Siddha/Yoga | String |
| Ingredients | List of components | Array[Object] |
| Preparation Method Code | Standardized prep method (80+ codes) | String |
| Bibliography Code | Reference citation | String |

### 5.3 TKRC Classification System

The Traditional Knowledge Resource Classification (TKRC) provides hierarchical classification compatible with IPC.

#### Main Sections

| Section | System | Description |
|---------|--------|-------------|
| A | Ayurveda | Traditional Indian medicine |
| B | Unani | Greco-Arabic medicine tradition |
| C | Siddha | Tamil traditional medicine |
| Y | Yoga | Mind-body practices |

#### Classes (Under Each Section)

| Class | Description |
|-------|-------------|
| 01 | Pharmaceutical preparations (Kalpana) |
| 02 | Personal Hygiene Preparations |
| 03 | Dietary (Food/Foodstuff/Beverages) |
| 04 | Biocides, Fumigatives (Dhoopana, Krimighna) |

#### Subclasses (Pharmaceutical Preparations - 01)

| Subclass | Description |
|----------|-------------|
| 01A | Based on Plants (Audbhida) |
| 01B | Based on Animals (Jangama) |
| 01C | Based on Minerals (Parthiva) |
| 01D | Characterised by Diseases (Roga) |
| 01E | Characterised by Actions (Karma) |
| 01F | Mode of Administration |
| 01G | Miscellaneous |

#### Example TKRC Code Structure

```
A01A-1/1326
|  |  | |
|  |  | +-- Subgroup (specific plant: Nigella sativa)
|  |  +---- Group (1/00 = Whole Medicinal Plants)
|  +------- Subclass (A = Plant-based)
+---------- Section-Class (A01 = Ayurveda Pharmaceutical)
```

**IPC Mapping**: A01A-1/1326 maps to IPC A61K 36/71

### 5.4 Preparation Method Codes

Approximately 80 standardized codes for preparation techniques:

| Code | Method | Description |
|------|--------|-------------|
| AM1 | Arista | Fermentation method |
| AM2 | Asava | Self-generated fermentation |
| AM3 | Churna | Powder preparation |
| AM4 | Kwatha | Decoction |
| AM5 | Taila | Oil-based preparation |
| ... | ... | ... |

### 5.5 IPC Integration

TKRC enabled expansion of IPC classification:

| Before TKRC | After TKRC |
|-------------|------------|
| Few subgroups under A61K 35/00 | ~200 subgroups under A61K 36/00 |

### 5.6 Search Features

| Search Type | Description |
|-------------|-------------|
| Single Word | Basic keyword search |
| Multiple Word | AND/OR combinations |
| Boolean Expression | Complex logical queries |
| Proximity Search | Terms within n words |
| Field Search | Specific field targeting |
| Phrase Search | Exact phrase matching |
| IPC Search | By patent classification |
| TKRC Search | By traditional classification |

---

## 6. Cross-Reference Standards

### 6.1 Chemical Identifiers

| Identifier | Standard | Usage |
|------------|----------|-------|
| PubChem CID | NCBI | Primary compound ID |
| CAS Number | Chemical Abstracts Service | Legacy identifier |
| InChI | IUPAC | Structure encoding |
| InChIKey | IUPAC | 27-character hash |
| SMILES | Daylight | Linear notation |
| ChEMBL ID | EMBL-EBI | Bioactivity database |

### 6.2 Protein/Gene Identifiers

| Identifier | Database | Usage |
|------------|----------|-------|
| HGNC Symbol | HUGO Gene Nomenclature | Gene symbols |
| UniProt Accession | UniProt | Protein sequences |
| Ensembl ID | Ensembl | Gene/transcript IDs |
| PDB ID | RCSB PDB | 3D structures |

### 6.3 Disease/Ontology Standards

| Ontology | Source | Coverage |
|----------|--------|----------|
| MeSH | NLM | Medical Subject Headings |
| DOID | Disease Ontology | Disease classifications |
| ICD-10 | WHO | International disease classification |
| SNOMED CT | SNOMED International | Clinical terminology |

### 6.4 Patent Classification

| System | Organization | Usage |
|--------|--------------|-------|
| IPC | WIPO | International Patent Classification |
| TKRC | India CSIR | Traditional Knowledge Classification |
| CPC | EPO/USPTO | Cooperative Patent Classification |

---

## 7. Data Access Summary

| Database | Format | API | Download |
|----------|--------|-----|----------|
| IMPPAT 2.0 | Web interface | No | Structure files (SDF, MOL, PDB) |
| NPACT | Web interface | No | MOL format structures |
| GRAYU | Neo4j graph | No | Export via interface |
| OSADHI | Web interface | No | Limited export |
| TKDL | Restricted access | No | Not publicly downloadable |

---

## 8. References

1. Vivek-Ananth RP, et al. "IMPPAT 2.0: An Enhanced and Expanded Phytochemical Atlas of Indian Medicinal Plants." ACS Omega. 2023;8(9):8827-8845. https://pubs.acs.org/doi/10.1021/acsomega.3c00156

2. Mangal M, et al. "NPACT: Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target database." Nucleic Acids Research. 2013;41(D1):D1124-D1129. https://academic.oup.com/nar/article/41/D1/D1124/1052661

3. GRAYU: Graph-based Database integrating Ayurvedic formulations, medicinal plants, phytochemicals and diseases. bioRxiv 2025. https://www.biorxiv.org/content/10.1101/2025.10.26.684703v1.full

4. OSADHI - An online structural and analytics based database for herbs of India. Computational Biology and Chemistry. 2023. https://pubmed.ncbi.nlm.nih.gov/36512929/

5. Traditional Knowledge Digital Library (TKDL). https://www.tkdl.res.in/

6. SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. Scientific Reports. 2017;7:42717. https://pmc.ncbi.nlm.nih.gov/articles/PMC5335600/

7. Yap CW. "PaDEL-descriptor: An open source software to calculate molecular descriptors and fingerprints." J Comput Chem. 2011;32(7):1466-74. https://pubmed.ncbi.nlm.nih.gov/21425294/

---

*Document generated: January 2026*
*Last updated: January 18, 2026*
