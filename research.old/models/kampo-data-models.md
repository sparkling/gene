# Kampo (Japanese Traditional Medicine) Database Data Models

This document provides comprehensive data models, schemas, and field descriptions for major Kampo medicine databases used in traditional Japanese medicine research.

---

## Table of Contents

1. [KampoDB](#1-kampodb)
2. [STORK](#2-stork)
3. [TradMPD](#3-tradmpd)
4. [KNApSAcK KAMPO](#4-knapsack-kampo)
5. [EKAT](#5-ekat-evidence-reports-of-kampo-treatment)
6. [Cross-Reference Systems](#6-cross-reference-systems)

---

## 1. KampoDB

**URL**: https://wakanmoview.inm.u-toyama.ac.jp/kampo/
**Institution**: Institute of Natural Medicine, University of Toyama
**License**: CC BY-SA 4.0

### 1.1 Overview

KampoDB is an integrated platform for mode-of-action analysis and repositioning of natural medicines. It provides comprehensive scientific resources on Kampo formulas, crude drugs, constituent compounds, and target proteins.

**Database Statistics (Current Version)**:
- Kampo Formulas: 298 entries
- Crude Drugs: 180 entries
- Natural Compounds: 3,002 entries
- Proteins/Genes: 62,906 entries
- Docking Simulation Results: 3,063,505 records

### 1.2 Four-Layer Hierarchical Structure

```
Layer 1: Kampo Medicine (e.g., "Kakkonto")
    |
    +-- Layer 2: Crude Drugs (e.g., "Ephedra Herb")
            |
            +-- Layer 3: Compounds (e.g., "Methylephedrine")
                    |
                    +-- Layer 4: Target Proteins (e.g., "ADRA1D")
```

### 1.3 Entity Schemas

#### 1.3.1 Formula Entity

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `id` | String | Unique alphanumeric identifier | "KT" |
| `name` | String | Romanized formula name | "Kakkonto" |
| `name_jp` | String (Unicode) | Japanese kanji name | "葛根湯" |
| `kana` | String (Unicode) | Japanese hiragana reading | "かっこんとう" |
| `synonyms` | String | Alternative English names | "" |
| `synonyms_jp` | String (Unicode) | Alternative Japanese names | "" |
| `link_dento` | String | TradMPD cross-reference URL | "kakkonto/" |

**API Endpoint**: `/api/formula/`

**Example Record**:
```json
{
  "id": "KT",
  "name": "Kakkonto",
  "name_jp": "\u845b\u6839\u6e6f",
  "kana": "\u304b\u3063\u3053\u3093\u3068\u3046",
  "synonyms": "",
  "synonyms_jp": "",
  "link_dento": "kakkonto/"
}
```

#### 1.3.2 Crude Drug Entity

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `id` | Integer | Unique numeric identifier | 1 |
| `name` | String | English common name | "Agarwood" |
| `name_jp` | String (Unicode) | Japanese name | "沈香" |
| `kana` | String (Unicode) | Japanese reading | "じんこう" |
| `orgn` | String | Botanical/zoological sources (comma-separated) | "Aquilaria agallocha, Aquilaria malaccensis" |
| `synonyms_jp` | String (Unicode) | Alternative Japanese names | "" |

**API Endpoint**: `/api/crude/`

**Example Record**:
```json
{
  "id": 1,
  "name": "Agarwood",
  "name_jp": "\u6c88\u9999",
  "kana": "\u3058\u3093\u3053\u3046",
  "orgn": "Aquilaria agallocha Roxburgh, Aquilaria malaccensis Lamarck, Aquilaria crassna Pierre",
  "synonyms_jp": ""
}
```

#### 1.3.3 Compound Entity

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `id` | Integer | PubChem CID or internal ID | 969516 |
| `name` | String | Chemical compound name | "Curcumin" |
| `knapsack_id` | String | KNApSAcK compound ID | "C00000123" |
| `pubchem_id` | Integer | PubChem compound ID | 969516 |
| `chembl_id` | String | ChEMBL compound ID | "CHEMBL116438" |

**API Endpoint**: `/api/compound/`

**Example Record**:
```json
{
  "id": 969516,
  "name": "Curcumin"
}
```

#### 1.3.4 Protein Entity

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `id` | Integer | KEGG GENES identifier | 1 |
| `name` | String | Gene symbol | "A1BG" |
| `aliases` | String | Alternative names (comma-separated) | "A1B, ABG, GAB, HYST2477" |
| `description` | String | Protein functional description | "alpha-1-B glycoprotein" |

**API Endpoint**: `/api/protein/`

**Example Record**:
```json
{
  "id": 1,
  "name": "A1BG",
  "aliases": "A1B, ABG, GAB, HYST2477",
  "description": "alpha-1-B glycoprotein"
}
```

### 1.4 Relationship Schemas

#### 1.4.1 Formula-Crude Drug Relationship

**API Endpoint**: `/api/formula/{id}/crude`

```json
[
  {"id": 40, "name": "Cinnamon Bark"},
  {"id": 62, "name": "Ephedra Herb"},
  {"id": 76, "name": "Ginger"},
  {"id": 78, "name": "Glycyrrhiza"},
  {"id": 90, "name": "Jujube"},
  {"id": 126, "name": "Peony Root"},
  {"id": 147, "name": "Pueraria Root"}
]
```

#### 1.4.2 Formula-Compound Relationship

**API Endpoint**: `/api/formula/{id}/compound`

Returns array of compounds with their PubChem CIDs and names.

### 1.5 Enrichment Analysis Schema

#### 1.5.1 Pathway Enrichment

**API Endpoint**: `/api/formula/{id}/pathway`

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `id` | String | KEGG pathway identifier | "hsa05417" |
| `name` | String | Pathway description | "Lipid and atherosclerosis" |
| `count` | Integer | Number of associated genes | 131 |
| `fdr` | Float | False discovery rate (adjusted p-value) | 2.2e-56 |

**Example Record**:
```json
{
  "id": "hsa05417",
  "name": "Lipid and atherosclerosis",
  "count": 131,
  "fdr": 2.2e-56
}
```

#### 1.5.2 Disease Enrichment

**API Endpoint**: `/api/formula/{id}/disease`

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `id` | String | Disease identifier | "disease_001" |
| `name` | String | Disease name | "Hepatocellular carcinoma" |
| `count` | Integer | Number of associated genes | 15 |
| `fdr` | Float | False discovery rate | 0.00027 |

### 1.6 Target Prediction Schema

#### 1.6.1 Ligand-Based Prediction (LBP)

**API Endpoint**: `/api/compound/{id}/lbp`

Uses TESS algorithm with KCF-S descriptors (475,692 features) and Jaccard correlation coefficient.

#### 1.6.2 Structure-Based Prediction (SBP)

**API Endpoint**: `/api/compound/{id}/sbp`

Docking simulation using AutoDock with protein structures from PDB and SAHG databases.

| Field | Data Type | Description |
|-------|-----------|-------------|
| `protein_id` | Integer | Target protein identifier |
| `binding_energy` | Float | Predicted binding free energy (kcal/mol) |
| `conformation` | Object | 3D binding pose coordinates |

#### 1.6.3 Docking Simulation

**API Endpoint**: `/api/docking/compound/{compound_id}/protein/{protein_id}`

Returns binding affinity predictions and optimal docking conformations.

### 1.7 Complete API Reference

| Endpoint Pattern | Description | Response |
|-----------------|-------------|----------|
| `/api/formula/` | List all formulas | Array of {id, name} |
| `/api/crude/` | List all crude drugs | Array of {id, name} |
| `/api/compound/` | List all compounds | Array of {id, name} |
| `/api/protein/` | List all proteins | Array of {id, name} |
| `/api/{entity}/{id}/info` | Entity details | Full entity object |
| `/api/formula/{id}/crude` | Formula's crude drugs | Array of crude drugs |
| `/api/formula/{id}/compound` | Formula's compounds | Array of compounds |
| `/api/formula/{id}/protein` | Formula's target proteins | Array of proteins |
| `/api/{entity}/{id}/pathway` | Pathway enrichment | Array of pathways |
| `/api/{entity}/{id}/disease` | Disease enrichment | Array of diseases |
| `/api/{entity}/{id}/process` | GO biological process | Array of GO terms |
| `/api/{entity}/{id}/function` | GO molecular function | Array of GO terms |
| `/api/compound/{id}/lbp` | Ligand-based prediction | Predicted targets |
| `/api/compound/{id}/sbp` | Structure-based prediction | Docking results (max 1000) |
| `/api/docking/compound/{cid}/protein/{pid}` | Specific docking | Binding data |

---

## 2. STORK

**URL**: http://mpdb.nibiohn.go.jp/stork/
**Full Name**: Standards of Reporting Kampo Products
**Institution**: National Institute of Health Sciences (NIHS) & NIBIOHN

### 2.1 Overview

STORK provides standardized reference information for 148 Kampo formulas approved in Japan, designed to ensure consistent reporting in clinical research articles.

### 2.2 Kampo Formula Schema

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `formula_id` | String | Unique formula identifier | "TJ-001" |
| `formula_name_jp` | String | Japanese name (kanji) | "葛根湯" |
| `formula_name_romaji` | String | Romanized name | "Kakkonto" |
| `formula_name_en` | String | English name | "Pueraria Decoction" |
| `manufacturer` | String | Pharmaceutical company | "Tsumura" |
| `approval_date` | Date | Ministry approval date | "1986-04-01" |

### 2.3 Crude Drug Composition Schema

| Field | Data Type | Description | Unit |
|-------|-----------|-------------|------|
| `crude_drug_name_jp` | String | Japanese name | - |
| `crude_drug_name_latin` | String | Latin pharmacopoeia name | - |
| `amount` | Float | Weight in formula | grams |
| `ratio` | Float | Percentage of total | % |
| `part_used` | String | Plant part specification | - |

**Example Composition Record**:
```json
{
  "formula_id": "TJ-001",
  "crude_drugs": [
    {
      "crude_drug_name_jp": "葛根",
      "crude_drug_name_latin": "Puerariae Radix",
      "amount": 4.0,
      "ratio": 25.0,
      "part_used": "root"
    },
    {
      "crude_drug_name_jp": "麻黄",
      "crude_drug_name_latin": "Ephedrae Herba",
      "amount": 3.0,
      "ratio": 18.75,
      "part_used": "aerial parts"
    }
  ]
}
```

### 2.4 Glycyrrhiza Dosing Data

Glycyrrhizae Radix (licorice) is the most frequently used crude drug, requiring special dosing attention due to pseudoaldosteronism risk.

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `formula_id` | String | Formula identifier | "TJ-068" |
| `glycyrrhiza_amount` | Float | Amount per daily dose (grams) | 4.0 |
| `glycyrrhizin_content` | Float | GL content (mg) | 120.5 |
| `risk_category` | String | Precaution level | "high" |

**Dosing Standards**:
- Initial: 0.5 g/day
- Increment: 0.25-0.5 g/day
- Maximum: 3.0 g/day

**Glycyrrhiza Content by Formula**:
| Formula | Glycyrrhizae Radix (g) | Risk Level |
|---------|------------------------|------------|
| Shakuyaku-kanzo-To | 4.0 | High |
| Sho-seiryu-To | 3.0 | Moderate |
| Hange-shashin-To | 2.5 | Moderate |

### 2.5 Product Information Schema

| Field | Data Type | Description |
|-------|-----------|-------------|
| `product_code` | String | Manufacturer product code |
| `dosage_form` | String | Extract granule, tablet, etc. |
| `daily_dose` | String | Recommended daily dosage |
| `package_insert_url` | URL | Link to official documentation |
| `nhi_price` | Float | National Health Insurance price (yen) |

---

## 3. TradMPD

**URL**: https://dentomed.toyama-wakan.net/
**Full Name**: Traditional Medical & Pharmaceutical Database
**Institution**: Institute of Natural Medicine, University of Toyama

### 3.1 Overview

TradMPD is the first repository of experimental data for natural medicines used in Kampo, including genetic data, LC-MS profiling, and biological activity data.

### 3.2 Compound Library Schema

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `compound_id` | String | Internal identifier | "SCH000679" |
| `compound_name` | String | Chemical name | "Capillene" |
| `molecular_formula` | String | Elemental composition | "C12H10" |
| `average_mass` | Float | Average molecular mass (Da) | 154.21 |
| `exact_mass` | Float | Monoisotopic mass (Da) | 154.07825 |
| `iupac_name` | String | Systematic chemical name | "hexa-2,4-diynylbenzene" |
| `inchi_key` | String | InChI hash | "WXQYRBLGGSLJHA-UHFFFAOYSA-N" |
| `pubchem_id` | Integer | PubChem CID | 3083613 |
| `cas_number` | String | CAS Registry Number | - |
| `nikkaji_id` | String | NIKKAJI database ID | - |

**Example Compound Record**:
```json
{
  "compound_id": "SCH000679",
  "compound_name": "Capillene",
  "molecular_formula": "C12H10",
  "average_mass": 154.21,
  "exact_mass": 154.07825,
  "iupac_name": "hexa-2,4-diynylbenzene",
  "inchi_key": "WXQYRBLGGSLJHA-UHFFFAOYSA-N",
  "pubchem_id": 3083613,
  "source_crude_drug": "Artemisia Capillaris Flower"
}
```

### 3.3 LC-MS Profiling Data Format

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `sample_id` | String | Sample identifier | "EH-001" |
| `crude_drug` | String | Source crude drug | "Ephedra Herb" |
| `ionization_mode` | String | MS ionization | "positive" / "negative" |
| `mh_positive` | Float | [M+H]+ m/z value | 155.086 |
| `mh_negative` | Float | [M-H]- m/z value | 153.07 |
| `retention_time` | Float | HPLC retention time (min) | 12.5 |
| `peak_area` | Float | Integrated peak area | 1234567 |

**Mass Spectrum Data Structure**:
```
# Format specification for LC-MS data files
Line 1: : Comment line (starts with colon)
Line 2: : Attribute headers (tab-separated)
Line 3+: m/z value [TAB] intensity_condition1 [TAB] intensity_condition2 ...
```

### 3.4 Crude Drug Schema

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `crude_drug_id` | String | Internal identifier | "SCC000076" |
| `market_name` | String (Unicode) | Japanese market name | "麻黄" |
| `latin_name` | String | Pharmacopoeia name | "Ephedrae Herba" |
| `tmpw_catalog` | String | TMPW catalog number | "18571" |
| `family` | String | Botanical family | "Ephedraceae" |
| `part_used` | String | Plant part description | "Dried green young branch" |
| `quality_criteria` | String | Quality specifications | "Freshness and green coloration" |

### 3.5 Biological Activity Data Structure

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `activity_id` | String | Activity record ID | "BA-001" |
| `crude_drug_id` | String | Source crude drug | "SCC000076" |
| `activity_category` | String | Activity classification | "Pharmacological" |
| `activity_type` | String | Specific activity | "Antitussive" |
| `mechanism` | String | Mode of action | "Sympathetic nerve excitation" |
| `test_system` | String | Assay method | "in vivo mouse model" |
| `reference` | String | Literature citation | "PMID:12345678" |

**Example Activity Record**:
```json
{
  "crude_drug_id": "SCC000076",
  "crude_drug_name": "Ephedra Herb",
  "activities": [
    {
      "category": "Pharmacological",
      "type": "Antitussive",
      "mechanism": "Beta-adrenergic stimulation"
    },
    {
      "category": "Pharmacological",
      "type": "Anti-inflammatory",
      "mechanism": "COX inhibition"
    }
  ],
  "clinical_indications": [
    "Diaphoretic",
    "Antifebrile",
    "Respiratory conditions"
  ]
}
```

### 3.6 Chemical Composition by Class

Compounds are organized by chemical class:

| Class | Description | Example Compounds |
|-------|-------------|-------------------|
| `alkaloids.isoquinoline` | Isoquinoline alkaloids | Ephedrine, Pseudoephedrine |
| `alkaloids.other` | Other alkaloids | Ephedroxane |
| `aliphatics` | Aliphatic compounds | Nonacosans, Tricosinol |
| `tannins` | Polyphenolic compounds | Catechin polymers |
| `flavonoids` | Flavonoid compounds | Quercetin glycosides |

### 3.7 Genetic Data Architecture

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `species` | String | Plant species | "Ephedra sinica" |
| `ipni_id` | String | IPNI identifier | "1234567-1" |
| `gene_region` | String | Sequenced region | "ITS", "rbcL", "matK" |
| `genome_type` | String | Nuclear/Chloroplast/Mitochondrial | "chloroplast" |
| `genbank_accession` | String | GenBank accession | "AB123456" |
| `sequence` | String | DNA sequence | "ATGCATGC..." |

---

## 4. KNApSAcK KAMPO

**URL**: http://www.knapsackfamily.com/kampo/top.php
**Institution**: NAIST (Nara Institute of Science and Technology)

### 4.1 Overview

KNApSAcK is a comprehensive species-metabolite relationship database. The KAMPO module focuses on traditional Japanese medicine formulations.

**Database Statistics**:
- Total Metabolites: 63,715 entries
- Species-Metabolite Pairs: 159,095 relationships
- KAMPO Formulas: 336 entries
- Medicinal Plants (KAMPO): 278 species

### 4.2 Compound ID System (C_ID)

The KNApSAcK compound identifier follows a standardized format:

| Component | Format | Description | Example |
|-----------|--------|-------------|---------|
| Prefix | "C" | Compound indicator | C |
| Number | 8 digits | Sequential identifier | 00001234 |
| Full ID | C + 8 digits | Complete identifier | C00001234 |

**ID Range**: C00000001 to C99999999

### 4.3 Metabolite Entity Schema

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `c_id` | String | KNApSAcK compound ID | "C00001234" |
| `metabolite_name` | String | Common name | "Quercetin" |
| `molecular_formula` | String | Chemical formula | "C15H10O7" |
| `molecular_weight` | Float | Molecular mass (Da) | 302.236 |
| `cas_number` | String | CAS Registry Number | "117-39-5" |
| `inchi` | String | InChI string | "InChI=1S/C15H10O7/..." |
| `smiles` | String | SMILES notation | "O=c1c(O)c(-c2ccc..." |

### 4.4 Species-Metabolite Relationship Schema

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `relationship_id` | Integer | Unique relationship ID | 12345 |
| `c_id` | String | Compound identifier | "C00001234" |
| `species_name` | String | Scientific name | "Pueraria lobata" |
| `family` | String | Botanical family | "Fabaceae" |
| `part` | String | Plant part | "root" |
| `reference` | String | Literature source | "PMID:12345678" |

### 4.5 KAMPO Formula Composition Schema

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `formula_id` | String | Formula identifier | "KP-001" |
| `formula_name_jp` | String | Japanese name | "葛根湯" |
| `formula_name_romaji` | String | Romanized name | "Kakkonto" |
| `medicinal_plants` | Array | List of constituent plants | ["Pueraria lobata", "Ephedra sinica"] |
| `total_crude_drugs` | Integer | Number of crude drugs | 7 |

### 4.6 Search Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `metabolite_name` | String | Search by compound name | "quercetin" |
| `organism` | String | Scientific name search | "Glycyrrhiza" |
| `molecular_weight` | Float | MW with margin | 302.0 |
| `margin` | Float | MW search tolerance | 1.0 |
| `molecular_formula` | String | Exact formula match | "C15H10O7" |

### 4.7 KNApSAcK Family Databases

| Database | Description | Size |
|----------|-------------|------|
| KNApSAcK Core | Species-metabolite relationships | 101,500 pairs |
| KNApSAcK KAMPO | Japanese herbal formulas | 336 formulas |
| KNApSAcK JAMU | Indonesian traditional medicine | 5,310 formulas |
| KNApSAcK WorldMap | Geographic distribution | 41,548 GZ-plant pairs |
| KNApSAcK 3D | 3D metabolite structures | 50,048 structures |
| KNApSAcK Activity | Metabolite-activity relationships | 9,584 triplets |

### 4.8 Cross-Reference IDs

| External Database | ID Field | Example |
|-------------------|----------|---------|
| PubChem | `pubchem_cid` | 5280343 |
| ChEBI | `chebi_id` | "CHEBI:16243" |
| KEGG | `kegg_id` | "C00389" |
| CAS | `cas_number` | "117-39-5" |
| ChemSpider | `chemspider_id` | 4444051 |

---

## 5. EKAT (Evidence Reports of Kampo Treatment)

**URL**: http://www.jsom.or.jp/medical/ebm/ere/index.html
**Institution**: Japan Society for Oriental Medicine (JSOM)

### 5.1 Overview

EKAT is a compilation of structured abstracts of randomized controlled trials (RCTs) for Kampo medicines, published since 2007.

**Database Statistics (EKAT 2016 Appendix 2018)**:
- Total RCTs: 488 studies
- Meta-analyses: 5 studies
- Total Publications: 595 papers

### 5.2 Twelve-Component Structured Abstract Format

The EKAT structured abstract consists of 12 mandatory components, based on the Altman-Gardner format (1987) with extensions:

| Component | Field Name | Description | Data Type |
|-----------|------------|-------------|-----------|
| 1 | `objective` | Research question/hypothesis | Text |
| 2 | `design` | Study design type | Enum |
| 3 | `setting` | Clinical setting/location | Text |
| 4 | `participants` | Subject characteristics | Text |
| 5 | `intervention` | Kampo treatment details | Text |
| 6 | `control` | Comparison group treatment | Text |
| 7 | `main_outcome` | Primary outcome measure | Text |
| 8 | `results` | Key findings with statistics | Text |
| 9 | `conclusion` | Authors' conclusions | Text |
| 10 | `kampo_discussion` | Kampo-specific interpretation | Text/Boolean |
| 11 | `abstractor_comments` | Third-party review notes | Text |
| 12 | `reference` | Full bibliographic citation | Text |

### 5.3 Study Design Classification

| Code | Design Type | Description |
|------|-------------|-------------|
| `RCT` | Randomized Controlled Trial | Standard parallel-group RCT |
| `DB-RCT` | Double-Blind RCT | Blinded to both participant and investigator |
| `quasi-RCT` | Quasi-Randomized Trial | Alternation or other non-random allocation |
| `RCT-envelope` | Envelope Randomization | Sealed envelope method |
| `RCT-crossover` | Crossover RCT | Within-subject crossover design |
| `meta-analysis` | Meta-Analysis | Systematic review with quantitative synthesis |

### 5.4 Data Source Classification

| Code | Source | Description |
|------|--------|-------------|
| `C` | CENTRAL | Cochrane Library (includes PubMed/MEDLINE & EMBASE) |
| `I` | Ichushi | Igaku Chuo Zasshi (Japanese medical literature) |
| `N` | JKMA | Japan Kampo Medicines Manufacturers' Association |

### 5.5 ICD-10 Disease Classification

Studies are organized by ICD-10 categories:

| Category | Description | RCT Count |
|----------|-------------|-----------|
| I | Infectious diseases | Variable |
| II | Neoplasms | Variable |
| IV | Endocrine/metabolic | Variable |
| V | Mental/behavioral | Variable |
| VI | Nervous system | Variable |
| IX | Circulatory system | Variable |
| X | Respiratory system | Variable |
| XI | Digestive system | Variable |
| XIII | Musculoskeletal | Variable |
| XIV | Genitourinary | Variable |

### 5.6 EKAT Record Schema

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `ekat_id` | String | Unique EKAT identifier | "EKAT-2016-001" |
| `icd10_code` | String | Disease classification | "K21" |
| `kampo_formula` | String | Formula studied | "Rikkunshito" |
| `study_design` | Enum | Design classification | "DB-RCT" |
| `data_source` | Enum | Database source | "C" |
| `language` | Enum | Original publication language | "ja" / "en" |
| `publication_year` | Integer | Year published | 2015 |
| `reference_pmid` | String | PubMed ID if available | "12345678" |
| `pdf_url` | URL | Link to full structured abstract | "/ekat/2016/001.pdf" |

### 5.7 Inclusion Criteria

Studies must meet all criteria:
1. Uses Kampo formulae approved for manufacture/sale in Japan
2. Study design: RCT, quasi-RCT, crossover, or meta-analysis
3. Publication date: 1986 or later

### 5.8 Example Structured Abstract

```json
{
  "ekat_id": "EKAT-2016-001",
  "icd10_code": "K21",
  "objective": "To evaluate the efficacy of rikkunshito for functional dyspepsia",
  "design": "DB-RCT",
  "setting": "Multicenter, outpatient clinics in Japan",
  "participants": "Adult patients with functional dyspepsia (n=247)",
  "intervention": "Rikkunshito extract granules 7.5g/day for 8 weeks",
  "control": "Placebo granules 7.5g/day for 8 weeks",
  "main_outcome": "Global patient assessment at week 8",
  "results": "Responder rate: 33.6% (rikkunshito) vs 23.8% (placebo), p=0.04",
  "conclusion": "Rikkunshito is effective for functional dyspepsia",
  "kampo_discussion": true,
  "abstractor_comments": "Well-designed study with adequate blinding",
  "reference": "Author A et al. J Gastroenterol 2015;50:xxx-xxx",
  "data_source": "C",
  "language": "en"
}
```

---

## 6. Cross-Reference Systems

### 6.1 Database Interconnections

```
                    +----------------+
                    |    KampoDB     |
                    +----------------+
                           |
    +----------------------+----------------------+
    |                      |                      |
    v                      v                      v
+--------+          +------------+          +----------+
| STORK  |          | KNApSAcK   |          | TradMPD  |
+--------+          +------------+          +----------+
    |                      |                      |
    +----------+-----------+----------+-----------+
               |                      |
               v                      v
        +------------+          +----------+
        |   EKAT     |          |  PubChem |
        +------------+          +----------+
```

### 6.2 Identifier Cross-Reference Table

| Database | Primary ID | External References |
|----------|------------|---------------------|
| KampoDB Formula | String (e.g., "KT") | TradMPD link, STORK code |
| KampoDB Compound | PubChem CID | KNApSAcK C_ID, ChEMBL ID |
| KampoDB Protein | KEGG GENES ID | UniProt ID, Gene Symbol |
| STORK | Formula code | Manufacturer product codes |
| TradMPD | SCH/SCC codes | PubChem CID, NIKKAJI ID |
| KNApSAcK | C_ID | PubChem CID, CAS, KEGG |
| EKAT | EKAT-YEAR-NUM | PubMed ID, DOI |

### 6.3 Common External Databases

| Database | ID Format | Usage |
|----------|-----------|-------|
| PubChem | Numeric CID | Compound identification |
| ChEMBL | CHEMBL + number | Bioactivity data |
| KEGG | Alphanumeric | Pathways, genes, compounds |
| UniProt | Alphanumeric | Protein sequences |
| PDB | 4-character | Protein structures |
| CAS | XXX-XX-X | Chemical registry |
| PubMed | 8-digit PMID | Literature references |

### 6.4 Data Integration Workflow

```
1. Kampo Formula Identification
   KampoDB ID -> STORK verification -> TradMPD linkage

2. Crude Drug Resolution
   Japanese name -> Latin name -> Botanical source (IPNI)

3. Compound Mapping
   KNApSAcK C_ID -> PubChem CID -> ChEMBL ID

4. Target Identification
   Gene Symbol -> KEGG GENES -> UniProt -> PDB

5. Evidence Linking
   EKAT study -> PubMed -> DOI
```

---

## References

1. Sawada R, Iwata M, et al. KampoDB, database of predicted targets and functional annotations of natural medicines. Scientific Reports. 2018;8:11216.

2. Arai M, et al. Standards of Reporting Kampo Products (STORK) in research articles. J Ethnopharmacol. 2017;207:89-95.

3. Afendi FM, et al. KNApSAcK Family Databases: Integrated Metabolite-Plant Species Databases for Multifaceted Plant Research. Plant Cell Physiol. 2012;53(2):e1.

4. Motoo Y, et al. Review of the first 20 years of the Evidence-Based Medicine Committee of the Japan Society for Oriental Medicine. Trad Kampo Med. 2021;8:65-71.

5. University of Toyama. Traditional Medical & Pharmaceutical Database (TradMPD). https://dentomed.toyama-wakan.net/

---

## Appendix A: API Quick Reference (KampoDB)

```
Base URL: https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/

# List endpoints
GET /formula/           # All formulas
GET /crude/             # All crude drugs
GET /compound/          # All compounds
GET /protein/           # All proteins

# Info endpoints
GET /formula/{id}/info
GET /crude/{id}/info
GET /compound/{id}/info
GET /protein/{id}/info

# Relationship endpoints
GET /formula/{id}/crude
GET /formula/{id}/compound
GET /formula/{id}/protein
GET /crude/{id}/formula
GET /crude/{id}/compound
GET /compound/{id}/crude
GET /protein/{id}/compound

# Enrichment analysis
GET /{entity}/{id}/pathway
GET /{entity}/{id}/disease
GET /{entity}/{id}/process
GET /{entity}/{id}/function

# Target prediction
GET /compound/{id}/lbp     # Ligand-based
GET /compound/{id}/sbp     # Structure-based (max 1000)
GET /protein/{id}/lbp
GET /protein/{id}/sbp

# Docking
GET /docking/compound/{cid}/protein/{pid}
```

---

## Appendix B: Data Type Definitions

```typescript
// KampoDB Types
interface KampoFormula {
  id: string;           // e.g., "KT"
  name: string;         // Romanized name
  name_jp: string;      // Japanese kanji
  kana: string;         // Hiragana reading
  synonyms: string;
  synonyms_jp: string;
  link_dento: string;   // TradMPD link
}

interface CrudeDrug {
  id: number;
  name: string;
  name_jp: string;
  kana: string;
  orgn: string;         // Botanical sources
  synonyms_jp: string;
}

interface Compound {
  id: number;           // PubChem CID
  name: string;
}

interface Protein {
  id: number;           // KEGG ID
  name: string;         // Gene symbol
  aliases: string;
  description: string;
}

interface PathwayEnrichment {
  id: string;           // KEGG pathway ID
  name: string;
  count: number;
  fdr: number;
}

interface DiseaseEnrichment {
  id: string;
  name: string;
  count: number;
  fdr: number;
}

// KNApSAcK Types
interface KNApSAcKCompound {
  c_id: string;         // Format: C00000000
  metabolite_name: string;
  molecular_formula: string;
  molecular_weight: number;
  cas_number: string;
  inchi: string;
  smiles: string;
}

// EKAT Types
interface EKATRecord {
  ekat_id: string;
  icd10_code: string;
  objective: string;
  design: 'RCT' | 'DB-RCT' | 'quasi-RCT' | 'RCT-envelope' | 'RCT-crossover' | 'meta-analysis';
  setting: string;
  participants: string;
  intervention: string;
  control: string;
  main_outcome: string;
  results: string;
  conclusion: string;
  kampo_discussion: boolean;
  abstractor_comments: string;
  reference: string;
  data_source: 'C' | 'I' | 'N';
  language: 'ja' | 'en' | 'ko';
}
```

---

*Document generated: 2026-01-18*
*Version: 1.0*
