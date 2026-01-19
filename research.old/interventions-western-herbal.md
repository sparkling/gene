# Western Herbal Medicine and Supplements Data Sources

**Research Date:** January 2026
**Purpose:** Comprehensive inventory of databases for genetics/health knowledge base platform

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Tier 1: Primary Open-Access Databases with APIs](#tier-1-primary-open-access-databases-with-apis)
3. [Tier 2: Authoritative Subscription/Licensed Sources](#tier-2-authoritative-subscriptionlicensed-sources)
4. [Tier 3: Monograph Collections & Reference Works](#tier-3-monograph-collections--reference-works)
5. [Tier 4: Phytochemical & Natural Products Databases](#tier-4-phytochemical--natural-products-databases)
6. [Tier 5: Supporting Databases](#tier-5-supporting-databases)
7. [Data Integration Recommendations](#data-integration-recommendations)
8. [Appendix: Quick Reference Table](#appendix-quick-reference-table)

---

## Executive Summary

This document catalogs **25+ databases** relevant to Western herbal medicine and dietary supplements. Sources are organized by accessibility and data quality tiers:

| Tier | Description | Best For |
|------|-------------|----------|
| **Tier 1** | Open APIs, bulk download, CC0/open license | Primary data integration |
| **Tier 2** | Subscription/commercial with APIs | Clinical decision support, efficacy data |
| **Tier 3** | Reference monographs, limited programmatic access | Evidence synthesis, safety data |
| **Tier 4** | Phytochemical/compound databases | Mechanism research, drug interactions |
| **Tier 5** | Supporting/specialized resources | Gap filling, validation |

**Top Recommendations for Integration:**
1. **DSLD** (NIH) - 200,000+ supplement labels, free API, CC0 license
2. **Dr. Duke's Phytochemical Database** - CC0 license, CSV bulk download
3. **Health Canada LNHPD** - Free REST API, daily updates
4. **LOTUS** - 750,000+ natural product occurrences, open source
5. **Natural Products Atlas** - RESTful API, CC-BY-NC license

---

## Tier 1: Primary Open-Access Databases with APIs

### 1.1 Dietary Supplement Label Database (DSLD)

| Attribute | Details |
|-----------|---------|
| **URL** | https://dsld.od.nih.gov |
| **API** | https://api.ods.od.nih.gov/dsld/v9/ |
| **Maintainer** | NIH Office of Dietary Supplements |
| **Coverage** | 200,000+ dietary supplement labels sold in USA |
| **Content** | Product names, ingredients (medicinal/non-medicinal), amounts, %DV, manufacturer, claims, warnings, label images |
| **Access Method** | REST API (v9.2.0), bulk download available |
| **Data Format** | JSON (primary), CSV, XLSX exports |
| **Licensing** | CC0 1.0 Universal (Public Domain) |
| **Full Content** | Yes - complete label information |
| **Evidence Quality** | Label claims only (no efficacy ratings) |

**API Endpoints:**
```
GET /v9/brand-products
GET /v9/browse-brands
GET /v9/ingredient-groups
GET /v9/browse-products
GET /v9/label/{id}
GET /v9/search-filter
GET /v9/search-filter-histogram
```

**Schema/Key Fields:**
- `dsldId` - Unique identifier
- `productName` - Brand product name
- `brandName` - Manufacturer brand
- `servingSize`, `servingUnit`
- `ingredients[]` - Array with name, amount, unit, %DV
- `claims[]` - Health claims
- `warnings[]` - Safety warnings
- `labelImages[]` - URLs to label scans

**Contact:** ODScomments@mail.nih.gov (for full database download)

---

### 1.2 NIH Office of Dietary Supplements (ODS) API

| Attribute | Details |
|-----------|---------|
| **URL** | https://ods.od.nih.gov/api/ |
| **Maintainer** | NIH Office of Dietary Supplements |
| **Coverage** | 80+ fact sheets on vitamins, minerals, botanicals |
| **Content** | Evidence summaries, health professional info, consumer info |
| **Access Method** | REST API (static URLs per fact sheet) |
| **Data Format** | XML (with XSD schema), stripped HTML |
| **Licensing** | Public domain (US Government work) |
| **Full Content** | Yes - complete fact sheet text |
| **Evidence Quality** | High - references scientific literature |

**API Structure:**
```
GET /api/index.aspx?resourcename={name}&readinglevel={level}&outputformat={format}

Parameters:
- resourcename: Ashwagandha, Biotin, Calcium, Vitamin D, etc.
- readinglevel: Consumer, Health Professional, Datos en espa�ol
- outputformat: HTML, XML
```

**Note:** Search functionality planned for future release. Currently one fact sheet per request.

---

### 1.3 Dr. Duke's Phytochemical and Ethnobotanical Databases

| Attribute | Details |
|-----------|---------|
| **URL** | https://phytochem.nal.usda.gov |
| **Bulk Download** | https://data.nal.usda.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases |
| **Maintainer** | USDA National Agricultural Library |
| **Coverage** | Extensive (thousands of plants, compounds, activities) |
| **Content** | Phytochemicals, biological activities, ethnobotanical uses, LD50 toxicity data |
| **Access Method** | Web search + CSV bulk download |
| **Data Format** | CSV (Duke-Source-CSV.zip) |
| **Licensing** | CC0 Public Domain |
| **Full Content** | Yes - complete database tables |
| **Evidence Quality** | Literature-based, with references |

**Database Tables (from Data Dictionary):**
- Plants (scientific/common names, families)
- Chemicals (compounds, CAS numbers)
- Activities (biological effects)
- Ethnobotanical uses
- Chemical-plant relationships
- Activity-chemical relationships

**Key Files:**
- `Duke-Source-CSV.zip` - Raw database tables
- `DrDukesDatabaseDataDictionary-prelim.csv` - Column descriptions

---

### 1.4 Health Canada Licensed Natural Health Products Database (LNHPD)

| Attribute | Details |
|-----------|---------|
| **URL** | https://health-products.canada.ca/lnhpd-bdpsnh/ |
| **API** | https://health-products.canada.ca/api/natural-licences/ |
| **Maintainer** | Health Canada (NNHPD) |
| **Coverage** | All licensed NHPs in Canada |
| **Content** | NPN, product name, medicinal/non-medicinal ingredients, recommended uses, warnings, dosage |
| **Access Method** | REST API (JSON/XML), daily updates |
| **Data Format** | JSON, XML |
| **Licensing** | Open Government License - Canada |
| **Full Content** | Yes - complete licensing information |
| **Evidence Quality** | Products assessed for safety, efficacy, quality |

**API Endpoints:**
```
Base: https://health-products.canada.ca/api/natural-licences/

GET /medicinalingredient/?lang=en&type=json
GET /nonmedicinalingredient/?lang=en&type=json
GET /productlicence/?lang=en&type=json

Parameters: id, lang (en/fr), type (json/xml)
```

**Related Resource:**
- **NHPID** (Natural Health Products Ingredients Database): Catalog of approved ingredients for NHP formulations

---

### 1.5 LOTUS - Natural Products Occurrence Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://lotus.naturalproducts.net |
| **Download** | https://lotus.naturalproducts.net/download/mongo |
| **GitHub** | https://github.com/lotusnprod/ |
| **Maintainer** | Open science community project |
| **Coverage** | 750,000+ referenced structure-organism pairs |
| **Content** | Chemical structures, source organisms, literature references |
| **Access Method** | MongoDB dump download, Wikidata integration |
| **Data Format** | MongoDB, JSON |
| **Licensing** | Open (Wikidata integration) |
| **Full Content** | Yes - full database |
| **Evidence Quality** | Literature-referenced |

**Key Features:**
- Structural search capability
- Taxonomy-oriented queries
- Integration with Wikidata for community curation
- Exports: flat tables, structures

---

### 1.6 Natural Products Atlas

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.npatlas.org |
| **Download** | https://www.npatlas.org/download |
| **Maintainer** | Simon Fraser University / Academic consortium |
| **Coverage** | 36,545 microbially-derived natural products (v3.0) |
| **Content** | Compound names, structures, origins, taxonomy, ClassyFire annotations |
| **Access Method** | RESTful API, bulk downloads |
| **Data Format** | TSV, Excel, JSON, SDF, GraphML |
| **Licensing** | CC-BY-NC 4.0 |
| **Full Content** | Yes - complete compound data |
| **Evidence Quality** | Curated from 1,347+ papers |

**API Resources:**
- `/compounds` - Compound entries
- `/references` - Literature references
- `/taxa` - Taxonomic information
- `/networks` - Compound networks

**Note:** Focus on microbial natural products (bacteria, fungi) - complements plant-focused databases.

---

### 1.7 EMA Herbal Medicines Data

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ema.europa.eu/en/medicines/download-medicine-data |
| **Maintainer** | European Medicines Agency (HMPC) |
| **Coverage** | EU assessed herbal medicinal products |
| **Content** | Assessment status, monographs, list entries |
| **Access Method** | Excel download (updated daily), JSON for other medicines |
| **Data Format** | XLSX (herbal medicines), JSON (general medicines) |
| **Licensing** | Open (EU public institution) |
| **Full Content** | Assessment outcomes; full monographs as PDF |
| **Evidence Quality** | High - HMPC scientific assessment |

**Note:** JSON endpoint for herbal medicines specifically may not be available yet - table format (Excel) is primary access method. Full monographs available as PDFs.

---

### 1.8 FooDB - Food Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://foodb.ca |
| **API** | https://foodb.ca/api_doc (Beta) |
| **Downloads** | https://foodb.ca/downloads |
| **Maintainer** | Wishart Research Group, University of Alberta |
| **Coverage** | Largest food constituent database |
| **Content** | Food compounds, chemical properties, biological activities, dietary sources, health effects |
| **Access Method** | API (beta), bulk download |
| **Data Format** | API (JSON), Downloads (various) |
| **Licensing** | Free for non-commercial; commercial requires permission |
| **Full Content** | Yes - 100+ data fields per compound |
| **Evidence Quality** | Integrated from multiple sources (HMDB, PubChem, etc.) |

**Content Types:**
- Food Browse (foods by chemical composition)
- Compound Browse (chemicals by food sources)
- Links to HMDB, PubChem, CHEBI, KEGG

---

### 1.9 PubChem Natural Products Data

| Attribute | Details |
|-----------|---------|
| **URL** | https://pubchem.ncbi.nlm.nih.gov |
| **API** | https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest |
| **Maintainer** | NIH/NCBI |
| **Coverage** | World's largest chemical database |
| **Content** | Chemical structures, properties, bioactivities, safety/toxicity |
| **Access Method** | PUG REST API, PUG SOAP, downloads |
| **Data Format** | JSON, XML, CSV, SDF |
| **Licensing** | Public domain |
| **Full Content** | Yes |
| **Evidence Quality** | Aggregated from 800+ data sources |

**Natural Products Data Sources in PubChem:**
- LOTUS - Natural products occurrence
- Natural Products Atlas - Microbial NPs
- NPASS - Natural Product Activity and Species Source

**API Example:**
```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/curcumin/JSON
```

---

## Tier 2: Authoritative Subscription/Licensed Sources

### 2.1 Natural Medicines Database (NatMed Pro)

| Attribute | Details |
|-----------|---------|
| **URL** | https://naturalmedicines.therapeuticresearch.com |
| **Maintainer** | TRC Healthcare |
| **Coverage** | Comprehensive natural medicines resource |
| **Content** | Efficacy ratings, safety, drug interactions, dosing, mechanisms |
| **Access Method** | RESTful API (enterprise), web subscription |
| **Data Format** | API (JSON assumed) |
| **Pricing** | Individual: $182/year; API: custom enterprise quote |
| **Licensing** | Commercial subscription |
| **Full Content** | Yes (with subscription) |
| **Evidence Quality** | **Highest** - systematic evidence reviews, graded efficacy |

**Key Features:**
- Effectiveness ratings (A-F scale)
- Interaction checker
- Nutrient depletion lookup
- Condition-based searches

**API Capabilities (Enterprise):**
- Integration with clinical systems
- Embedding in professional/consumer platforms
- Contact TRC Healthcare for API licensing

---

### 2.2 Examine.com

| Attribute | Details |
|-----------|---------|
| **URL** | https://examine.com |
| **API Request** | https://examine.com/api-requests/ |
| **Maintainer** | Examine Inc. |
| **Coverage** | 800+ supplements and health interventions |
| **Content** | Efficacy grades (A-F), RCT summaries, dosing, safety |
| **Access Method** | API (in development), web subscription |
| **Data Format** | TBD (API not yet launched) |
| **Pricing** | Free FAQs; Examine+ for full database access |
| **Licensing** | Commercial; data licensing available |
| **Full Content** | Partial free, full with subscription |
| **Evidence Quality** | **High** - graded evidence from RCTs and meta-analyses |

**Efficacy Grading System:**
- A: Strong evidence of benefit
- B: Good evidence
- C: Moderate evidence
- D: Weak evidence
- F: Evidence of harm/no benefit

**To Request API/Data Access:** Submit form at https://examine.com/api-requests/

---

### 2.3 ConsumerLab.com

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.consumerlab.com |
| **Maintainer** | ConsumerLab.com, LLC |
| **Coverage** | 7,000+ products tested, 1,000+ brands |
| **Content** | Independent lab testing, quality verification, contamination data |
| **Access Method** | Web subscription only (no public API) |
| **Data Format** | Web reports |
| **Pricing** | Individual: $69/year; Institutional: contact for quote |
| **Licensing** | Commercial; licensing available for proprietary data |
| **Full Content** | Yes (with subscription) |
| **Evidence Quality** | **High** - independent laboratory verification |

**Testing Includes:**
- Label accuracy
- Contamination (heavy metals, pesticides)
- Disintegration/dissolution
- Product purity

**Data Access:** Contact ConsumerLab directly for licensing arrangements - no public API.

---

### 2.4 HerbMed / HerbMedPro

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.herbmed.org (free) |
| **Pro Version** | Via American Botanical Council |
| **Maintainer** | American Botanical Council |
| **Coverage** | 211 herbs (HerbMedPro) |
| **Content** | Categorized research summaries with PubMed links |
| **Access Method** | Web interface; institutional licensing |
| **Data Format** | Web |
| **Pricing** | Free: 20 herbs; Pro: ABC subscription/licensing |
| **Licensing** | Mixed (free limited, commercial pro) |
| **Full Content** | Limited free; full with subscription |
| **Evidence Quality** | Literature-based with categorization |

**Categories:**
- Human clinical trials
- Animal studies
- In vitro research
- Adverse effects
- Drug interactions
- Reviews

---

### 2.5 USP Dietary Supplements Compendium

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.usp.org/dietary-supplements-herbal-medicines |
| **Maintainer** | United States Pharmacopeia |
| **Coverage** | 650 ingredient/supplement monographs |
| **Content** | Quality specifications, acceptance criteria, testing methods |
| **Access Method** | Subscription (USP-NF, DSC) |
| **Data Format** | Proprietary (USP-NF format) |
| **Pricing** | Institutional subscription required |
| **Licensing** | Commercial |
| **Full Content** | Yes (with subscription) |
| **Evidence Quality** | **Highest** - pharmacopeial standards |

**Content Includes:**
- 650 dietary ingredient/supplement monographs
- 210 general chapters (testing guidance)
- 120 admission evaluation summaries (safety reviews)

**Free Resource:** USP Herbal Medicines Compendium (HMC) - standards for traditional herbal medicines

---

### 2.6 Memorial Sloan Kettering - About Herbs

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs |
| **App** | iOS and Android (About Herbs) |
| **Maintainer** | Memorial Sloan Kettering Cancer Center |
| **Coverage** | 290 herbs, botanicals, supplements |
| **Content** | Purported benefits, side effects, drug interactions, mechanisms |
| **Access Method** | Free web, mobile app; no public API |
| **Data Format** | Web/app |
| **Licensing** | Free public access |
| **Full Content** | Yes |
| **Evidence Quality** | High - evidence-based, scientific references |

**Dual Format:** Content formatted for both healthcare professionals and consumers.

**Statistics:** 26+ million hits since 2002 launch.

---

## Tier 3: Monograph Collections & Reference Works

### 3.1 German Commission E Monographs

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.herbalgram.org/resources/commission-e-monographs/ |
| **Browse** | http://cms.herbalgram.org/commissione/index.html |
| **Maintainer** | American Botanical Council (English translation) |
| **Coverage** | 380 monographs (300+ herbs evaluated) |
| **Content** | Approved uses, contraindications, side effects, dosage, drug interactions |
| **Access Method** | Web (free browsing), book purchase |
| **Data Format** | Web/PDF |
| **Pricing** | Online: free browse; Book: ~$189 |
| **Licensing** | Copyrighted (ABC translation) |
| **Full Content** | Yes (web/book) |
| **Evidence Quality** | **Regulatory standard** - German government expert committee |

**Historical Context:** Established 1978 by German government to evaluate herbal medicine safety/efficacy. Gold standard for herbal monographs.

**Alternative Access:** Internet Archive (borrowable)

---

### 3.2 ESCOP Monographs

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.escop.com |
| **Monographs** | https://www.escop.com/escop-products/online-monographs/ |
| **Maintainer** | European Scientific Cooperative on Phytotherapy |
| **Coverage** | 107 monographs (86 currently available) |
| **Content** | Quality, safety, efficacy data; SPC format |
| **Access Method** | Online subscription, individual purchase |
| **Data Format** | PDF |
| **Pricing** | Full access: EUR30/year; free for ESCOP member societies |
| **Licensing** | Commercial |
| **Full Content** | Yes (with subscription) |
| **Evidence Quality** | **High** - submitted to EMA, scientific committee review |

**Monograph Structure (SPC Format):**
- Therapeutic indications
- Posology and method of administration
- Contraindications
- Special warnings
- Interactions
- Pregnancy/lactation
- Pharmacological properties
- Preclinical safety data

---

### 3.3 WHO Monographs on Selected Medicinal Plants

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.who.int/publications/i/item/9241545178 |
| **Download** | https://iris.who.int/handle/10665/42052 |
| **Maintainer** | World Health Organization |
| **Coverage** | 147 monographs across 5 volumes |
| **Content** | Quality assurance, identity tests, purity, clinical applications, pharmacology |
| **Access Method** | Free PDF download |
| **Data Format** | PDF |
| **Licensing** | WHO copyright (free access) |
| **Full Content** | Yes |
| **Evidence Quality** | **High** - WHO expert committee review |

**Volumes:**
| Volume | Year | Monographs | Pages |
|--------|------|------------|-------|
| Vol. 1 | 1999 | 28 | 295 |
| Vol. 2 | 2003 | 30 | 357 |
| Vol. 3 | 2007 | 31 | 390 |
| Vol. 4 | 2009 | 28 | 456 |
| NIS Special | 2010 | 30 | 450 |

**Monograph Structure:**
1. **Part 1:** Pharmacopoeial summary (botanical features, identity tests, purity, chemical assays)
2. **Part 2:** Clinical applications (pharmacology, contraindications, warnings, adverse reactions)

---

### 3.4 Cochrane Complementary Medicine Reviews

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.cochranelibrary.com |
| **CAM Portal** | https://cam.cochrane.org |
| **Maintainer** | Cochrane Collaboration |
| **Coverage** | Systematic reviews of herbal/supplement interventions |
| **Content** | Meta-analyses, evidence synthesis |
| **Access Method** | Subscription (some free access via national provision) |
| **Data Format** | Web, PDF |
| **Licensing** | Subscription (Wiley) |
| **Full Content** | Yes (with access) |
| **Evidence Quality** | **Gold standard** - systematic review methodology |

**Relevant Review Topics:**
- Chinese herbal medicines
- St. John's wort for depression
- Selenium supplements
- Riboflavin for blood pressure
- Various botanical interventions

---

## Tier 4: Phytochemical & Natural Products Databases

### 4.1 NAPRALERT

| Attribute | Details |
|-----------|---------|
| **URL** | https://napralert.org |
| **Maintainer** | University of Illinois Chicago |
| **Coverage** | 200,000+ published studies on natural products |
| **Content** | Ethnobotany, chemistry, pharmacology, toxicology, clinical reports |
| **Access Method** | Web search (limited free, fee-based full) |
| **Data Format** | Web |
| **Pricing** | Limited free; fee per citation for expanded searches |
| **Licensing** | Mixed (limited free, commercial expanded) |
| **Full Content** | Fee-based |
| **Evidence Quality** | Comprehensive literature indexing since 1975 |

**Query Types:**
- Organism-based
- Pharmacology-based
- Compound-based
- Author-based

**Note:** No public API identified. Contact UIC for data licensing inquiries.

---

### 4.2 IMPPAT - Indian Medicinal Plants Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://cb.imsc.res.in/imppat/ |
| **Maintainer** | Institute of Mathematical Sciences, Chennai |
| **Coverage** | 4,010 Indian medicinal plants, 17,967 phytochemicals |
| **Content** | 2D/3D structures, ADMET properties, therapeutic uses |
| **Access Method** | Web search, downloads |
| **Data Format** | Web, downloadable structures |
| **Licensing** | Academic/research |
| **Full Content** | Yes |
| **Evidence Quality** | Curated from traditional medicine literature |

---

### 4.3 PhytoHub

| Attribute | Details |
|-----------|---------|
| **URL** | https://phytohub.eu |
| **Maintainer** | European research consortium |
| **Coverage** | 1,200 polyphenols, terpenoids, alkaloids + 560 metabolites |
| **Content** | Plant secondary metabolites, dietary sources, human metabolites |
| **Access Method** | Web search |
| **Data Format** | Web |
| **Licensing** | Research use |
| **Full Content** | Yes |
| **Evidence Quality** | Curated from FooDB, Phenol-Explorer, literature |

---

### 4.4 MPD3 - Medicinal Plant Database for Drug Designing

| Attribute | Details |
|-----------|---------|
| **URL** | https://mpd3.com |
| **Maintainer** | Academic |
| **Coverage** | 632 genera, 1,022 plants, 2,295 phytochemicals |
| **Content** | Phytochemicals, activities, targets, structures |
| **Access Method** | Web search, structure library download |
| **Data Format** | Web, downloadable SDF |
| **Licensing** | Free |
| **Full Content** | Yes |
| **Evidence Quality** | Literature-referenced with PubMed IDs |

**Downloadable:** 2,295 phytochemical structures for molecular docking

---

## Tier 5: Supporting Databases

### 5.1 Linus Pauling Institute Micronutrient Information Center

| Attribute | Details |
|-----------|---------|
| **URL** | https://lpi.oregonstate.edu/mic |
| **Maintainer** | Oregon State University |
| **Coverage** | All essential vitamins, minerals, other nutrients |
| **Content** | Evidence summaries, RDAs, health effects, drug interactions |
| **Access Method** | Free web access |
| **Data Format** | Web (HTML) |
| **Licensing** | Free |
| **Full Content** | Yes |
| **Evidence Quality** | High - scientific literature review |

**Languages:** English, Spanish, Japanese

---

### 5.2 NCCIH Herbs at a Glance

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.nccih.nih.gov/health/herbsataglance |
| **Maintainer** | NIH National Center for Complementary and Integrative Health |
| **Coverage** | Common herbs and botanicals |
| **Content** | Brief fact sheets - uses, evidence, safety |
| **Access Method** | Free web access |
| **Data Format** | Web (HTML) |
| **Licensing** | Public domain |
| **Full Content** | Summary level |
| **Evidence Quality** | NIH quality standards |

---

### 5.3 MedlinePlus Herbs and Supplements

| Attribute | Details |
|-----------|---------|
| **URL** | https://medlineplus.gov/druginfo/herb_All.html |
| **Maintainer** | NIH National Library of Medicine |
| **Coverage** | Herbs, supplements, vitamins |
| **Content** | Effectiveness, dosing, interactions, safety |
| **Access Method** | Free web; MedlinePlus Connect API |
| **Data Format** | Web, API (for Connect) |
| **Licensing** | Public domain |
| **Full Content** | Yes |
| **Evidence Quality** | High - curated by NLM |

**Note:** As of July 2025, some TRC Natural Medicines content may be unavailable.

---

### 5.4 Drug Interaction Databases

#### RxNav Interaction API (NIH/NLM)
- **URL:** https://rxnav.nlm.nih.gov/InteractionAPIs.html
- **Access:** Free API, no license required
- **Coverage:** Drug-drug interactions (including some supplements via DrugBank)

#### PHYDGI Database
- **URL:** Herb-drug interaction database
- **Coverage:** 58 plants, 114 drugs
- **Content:** Pharmacokinetic interaction strength levels

#### First Databank (FDB)
- **Access:** Commercial subscription
- **Coverage:** Drug-alternative therapy interactions
- **Features:** Severity levels, alternative therapy category

---

### 5.5 Australian TGA ARTG

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.tga.gov.au/products/australian-register-therapeutic-goods-artg |
| **Maintainer** | Therapeutic Goods Administration (Australia) |
| **Coverage** | All licensed therapeutic goods in Australia including complementary medicines |
| **Access Method** | Web search, MedSearch app |
| **Data Format** | Web |
| **Licensing** | Open |

---

## Data Integration Recommendations

### Priority 1: Immediate Integration (Open APIs, CC0/Open License)

| Database | Priority Reason |
|----------|-----------------|
| **DSLD** | 200k+ products, full API, CC0 license |
| **Dr. Duke's** | CC0 license, comprehensive phytochemical data |
| **Health Canada LNHPD** | Full REST API, daily updates |
| **LOTUS** | 750k+ records, open source, MongoDB |
| **ODS Fact Sheets** | Authoritative, XML API |

### Priority 2: Consider for Licensed Integration

| Database | Value Proposition |
|----------|-------------------|
| **NatMed Pro** | Best efficacy ratings, drug interactions |
| **Examine.com** | RCT-based grades, API in development |
| **USP** | Quality specifications (if formulation focus) |

### Priority 3: Reference/Validation Sources

| Database | Use Case |
|----------|----------|
| **WHO Monographs** | Safety/efficacy validation |
| **German Commission E** | Regulatory precedent |
| **Cochrane** | Evidence synthesis |
| **MSK About Herbs** | Cancer-specific interactions |

### Data Mapping Considerations

**Common Identifiers to Link Across Databases:**
- Plant scientific names (Latin binomials)
- CAS Registry Numbers (chemicals)
- InChIKey (chemical structures)
- PubChem CID
- UNII (FDA Substance Registration System)

**Key Relationships to Model:**
1. Plant/Organism -> Contains -> Compound
2. Compound -> Has -> Biological Activity
3. Supplement Product -> Contains -> Ingredient(s)
4. Ingredient -> Interacts With -> Drug
5. Intervention -> Evidence For -> Health Outcome

---

## Appendix: Quick Reference Table

| Database | API | Bulk Download | License | Efficacy Data | Drug Interactions |
|----------|-----|---------------|---------|---------------|-------------------|
| DSLD | Yes | Yes | CC0 | No | No |
| ODS API | Yes | No | Public | Partial | Partial |
| Dr. Duke's | No | Yes (CSV) | CC0 | Partial | No |
| Health Canada LNHPD | Yes | Yes | Open Gov | No | No |
| LOTUS | No | Yes | Open | No | No |
| Natural Products Atlas | Yes | Yes | CC-BY-NC | No | No |
| EMA Herbal | No | Yes (Excel) | Open | Yes | Yes |
| NatMed Pro | Yes* | No | Commercial | **Yes** | **Yes** |
| Examine | Soon | No | Commercial | **Yes** | Partial |
| ConsumerLab | No | No | Commercial | Quality only | No |
| Commission E | No | No | Copyright | **Yes** | **Yes** |
| ESCOP | No | PDF | Commercial | **Yes** | **Yes** |
| WHO Monographs | No | PDF | WHO | **Yes** | **Yes** |
| NAPRALERT | No | No | Mixed | Partial | Partial |
| FooDB | Yes* | Yes | Mixed | No | No |
| PubChem | Yes | Yes | Public | Via BioAssay | No |
| MSK About Herbs | No | No | Free | Partial | **Yes** |
| Cochrane | No | No | Subscription | **Yes** | No |
| USP | No | No | Commercial | Standards | No |
| LPI MIC | No | No | Free | **Yes** | **Yes** |

*API with restrictions or enterprise licensing

---

## Sources and References

### Primary Sources Consulted
- [NIH Office of Dietary Supplements API](https://ods.od.nih.gov/api/)
- [Dietary Supplement Label Database](https://dsld.od.nih.gov)
- [Dr. Duke's Phytochemical Database](https://phytochem.nal.usda.gov)
- [Health Canada LNHPD API Documentation](https://health-products.canada.ca/api/documentation/lnhpd-documentation-en.html)
- [LOTUS Natural Products Online](https://lotus.naturalproducts.net)
- [Natural Products Atlas](https://www.npatlas.org)
- [EMA Download Medicine Data](https://www.ema.europa.eu/en/medicines/download-medicine-data)
- [TRC Healthcare NatMed Pro](https://trchealthcare.com/product/natmed-pro/)
- [Examine.com API Requests](https://examine.com/api-requests/)
- [ConsumerLab About](https://www.consumerlab.com/about/)
- [American Botanical Council Commission E](https://www.herbalgram.org/resources/commission-e-monographs/)
- [ESCOP Online Monographs](https://www.escop.com/escop-products/online-monographs/)
- [WHO Monographs on Medicinal Plants](https://www.who.int/publications/i/item/9241545178)
- [NAPRALERT UIC](https://pharmacognosy.pharmacy.uic.edu/napralert/)
- [FooDB](https://foodb.ca)
- [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
- [Memorial Sloan Kettering About Herbs](https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs)
- [Cochrane Complementary Medicine](https://cam.cochrane.org)
- [USP Dietary Supplements](https://www.usp.org/dietary-supplements-herbal-medicines)
- [Linus Pauling Institute MIC](https://lpi.oregonstate.edu/mic)

---

*Document last updated: January 2026*
