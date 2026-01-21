# Western Herbal Medicine and Supplements Data Sources

**Document ID:** 43-24-WESTERN-HERBAL
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-00-INDEX.md](./43-00-INDEX.md)

---

## TL;DR

This document catalogs 25+ databases relevant to Western herbal medicine and dietary supplements, organized by accessibility and data quality tiers. Primary integration targets include DSLD (200K+ supplement labels, CC0), Dr. Duke's Phytochemical Database (CC0, bulk CSV), Health Canada LNHPD (REST API), and LOTUS (750K+ natural product occurrences). Commercial sources like NatMed Pro and Examine.com provide highest-quality efficacy ratings for clinical decision support.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary supplement source | DSLD (NIH) | 200K+ labels, full API, CC0 license | Jan 2026 |
| Primary phytochemical source | Dr. Duke's | CC0 license, comprehensive, bulk CSV download | Jan 2026 |
| Natural products aggregation | LOTUS | 750K+ records, open source, MongoDB format | Jan 2026 |
| Efficacy data approach | NatMed Pro + Examine.com | Best evidence grading, requires licensing | Jan 2026 |
| Monograph reference | WHO + Commission E + ESCOP | Regulatory standards, free/low-cost access | Jan 2026 |

---

## Database Catalog

### Tier 1: Primary Open-Access Databases with APIs

#### 1.1 Dietary Supplement Label Database (DSLD)

| Field | Value |
|-------|-------|
| **URL** | https://dsld.od.nih.gov |
| **API** | https://api.ods.od.nih.gov/dsld/v9/ |
| **Content** | Product names, ingredients (medicinal/non-medicinal), amounts, %DV, manufacturer, claims, warnings, label images |
| **Records** | 200,000+ dietary supplement labels sold in USA |
| **License** | CC0 1.0 Universal (Public Domain) |
| **API** | REST API v9.2.0, bulk download available |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5 GB (labels + images) |
| **Data Format** | JSON (primary), CSV, XLSX exports |

**API Endpoints:**
- `GET /v9/brand-products` - Brand product listings
- `GET /v9/browse-brands` - Browse by brand
- `GET /v9/ingredient-groups` - Ingredient categorization
- `GET /v9/browse-products` - Product browsing
- `GET /v9/label/{id}` - Individual label details
- `GET /v9/search-filter` - Filtered search
- `GET /v9/search-filter-histogram` - Search analytics

**Key Schema Fields:**
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

#### 1.2 NIH Office of Dietary Supplements (ODS) API

| Field | Value |
|-------|-------|
| **URL** | https://ods.od.nih.gov/api/ |
| **Content** | Evidence summaries, health professional info, consumer info |
| **Records** | 80+ fact sheets on vitamins, minerals, botanicals |
| **License** | Public domain (US Government work) |
| **API** | REST API (static URLs per fact sheet) |
| **Update Frequency** | As needed |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 MB |
| **Data Format** | XML (with XSD schema), stripped HTML |

**API Structure:**
```
GET /api/index.aspx?resourcename={name}&readinglevel={level}&outputformat={format}

Parameters:
- resourcename: Ashwagandha, Biotin, Calcium, Vitamin D, etc.
- readinglevel: Consumer, Health Professional, Datos en espanol
- outputformat: HTML, XML
```

**Note:** Search functionality planned for future release. Currently one fact sheet per request.

---

#### 1.3 Dr. Duke's Phytochemical and Ethnobotanical Databases

| Field | Value |
|-------|-------|
| **URL** | https://phytochem.nal.usda.gov |
| **Download** | https://data.nal.usda.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases |
| **Content** | Phytochemicals, biological activities, ethnobotanical uses, LD50 toxicity data |
| **Records** | Thousands of plants, compounds, activities |
| **License** | CC0 Public Domain |
| **API** | Web search + CSV bulk download |
| **Update Frequency** | Static (historical dataset) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |
| **Data Format** | CSV (Duke-Source-CSV.zip) |

**Database Tables:**
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

#### 1.4 Health Canada Licensed Natural Health Products Database (LNHPD)

| Field | Value |
|-------|-------|
| **URL** | https://health-products.canada.ca/lnhpd-bdpsnh/ |
| **API** | https://health-products.canada.ca/api/natural-licences/ |
| **Content** | NPN, product name, medicinal/non-medicinal ingredients, recommended uses, warnings, dosage |
| **Records** | All licensed NHPs in Canada |
| **License** | Open Government License - Canada |
| **API** | REST API (JSON/XML) |
| **Update Frequency** | Daily |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB |
| **Data Format** | JSON, XML |

**API Endpoints:**
```
Base: https://health-products.canada.ca/api/natural-licences/

GET /medicinalingredient/?lang=en&type=json
GET /nonmedicinalingredient/?lang=en&type=json
GET /productlicence/?lang=en&type=json

Parameters: id, lang (en/fr), type (json/xml)
```

**Related Resource:** NHPID (Natural Health Products Ingredients Database) - Catalog of approved ingredients for NHP formulations

---

#### 1.5 LOTUS - Natural Products Occurrence Database

| Field | Value |
|-------|-------|
| **URL** | https://lotus.naturalproducts.net |
| **Download** | https://lotus.naturalproducts.net/download/mongo |
| **GitHub** | https://github.com/lotusnprod/ |
| **Content** | Chemical structures, source organisms, literature references |
| **Records** | 750,000+ referenced structure-organism pairs |
| **License** | Open (Wikidata integration) |
| **API** | MongoDB dump download, Wikidata integration |
| **Update Frequency** | Community-driven |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~3 GB |
| **Data Format** | MongoDB, JSON |

**Key Features:**
- Structural search capability
- Taxonomy-oriented queries
- Integration with Wikidata for community curation
- Exports: flat tables, structures

---

#### 1.6 Natural Products Atlas

| Field | Value |
|-------|-------|
| **URL** | https://www.npatlas.org |
| **Download** | https://www.npatlas.org/download |
| **Content** | Compound names, structures, origins, taxonomy, ClassyFire annotations |
| **Records** | 36,545 microbially-derived natural products (v3.0) |
| **License** | CC-BY-NC 4.0 |
| **API** | RESTful API, bulk downloads |
| **Update Frequency** | Periodic releases |
| **Priority** | Tier 1 |
| **Storage Estimate** | ~500 MB |
| **Data Format** | TSV, Excel, JSON, SDF, GraphML |

**API Resources:**
- `/compounds` - Compound entries
- `/references` - Literature references
- `/taxa` - Taxonomic information
- `/networks` - Compound networks

**Note:** Focus on microbial natural products (bacteria, fungi) - complements plant-focused databases.

---

#### 1.7 EMA Herbal Medicines Data

| Field | Value |
|-------|-------|
| **URL** | https://www.ema.europa.eu/en/medicines/download-medicine-data |
| **Content** | Assessment status, monographs, list entries |
| **Records** | EU assessed herbal medicinal products |
| **License** | Open (EU public institution) |
| **API** | Excel download (updated daily), JSON for other medicines |
| **Update Frequency** | Daily |
| **Priority** | Tier 1 |
| **Storage Estimate** | ~100 MB |
| **Data Format** | XLSX (herbal medicines), PDF (full monographs) |

**Note:** JSON endpoint for herbal medicines specifically may not be available yet - table format (Excel) is primary access method. Full monographs available as PDFs.

---

#### 1.8 FooDB - Food Database

| Field | Value |
|-------|-------|
| **URL** | https://foodb.ca |
| **API** | https://foodb.ca/api_doc (Beta) |
| **Downloads** | https://foodb.ca/downloads |
| **Content** | Food compounds, chemical properties, biological activities, dietary sources, health effects |
| **Records** | Largest food constituent database (100+ data fields per compound) |
| **License** | Free for non-commercial; commercial requires permission |
| **API** | API (beta), bulk download |
| **Update Frequency** | Periodic |
| **Priority** | Tier 1 |
| **Storage Estimate** | ~2 GB |
| **Data Format** | JSON (API), various (downloads) |

**Content Types:**
- Food Browse (foods by chemical composition)
- Compound Browse (chemicals by food sources)
- Links to HMDB, PubChem, CHEBI, KEGG

---

#### 1.9 PubChem Natural Products Data

| Field | Value |
|-------|-------|
| **URL** | https://pubchem.ncbi.nlm.nih.gov |
| **API** | https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest |
| **Content** | Chemical structures, properties, bioactivities, safety/toxicity |
| **Records** | World's largest chemical database |
| **License** | Public domain |
| **API** | PUG REST API, PUG SOAP, downloads |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 |
| **Storage Estimate** | Variable (query-based) |
| **Data Format** | JSON, XML, CSV, SDF |

**Natural Products Data Sources in PubChem:**
- LOTUS - Natural products occurrence
- Natural Products Atlas - Microbial NPs
- NPASS - Natural Product Activity and Species Source

**API Example:**
```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/curcumin/JSON
```

---

### Tier 2: Authoritative Subscription/Licensed Sources

#### 2.1 Natural Medicines Database (NatMed Pro)

| Field | Value |
|-------|-------|
| **URL** | https://naturalmedicines.therapeuticresearch.com |
| **Content** | Efficacy ratings, safety, drug interactions, dosing, mechanisms |
| **Records** | Comprehensive natural medicines resource |
| **License** | Commercial subscription |
| **API** | RESTful API (enterprise), web subscription |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | N/A (subscription access) |
| **Data Format** | JSON (API assumed) |
| **Pricing** | Individual: $182/year; API: custom enterprise quote |

**Key Features:**
- Effectiveness ratings (A-F scale)
- Interaction checker
- Nutrient depletion lookup
- Condition-based searches

**Evidence Quality:** **Highest** - systematic evidence reviews, graded efficacy

**API Capabilities (Enterprise):**
- Integration with clinical systems
- Embedding in professional/consumer platforms
- Contact TRC Healthcare for API licensing

---

#### 2.2 Examine.com

| Field | Value |
|-------|-------|
| **URL** | https://examine.com |
| **API Request** | https://examine.com/api-requests/ |
| **Content** | Efficacy grades (A-F), RCT summaries, dosing, safety |
| **Records** | 800+ supplements and health interventions |
| **License** | Commercial; data licensing available |
| **API** | API (in development) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | N/A (subscription access) |
| **Data Format** | TBD (API not yet launched) |
| **Pricing** | Free FAQs; Examine+ for full database access |

**Efficacy Grading System:**
- A: Strong evidence of benefit
- B: Good evidence
- C: Moderate evidence
- D: Weak evidence
- F: Evidence of harm/no benefit

**Evidence Quality:** **High** - graded evidence from RCTs and meta-analyses

**To Request API/Data Access:** Submit form at https://examine.com/api-requests/

---

#### 2.3 HerbMed / HerbMedPro

| Field | Value |
|-------|-------|
| **URL** | http://www.herbmed.org (free) |
| **Pro Version** | Via American Botanical Council |
| **Content** | Categorized research summaries with PubMed links |
| **Records** | 211 herbs (HerbMedPro), 20 herbs (free) |
| **License** | Mixed (free limited, commercial pro) |
| **API** | Web interface; institutional licensing |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | N/A (subscription access) |
| **Data Format** | Web |
| **Pricing** | Free: 20 herbs; Pro: ABC subscription/licensing |

**Categories:**
- Human clinical trials
- Animal studies
- In vitro research
- Adverse effects
- Drug interactions
- Reviews

---

#### 2.4 USP Dietary Supplements Compendium

| Field | Value |
|-------|-------|
| **URL** | https://www.usp.org/dietary-supplements-herbal-medicines |
| **Content** | Quality specifications, acceptance criteria, testing methods |
| **Records** | 650 ingredient/supplement monographs |
| **License** | Commercial |
| **API** | Subscription (USP-NF, DSC) |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | N/A (subscription access) |
| **Data Format** | Proprietary (USP-NF format) |
| **Pricing** | Institutional subscription required |

**Content Includes:**
- 650 dietary ingredient/supplement monographs
- 210 general chapters (testing guidance)
- 120 admission evaluation summaries (safety reviews)

**Evidence Quality:** **Highest** - pharmacopeial standards

**Free Resource:** USP Herbal Medicines Compendium (HMC) - standards for traditional herbal medicines

---

#### 2.5 Memorial Sloan Kettering - About Herbs

| Field | Value |
|-------|-------|
| **URL** | https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs |
| **App** | iOS and Android (About Herbs) |
| **Content** | Purported benefits, side effects, drug interactions, mechanisms |
| **Records** | 290 herbs, botanicals, supplements |
| **License** | Free public access |
| **API** | No public API |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~200 MB (web scraping) |
| **Data Format** | Web/app |

**Evidence Quality:** High - evidence-based, scientific references

**Features:**
- Dual format for healthcare professionals and consumers
- 26+ million hits since 2002 launch

---

### Tier 3: Monograph Collections and Reference Works

#### 3.1 German Commission E Monographs

| Field | Value |
|-------|-------|
| **URL** | https://www.herbalgram.org/resources/commission-e-monographs/ |
| **Browse** | http://cms.herbalgram.org/commissione/index.html |
| **Content** | Approved uses, contraindications, side effects, dosage, drug interactions |
| **Records** | 380 monographs (300+ herbs evaluated) |
| **License** | Copyrighted (ABC translation) |
| **API** | None (web/PDF access) |
| **Update Frequency** | Historical (1978-1994) |
| **Priority** | Tier 3 (Reference) |
| **Storage Estimate** | ~100 MB |
| **Data Format** | Web/PDF |
| **Pricing** | Online: free browse; Book: ~$189 |

**Evidence Quality:** **Regulatory standard** - German government expert committee

**Historical Context:** Established 1978 by German government to evaluate herbal medicine safety/efficacy. Gold standard for herbal monographs.

**Alternative Access:** Internet Archive (borrowable)

---

#### 3.2 ESCOP Monographs

| Field | Value |
|-------|-------|
| **URL** | https://www.escop.com |
| **Monographs** | https://www.escop.com/escop-products/online-monographs/ |
| **Content** | Quality, safety, efficacy data; SPC format |
| **Records** | 107 monographs (86 currently available) |
| **License** | Commercial |
| **API** | None (PDF access) |
| **Update Frequency** | Periodic |
| **Priority** | Tier 3 (Reference) |
| **Storage Estimate** | ~50 MB |
| **Data Format** | PDF |
| **Pricing** | Full access: EUR30/year; free for ESCOP member societies |

**Evidence Quality:** **High** - submitted to EMA, scientific committee review

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

#### 3.3 WHO Monographs on Selected Medicinal Plants

| Field | Value |
|-------|-------|
| **URL** | https://www.who.int/publications/i/item/9241545178 |
| **Download** | https://iris.who.int/handle/10665/42052 |
| **Content** | Quality assurance, identity tests, purity, clinical applications, pharmacology |
| **Records** | 147 monographs across 5 volumes |
| **License** | WHO copyright (free access) |
| **API** | None (PDF download) |
| **Update Frequency** | Historical (completed) |
| **Priority** | Tier 3 (Reference) |
| **Storage Estimate** | ~200 MB |
| **Data Format** | PDF |

**Evidence Quality:** **High** - WHO expert committee review

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

#### 3.4 Cochrane Complementary Medicine Reviews

| Field | Value |
|-------|-------|
| **URL** | https://www.cochranelibrary.com |
| **CAM Portal** | https://cam.cochrane.org |
| **Content** | Meta-analyses, evidence synthesis |
| **Records** | Systematic reviews of herbal/supplement interventions |
| **License** | Subscription (Wiley) |
| **API** | None |
| **Update Frequency** | Continuous |
| **Priority** | Tier 3 (Reference) |
| **Storage Estimate** | N/A (subscription access) |
| **Data Format** | Web, PDF |

**Evidence Quality:** **Gold standard** - systematic review methodology

**Relevant Review Topics:**
- Chinese herbal medicines
- St. John's wort for depression
- Selenium supplements
- Riboflavin for blood pressure
- Various botanical interventions

---

### Tier 4: Phytochemical and Natural Products Databases

#### 4.1 NAPRALERT

| Field | Value |
|-------|-------|
| **URL** | https://napralert.org |
| **Content** | Ethnobotany, chemistry, pharmacology, toxicology, clinical reports |
| **Records** | 200,000+ published studies on natural products |
| **License** | Mixed (limited free, commercial expanded) |
| **API** | None (web search) |
| **Update Frequency** | Continuous (since 1975) |
| **Priority** | Tier 4 |
| **Storage Estimate** | N/A (fee-based access) |
| **Data Format** | Web |
| **Pricing** | Limited free; fee per citation for expanded searches |

**Query Types:**
- Organism-based
- Pharmacology-based
- Compound-based
- Author-based

**Note:** No public API identified. Contact UIC for data licensing inquiries.

---

#### 4.2 IMPPAT - Indian Medicinal Plants Database

| Field | Value |
|-------|-------|
| **URL** | https://cb.imsc.res.in/imppat/ |
| **Content** | 2D/3D structures, ADMET properties, therapeutic uses |
| **Records** | 4,010 Indian medicinal plants, 17,967 phytochemicals |
| **License** | Academic/research |
| **API** | Web search, downloads |
| **Update Frequency** | Periodic |
| **Priority** | Tier 4 |
| **Storage Estimate** | ~1 GB |
| **Data Format** | Web, downloadable structures |

---

#### 4.3 PhytoHub

| Field | Value |
|-------|-------|
| **URL** | https://phytohub.eu |
| **Content** | Plant secondary metabolites, dietary sources, human metabolites |
| **Records** | 1,200 polyphenols, terpenoids, alkaloids + 560 metabolites |
| **License** | Research use |
| **API** | Web search |
| **Update Frequency** | Periodic |
| **Priority** | Tier 4 |
| **Storage Estimate** | ~100 MB |
| **Data Format** | Web |

**Note:** Curated from FooDB, Phenol-Explorer, literature

---

#### 4.4 MPD3 - Medicinal Plant Database for Drug Designing

| Field | Value |
|-------|-------|
| **URL** | https://mpd3.com |
| **Content** | Phytochemicals, activities, targets, structures |
| **Records** | 632 genera, 1,022 plants, 2,295 phytochemicals |
| **License** | Free |
| **API** | Web search, structure library download |
| **Update Frequency** | Static |
| **Priority** | Tier 4 |
| **Storage Estimate** | ~50 MB |
| **Data Format** | Web, downloadable SDF |

**Downloadable:** 2,295 phytochemical structures for molecular docking

---

### Tier 5: Supporting Databases

#### 5.1 Linus Pauling Institute Micronutrient Information Center

| Field | Value |
|-------|-------|
| **URL** | https://lpi.oregonstate.edu/mic |
| **Content** | Evidence summaries, RDAs, health effects, drug interactions |
| **Records** | All essential vitamins, minerals, other nutrients |
| **License** | Free |
| **API** | None (web access) |
| **Update Frequency** | Periodic |
| **Priority** | Tier 5 |
| **Storage Estimate** | ~100 MB (web scraping) |
| **Data Format** | Web (HTML) |

**Languages:** English, Spanish, Japanese

---

#### 5.2 NCCIH Herbs at a Glance

| Field | Value |
|-------|-------|
| **URL** | https://www.nccih.nih.gov/health/herbsataglance |
| **Content** | Brief fact sheets - uses, evidence, safety |
| **Records** | Common herbs and botanicals |
| **License** | Public domain |
| **API** | None (web access) |
| **Update Frequency** | As needed |
| **Priority** | Tier 5 |
| **Storage Estimate** | ~50 MB |
| **Data Format** | Web (HTML) |

---

#### 5.3 MedlinePlus Herbs and Supplements

| Field | Value |
|-------|-------|
| **URL** | https://medlineplus.gov/druginfo/herb_All.html |
| **Content** | Effectiveness, dosing, interactions, safety |
| **Records** | Herbs, supplements, vitamins |
| **License** | Public domain |
| **API** | MedlinePlus Connect API |
| **Update Frequency** | Periodic |
| **Priority** | Tier 5 |
| **Storage Estimate** | ~100 MB |
| **Data Format** | Web, API (for Connect) |

**Note:** As of July 2025, some TRC Natural Medicines content may be unavailable.

---

#### 5.4 Drug Interaction Databases

**RxNav Interaction API (NIH/NLM)**
| Field | Value |
|-------|-------|
| **URL** | https://rxnav.nlm.nih.gov/InteractionAPIs.html |
| **License** | Free API, no license required |
| **Content** | Drug-drug interactions (including some supplements via DrugBank) |

**PHYDGI Database**
| Field | Value |
|-------|-------|
| **Content** | Pharmacokinetic interaction strength levels |
| **Records** | 58 plants, 114 drugs |

**First Databank (FDB)**
| Field | Value |
|-------|-------|
| **License** | Commercial subscription |
| **Content** | Drug-alternative therapy interactions |
| **Features** | Severity levels, alternative therapy category |

---

#### 5.5 Australian TGA ARTG

| Field | Value |
|-------|-------|
| **URL** | https://www.tga.gov.au/products/australian-register-therapeutic-goods-artg |
| **Content** | All licensed therapeutic goods in Australia including complementary medicines |
| **License** | Open |
| **API** | Web search, MedSearch app |
| **Update Frequency** | Continuous |
| **Priority** | Tier 5 |
| **Data Format** | Web |

---

## Quick Reference Table

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

## Integration Recommendations

### Priority 1: Immediate Integration (Open APIs, CC0/Open License)

| Database | Priority Reason | Storage |
|----------|-----------------|---------|
| **DSLD** | 200k+ products, full API, CC0 license | ~5 GB |
| **Dr. Duke's** | CC0 license, comprehensive phytochemical data | ~500 MB |
| **Health Canada LNHPD** | Full REST API, daily updates | ~2 GB |
| **LOTUS** | 750k+ records, open source, MongoDB | ~3 GB |
| **ODS Fact Sheets** | Authoritative, XML API | ~50 MB |

### Priority 2: Consider for Licensed Integration

| Database | Value Proposition | Cost |
|----------|-------------------|------|
| **NatMed Pro** | Best efficacy ratings, drug interactions | $182/yr individual, enterprise API varies |
| **Examine.com** | RCT-based grades, API in development | Examine+ subscription |
| **USP** | Quality specifications (if formulation focus) | Institutional subscription |

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

## Storage Estimates Summary

| Tier | Sources | Estimated Storage |
|------|---------|-------------------|
| Tier 1 (Open) | DSLD, Duke's, LNHPD, LOTUS, NPA, EMA, FooDB, PubChem | ~15 GB |
| Tier 2 (Licensed) | NatMed, Examine, HerbMed, USP, MSK | Subscription-based |
| Tier 3 (Monographs) | Commission E, ESCOP, WHO, Cochrane | ~400 MB |
| Tier 4 (Phytochemical) | NAPRALERT, IMPPAT, PhytoHub, MPD3 | ~1.2 GB |
| Tier 5 (Supporting) | LPI, NCCIH, MedlinePlus, TGA | ~400 MB |
| **Total (Open Sources)** | | **~17 GB** |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX.md](./43-00-INDEX.md) | Parent index |
| [43-20-TRADITIONAL-OVERVIEW.md](./43-20-TRADITIONAL-OVERVIEW.md) | Section context |
| [43-32-SUPPLEMENTS.md](./43-32-SUPPLEMENTS.md) | Overlapping supplement data |
| [43-52-NATURAL-PRODUCTS.md](./43-52-NATURAL-PRODUCTS.md) | Related natural products |
| [43-53-TARGETS-LINKING.md](./43-53-TARGETS-LINKING.md) | Compound-target relationships |
| [44-ARCHITECTURE.md](../44-ARCHITECTURE.md) | Informs database design |
| [45-DATA-MODEL.md](../45-DATA-MODEL.md) | Informs entity structure |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/interventions-western-herbal.md |
