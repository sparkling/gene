# Western Herbal Medicine and Supplement Database Data Models

This document provides comprehensive data models, schemas, and field descriptions for major Western herbal medicine and supplement databases.

---

## Table of Contents

1. [DSLD (Dietary Supplement Label Database)](#1-dsld-dietary-supplement-label-database)
2. [Dr. Duke's Phytochemical and Ethnobotanical Databases](#2-dr-dukes-phytochemical-and-ethnobotanical-databases)
3. [Health Canada LNHPD](#3-health-canada-lnhpd)
4. [LOTUS Initiative](#4-lotus-initiative)
5. [Natural Products Atlas](#5-natural-products-atlas)
6. [FooDB](#6-foodb)

---

## 1. DSLD (Dietary Supplement Label Database)

**Base URL:** `https://api.ods.od.nih.gov/dsld/`
**Current Version:** 9.4.0
**License:** CC0 1.0 Universal

### 1.1 API Endpoints

| Endpoint | Method | Description | Required Parameters |
|----------|--------|-------------|---------------------|
| `/v9/brand-products` | GET | Retrieve label info by brand | `q` (brand name) |
| `/v9/browse-brands` | GET | List brands by name/keyword/letter | `method` |
| `/v9/browse-products` | GET | List labels by product name | `method` |
| `/v9/ingredient-groups` | GET | Retrieve ingredient groups | `term`, `method` |
| `/v9/label/{id}` | GET | Get complete label data | `id` (path) |
| `/v9/search-filter` | GET | Complex multi-filter search | `q` |
| `/v9/search-filter-histogram` | GET | Histogram of labels over time | `q` |
| `/version` | GET | API version information | None |

### 1.2 Common Query Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `from` | Integer | 0 | Pagination start offset |
| `size` | Integer | 1000 | Number of records returned |
| `q` | String | - | Query term (use quotes for exact, * for wildcard) |
| `sort_by` | String | - | `_score`, `entryDate`, `fullName.keyword` |
| `sort_order` | String | - | Sort direction |
| `status` | Integer | - | 0=off-market, 1=on-market, 2=all |
| `date_start` | String | - | YYYY format |
| `date_end` | String | - | YYYY format |

### 1.3 Search Filter Parameters

| Parameter | Description | Example Values |
|-----------|-------------|----------------|
| `product_name` | Comma-separated names | - |
| `product_type` | LanguaL codes | a1305, a1306, a1326, a1310, a1302, a1299, a1316, a1315, a1317, a1309, a1325 |
| `ingredient_name` | Comma-separated names | - |
| `apply_synonyms` | Include alternate names | Boolean |
| `ingredient_category` | Category filter | amino_acid, botanical, enzyme, vitamin, mineral |
| `brand` | Comma-separated brands | - |
| `target_group` | Target demographic codes | p0250, p0192, p0266, p0253 |
| `supplement_form` | Form codes | e0164, e0159, e0161, e0155, e0176, e0165, e0174, e0162, e0172, e0177 |
| `claim_type` | Claim type codes | p0065, p0265, p0124, p0264, p0115, p0276 |
| `label_claim` | Dietary claims | Comma-separated |

### 1.4 Label Data Model

#### Primary Label Object

```json
{
  "id": "Integer - Unique DSLD identifier",
  "fullName": "String - Complete product name",
  "brandName": "String - Brand name",
  "upcSku": "String - UPC/SKU code",
  "bundleName": "String (optional) - Bundle product name",
  "nhanesId": "String (optional) - NHANES identifier",
  "servingsPerContainer": "String (optional)",
  "hasOuterCarton": "Boolean",
  "percentDvFootnote": "String - Daily value footnote",
  "targetGroups": "Array<String> - e.g., ['Vegetarian', 'Adult (18-50 Years)']",
  "entryDate": "String - Date entered (YYYY-MM-DD)",
  "offMarket": "Integer - 0 or 1"
}
```

#### Nested Objects

**Contact:**
```json
{
  "contactId": "Integer",
  "text": "String - Contact description",
  "types": "Array<String> - Contact types",
  "contactDetails": {
    "address": "String",
    "phone": "String",
    "email": "String",
    "web": "String"
  }
}
```

**netContents:**
```json
{
  "order": "BigDecimal",
  "quantity": "BigDecimal",
  "unit": "String - e.g., Tablet(s), mL",
  "display": "String - Formatted display text"
}
```

**servingSize:**
```json
{
  "order": "BigDecimal",
  "minQuantity": "BigDecimal",
  "maxQuantity": "BigDecimal",
  "minDailyServings": "BigDecimal",
  "maxDailyServings": "BigDecimal",
  "unit": "String",
  "inSFB": "Boolean"
}
```

**physicalState:**
```json
{
  "langualCode": "String - e.g., 'E0155'",
  "langualCodeDescription": "String - e.g., 'Tablet or Pill'"
}
```

**productType:**
```json
{
  "langualCode": "String - e.g., 'A1325'",
  "langualCodeDescription": "String - e.g., 'Other Combinations'"
}
```

### 1.5 Ingredient Structure

#### Label_ingredientRows

```json
{
  "order": "Integer - Position in ingredient list",
  "ingredientId": "Integer - Unique ingredient ID",
  "name": "String - Ingredient name",
  "category": "String - vitamin, mineral, botanical, enzyme, amino_acid",
  "ingredientGroup": "String - Classification group",
  "uniiCode": "String - Unique Ingredient Identifier",
  "alternateNames": "Array<String>",
  "quantity": "Array<Label_quantity>",
  "forms": "Array<forms>",
  "nestedRows": "Array<Label_nestedRows>"
}
```

#### Label_quantity

```json
{
  "servingSizeOrder": "Integer",
  "servingSizeQuantity": "Integer",
  "operator": "String - '=', '<', '>'",
  "quantity": "Integer - Amount",
  "unit": "String - mg, mcg, IU, etc.",
  "dailyValueTargetGroup": [
    {
      "name": "String - Target group name",
      "operator": "String",
      "percent": "Integer - Percent daily value"
    }
  ],
  "servingSizeUnit": "String - Tablet(s), Capsule(s), etc."
}
```

#### Label_otheringredients

```json
{
  "text": "String (optional) - Free text description",
  "ingredients": [
    {
      "ingredientId": "Integer",
      "name": "String",
      "alternateNames": "Array<String>"
    }
  ]
}
```

### 1.6 Additional Structures

**labelStatements:**
```json
{
  "type": "String - FDA Statement, General Statements, etc.",
  "notes": "String"
}
```

**claims:**
```json
{
  "langualCode": "String",
  "langualCodeDescription": "String"
}
```

**events:**
```json
{
  "date": "String",
  "type": "String - 'Off Market', 'Date entered into DSLD'"
}
```

### 1.7 Search Response Models

**SearchResult:**
```json
{
  "hits": [
    {
      "_index": "String - Elasticsearch index",
      "_type": "String - '_doc'",
      "_id": "String - Unique identifier",
      "_score": "BigDecimal - Match relevance",
      "_source": "Object - Full label data"
    }
  ],
  "stats": {
    "count": "Integer",
    "pct": "Number"
  }
}
```

**HistogramEntry:**
```json
{
  "key_as_string": "String - ISO 8601 date",
  "key": "Integer - Unix timestamp (ms)",
  "doc_count": "Integer - Labels added on date"
}
```

### 1.8 Example Response

```json
{
  "id": 25,
  "fullName": "Vitamins For The Hair",
  "brandName": "Nature's Bounty",
  "upcSku": "0 74312 02100 8",
  "servingsPerContainer": null,
  "hasOuterCarton": false,
  "entryDate": "2012-01-25",
  "offMarket": 0,
  "ingredientRows": [
    {
      "order": 1,
      "ingredientId": 280871,
      "name": "Pantothenic Acid",
      "category": "vitamin",
      "ingredientGroup": "Pantothenic Acid (Vitamin B5)",
      "uniiCode": "19F5HK2737",
      "quantity": [
        {
          "servingSizeOrder": 1,
          "servingSizeQuantity": 1,
          "operator": "=",
          "quantity": 100,
          "unit": "mg",
          "dailyValueTargetGroup": [
            {
              "name": "Adults and children 4 or more years of age",
              "operator": "=",
              "percent": 1000
            }
          ]
        }
      ]
    }
  ]
}
```

---

## 2. Dr. Duke's Phytochemical and Ethnobotanical Databases

**URL:** https://phytochem.nal.usda.gov
**Data Downloads:** https://data.nal.usda.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases
**License:** Creative Commons CCZero (CC0)

### 2.1 Overview

The database contains approximately 49,788 indexed entries covering:
- Plants and their chemical profiles
- Phytochemical compounds
- Biological activities
- Ethnobotanical uses

### 2.2 Available Data Files

| File | Description | Format |
|------|-------------|--------|
| `Duke-Source-CSV.zip` | Raw database tables archive | ZIP (CSV files) |
| `DrDukesDatabaseDataDictionary-prelim.csv` | Column definitions for all tables | CSV |

### 2.3 Core Database Tables (Inferred Schema)

#### PLANTS Table

| Column | Type | Description |
|--------|------|-------------|
| `FNID` | Integer | Unique plant identifier |
| `ScientificName` | String | Botanical name (Genus species) |
| `CommonName` | String | Common English name |
| `Family` | String | Plant family |
| `Part` | String | Plant part used |
| `Notes` | Text | Additional information |

#### CHEMICALS Table

| Column | Type | Description |
|--------|------|-------------|
| `ChemicalID` | Integer | Unique chemical identifier |
| `ChemicalName` | String | Compound name |
| `CASNumber` | String | Chemical Abstracts Service number |
| `MolecularFormula` | String | Chemical formula |
| `MolecularWeight` | Decimal | Molecular weight (g/mol) |
| `SMILES` | String | SMILES notation (where available) |

#### ACTIVITIES Table

| Column | Type | Description |
|--------|------|-------------|
| `ActivityID` | Integer | Unique activity identifier |
| `ActivityName` | String | Biological activity name |
| `ActivityType` | String | Category of activity |
| `Dosage` | String | Effective dosage |
| `Reference` | String | Literature citation |

#### PLANT_CHEMICAL (Relationship Table)

| Column | Type | Description |
|--------|------|-------------|
| `FNFKEY` | Integer | Foreign key to PLANTS |
| `ChemicalID` | Integer | Foreign key to CHEMICALS |
| `Concentration` | Decimal | Concentration value |
| `ConcentrationUnit` | String | Unit (ppm, %, mg/g) |
| `PlantPart` | String | Specific plant part |
| `Reference` | String | Source citation |

#### CHEMICAL_ACTIVITY (Relationship Table)

| Column | Type | Description |
|--------|------|-------------|
| `ChemicalID` | Integer | Foreign key to CHEMICALS |
| `ActivityID` | Integer | Foreign key to ACTIVITIES |
| `MinDosage` | Decimal | Minimum effective dose |
| `MaxDosage` | Decimal | Maximum effective dose |
| `DosageUnit` | String | Unit of measurement |
| `LD50` | Decimal | Lethal dose 50% (toxicity) |
| `Reference` | String | Literature source |

#### ETHNOBOTANY Table

| Column | Type | Description |
|--------|------|-------------|
| `EthnoID` | Integer | Unique ethnobotanical record ID |
| `FNID` | Integer | Foreign key to PLANTS |
| `Use` | String | Traditional use description |
| `Culture` | String | Cultural/geographic origin |
| `Reference` | String | Documentation source |

### 2.4 Data Relationships

```
PLANTS (1) ----< PLANT_CHEMICAL >---- (N) CHEMICALS
                                              |
                                              |
                         CHEMICAL_ACTIVITY >----< ACTIVITIES

PLANTS (1) ----< ETHNOBOTANY (N)
```

### 2.5 Search Capabilities

- Plant search by scientific or common name
- Chemical search with activity filtering
- Biological activity queries
- Ethnobotanical use searches
- LD50 toxicity data display

### 2.6 Export Formats

- CSV (comma-separated values)
- Excel spreadsheet
- PDF reports

---

## 3. Health Canada LNHPD

**Base URL:** `https://health-products.canada.ca/api/natural-licences/`
**Documentation:** https://health-products.canada.ca/api/documentation/lnhpd-documentation-en.html

### 3.1 API Endpoints

| Endpoint | Required Params | Optional Params | Description |
|----------|-----------------|-----------------|-------------|
| `/medicinalingredient/` | None | id, page, lang, type | Active ingredients |
| `/nonmedicinalingredient/` | id | lang, type | Inactive ingredients |
| `/productdose/` | id | lang, type | Dosage information |
| `/productlicence/` | id | lang, type | License details |
| `/productpurpose/` | None | id, page, lang, type | Intended uses |
| `/productrisk/` | None | id, page, lang, type | Risk information |
| `/productroute/` | id | lang, type | Administration routes |

### 3.2 Common Parameters

| Parameter | Type | Values | Default | Description |
|-----------|------|--------|---------|-------------|
| `lang` | String | en, fr | en | Response language |
| `type` | String | json, xml | json | Response format |
| `id` | Integer | - | - | LNHPD product ID |
| `page` | Integer | 1-n | 1 | Pagination (paginated endpoints) |

### 3.3 Response Schemas

#### Product Licence Object

```json
{
  "lnhpd_id": "Integer - Unique product identifier",
  "licence_number": "String - NPN or DIN-HM number",
  "licence_date": "String - Date of initial licensing",
  "revised_date": "String - Most recent revision date",
  "time_receipt": "String - Receipt timestamp",
  "date_start": "String - Product start date",
  "product_name_id": "Integer - Product name identifier",
  "product_name": "String - Product brand name",
  "dosage_form": "String - Capsule, Tablet, Liquid, etc.",
  "company_id": "Integer - Company identifier",
  "company_name_id": "Integer - Company name ID",
  "company_name": "String - Licence holder name",
  "sub_submission_type_code": "String - Submission type code",
  "sub_submission_type_desc": "String - Submission type description",
  "flag_primary_name": "Boolean - Primary name indicator",
  "flag_product_status": "Boolean - Active status",
  "flag_attested_monograph": "Boolean - Monograph attestation"
}
```

#### Medicinal Ingredient Object

```json
{
  "lnhpd_id": "Integer - Product identifier",
  "ingredient_name": "String - Ingredient name",
  "potency_amount": "Decimal - Potency value",
  "potency_unit_of_measure": "String - Unit (mg, mcg, IU)",
  "potency_constituent": "String - Active constituent name",
  "quantity": "Decimal - Amount per serving",
  "quantity_minimum": "Decimal - Minimum quantity",
  "quantity_maximum": "Decimal - Maximum quantity",
  "quantity_unit_of_measure": "String - Quantity unit",
  "ratio_numerator": "Integer - Extract ratio numerator",
  "ratio_denominator": "Integer - Extract ratio denominator",
  "dried_herb_equivalent": "Decimal - DHE value",
  "dhe_unit_of_measure": "String - DHE unit",
  "extract_type_desc": "String - Extract type description",
  "source_material": "String - Source material/plant part"
}
```

#### Non-Medicinal Ingredient Object

```json
{
  "lnhpd_id": "Integer - Product identifier",
  "ingredient_name": "String - Non-medicinal ingredient name"
}
```

#### Product Dose Object

```json
{
  "lnhpd_id": "Integer - Product identifier",
  "dose_id": "Integer - Dose record identifier",
  "population_type_desc": "String - Target population (Adults, Children)",
  "age": "Integer - Specific age (if applicable)",
  "age_minimum": "Integer - Minimum age",
  "age_maximum": "Integer - Maximum age",
  "uom_type_desc_age": "String - Age unit (years, months)",
  "quantity_dose": "Decimal - Dose quantity",
  "quantity_dose_minimum": "Decimal - Minimum dose",
  "quantity_dose_maximum": "Decimal - Maximum dose",
  "uom_type_desc_quantity_dose": "String - Dose unit",
  "frequency": "Integer - Doses per time period",
  "frequency_minimum": "Integer - Minimum frequency",
  "frequency_maximum": "Integer - Maximum frequency",
  "uom_type_desc_frequency": "String - Frequency unit (daily, weekly)"
}
```

#### Product Purpose Object

```json
{
  "text_id": "Integer - Purpose text identifier",
  "lnhpd_id": "Integer - Product identifier",
  "purpose": "String - Intended use/health claim text"
}
```

#### Product Risk Object

```json
{
  "lnhpd_id": "Integer - Product identifier",
  "risk_id": "Integer - Risk record identifier",
  "risk_type_desc": "String - Risk category",
  "sub_risk_type_desc": "String - Risk subcategory",
  "risk_text": "String - Warning/caution text"
}
```

#### Product Route Object

```json
{
  "lnhpd_id": "Integer - Product identifier",
  "route_id": "Integer - Route identifier",
  "route_type_desc": "String - Administration route (Oral, Topical)"
}
```

### 3.4 Pagination Structure

```json
{
  "limit": "Integer - Results per page",
  "page": "Integer - Current page number",
  "total": "Integer - Total result count",
  "next": "String/null - Next page URL",
  "previous": "String/null - Previous page URL"
}
```

### 3.5 Data Extract Files

Available compressed ASCII text files (~30 MB uncompressed):

| File | Description |
|------|-------------|
| `NHP_Products` | Core product information |
| `NHP_Products_purpose` | Intended uses/claims |
| `NHP_Prod_recommended_dose` | Dosage guidelines |
| `NHP_Products_risk` | Safety/risk information |
| `NHP_Companies` | Manufacturer details |
| `NHP_Medicinal_ingredients` | Active components |
| `NHP_Route` | Administration methods |
| `NHP_Nonmedicinal_ingredients` | Inactive components |
| `NHP_Product_names` | Product nomenclature |
| `All_NHP_files` | Combined package (27 MB) |

---

## 4. LOTUS Initiative

**Web Portal:** https://lotus.naturalproducts.net
**Primary Host:** Wikidata
**Publication:** eLife 2022;11:e70780 (DOI: 10.7554/eLife.70780)
**License:** CC0 (via Wikidata)

### 4.1 Overview

LOTUS (Linked Open Universal Tripartite Source) contains 750,000+ referenced structure-organism pairs, establishing relationships between chemical structures and their biological sources.

### 4.2 Core Data Model (Three Objects)

#### Chemical Structure Object

```json
{
  "smiles_canonical": "String - Canonical SMILES (Wikidata P233)",
  "smiles_isomeric": "String - Isomeric SMILES (Wikidata P2017)",
  "inchi": "String - International Chemical Identifier",
  "inchikey": "String - Hashed InChI (Wikidata P235)",
  "molecular_formula": "String - Chemical formula",
  "molecular_weight": "Decimal - Molecular weight",
  "structure_sanitized": "Boolean - Salt removal/dimer resolution applied"
}
```

#### Biological Organism Object

```json
{
  "taxon_name": "String - Standardized taxonomic name",
  "taxon_db": "String - Source database (NCBI, GBIF, etc.)",
  "taxon_id": "String - Identifier in source database",
  "wikidata_id": "String - Wikidata Q-number",
  "open_tree_of_life_id": "String - OTL ID (Wikidata P9157)",
  "rank": "String - Taxonomic rank"
}
```

#### Reference Object

```json
{
  "title": "String - Article title",
  "doi": "String - Digital Object Identifier",
  "pmid": "String - PubMed ID",
  "pmcid": "String - PubMed Central ID",
  "wikidata_reference_id": "String - Wikidata Q-number"
}
```

### 4.3 Structure-Organism Pair Format

```json
{
  "structure": {
    "smiles": "String",
    "inchikey": "String",
    "wikidata_id": "String"
  },
  "organism": {
    "taxon_name": "String",
    "taxon_id": "String",
    "wikidata_id": "String"
  },
  "reference": {
    "doi": "String",
    "wikidata_id": "String"
  },
  "validation_status": "String - manually_validated, automatically_validated"
}
```

### 4.4 Wikidata Property Mappings

| Property | Code | Description |
|----------|------|-------------|
| Found in taxon | P703 | Links compound to organism |
| Canonical SMILES | P233 | Canonical SMILES representation |
| Isomeric SMILES | P2017 | Stereochemistry-preserving SMILES |
| InChIKey | P235 | Hashed InChI identifier |
| Stated in | P248 | Reference to bibliographic source |
| Open Tree of Life ID | P9157 | OTL taxonomy identifier |

### 4.5 Download Files (Zenodo)

| File | Size | Description |
|------|------|-------------|
| `structure_metadata.tsv.gz` | 95.8 MB | Chemical structure data |
| `organism_metadata.tsv.gz` | 7.1 MB | Taxonomic data |
| `reference_metadata.tsv.gz` | 114.7 MB | Bibliographic references |

### 4.6 TSV File Schemas

#### structure_metadata.tsv

| Column | Type | Description |
|--------|------|-------------|
| structure_wikidata | String | Wikidata Q-number |
| structure_smiles_canonical | String | Canonical SMILES |
| structure_smiles_isomeric | String | Isomeric SMILES |
| structure_inchi | String | Full InChI |
| structure_inchikey | String | InChIKey |
| structure_molecular_formula | String | Molecular formula |
| structure_exact_mass | Decimal | Exact molecular mass |

#### organism_metadata.tsv

| Column | Type | Description |
|--------|------|-------------|
| organism_wikidata | String | Wikidata Q-number |
| organism_name | String | Scientific name |
| organism_taxonomy_db | String | Source taxonomy database |
| organism_taxonomy_id | String | ID in source database |
| organism_taxonomy_rank | String | Taxonomic rank |
| organism_taxonomy_parent | String | Parent taxon |

#### reference_metadata.tsv

| Column | Type | Description |
|--------|------|-------------|
| reference_wikidata | String | Wikidata Q-number |
| reference_doi | String | DOI |
| reference_pmid | String | PubMed ID |
| reference_pmcid | String | PMC ID |
| reference_title | String | Article title |

### 4.7 Data Processing Pipeline

1. **Harmonization** - Aggregates 40+ heterogeneous resources
2. **Processing** - Standardizes chemical, organism, and reference fields
3. **Validation** - Automatic filtering (~97% accuracy)
4. **Dissemination** - Wikidata + mirror portal

### 4.8 SPARQL Query Examples

```sparql
# Find compounds in a specific organism
SELECT ?compound ?compoundLabel ?smiles WHERE {
  ?compound p:P703 ?stmt .
  ?stmt ps:P703 wd:Q158695 .  # Example: Ginkgo biloba
  ?compound wdt:P233 ?smiles .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

## 5. Natural Products Atlas

**URL:** https://www.npatlas.org
**API Base:** https://www.npatlas.org/api/v1/
**API Docs:** https://www.npatlas.org/api/v1/docs
**License:** CC BY-NC 4.0
**Current Version:** 2024_09 (36,545 compounds)

### 5.1 API Rate Limits

- **Standard tier:** 20 requests/minute
- **API key tier:** Higher limits (contact support@npatlas.org)

### 5.2 API Endpoints

#### Compound Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/compound/{npaid}` | GET | Full compound record |
| `/compound/{npaid}/mol` | GET | Molecular structure file |
| `/compound/{npaid}/node` | GET | Network node data |
| `/compound/{npaid}/cluster` | GET | Cluster data |
| `/compound/{npaid}/npclassifier` | GET | NPClassifier taxonomy |
| `/compound/{npaid}/classyfire` | GET | ClassyFire taxonomy |
| `/compounds/` | GET | List abbreviated compounds |
| `/compounds/full/` | GET | List complete records |
| `/compounds/{prop}/` | GET | Get properties (names, inchikeys, npaids) |

#### Search Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/compounds/massSearch` | POST | Search by molecular mass |
| `/compounds/structureSearch` | POST | Structure-based search |
| `/compounds/taxonSearch` | POST | Search by taxonomy |
| `/compounds/basicSearch` | POST | Multi-parameter filtered search |
| `/compounds/basicSearchExport` | POST | Export search results |
| `/compounds/advancedSearch` | POST | Complex conditional queries |

#### Reference Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/reference/journals` | GET | List all journals |
| `/reference/{doi}` | GET | Full reference record |
| `/reference/{doi}/compounds/` | GET | Compounds from reference |

#### Taxonomy Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/taxon/` | GET | List taxonomy ranks |
| `/taxon/{rank}/` | GET | Get taxa for rank |
| `/taxon/{rank}/{taxon}` | GET | Full taxon description |
| `/taxon/{taxon_id}` | GET | Get taxon by ID |
| `/taxon/search` | POST | Search taxon by name |

### 5.3 Query Parameters

#### Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `skip` | Integer | 0 | Pagination offset |
| `limit` | Integer | 10 | Results per page (max 100) |
| `name` | String | - | SQL LIKE search (% or * wildcards) |
| `origin_type` | Enum | - | `bacterium` or `fungi` |
| `orderby` | Enum | - | npaid, name, origin_date, similarity |
| `ascending` | Boolean | true | Sort direction |

#### Mass Search Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mass` | Number | - | Target molecular mass |
| `type` | Enum | - | exact_mass, mol_weight, m_plus_h, m_plus_na |
| `range` | Number | 5.0 | Mass tolerance |
| `unit` | Enum | - | ppm, dalton |

#### Structure Search Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structure` | String | - | SMILES, InChI, formula, or SMARTS |
| `method` | Enum | - | sim (similarity), sub (substructure), full (exact) |
| `threshold` | Number | 0.8 | Similarity cutoff (0-1) |
| `stereo` | Boolean | false | Stereochemistry matching |

### 5.4 Response Schemas

#### CompoundBase

```json
{
  "id": "Integer - Internal database ID",
  "npaid": "String - NPA identifier (e.g., NPA012345)",
  "original_name": "String - Compound name from literature",
  "mol_formula": "String - Molecular formula",
  "mol_weight": "Decimal - Molecular weight",
  "exact_mass": "Decimal - Exact mass",
  "inchikey": "String - InChIKey",
  "smiles": "String - SMILES notation",
  "cluster_id": "Integer - Network cluster ID",
  "node_id": "Integer - Network node ID",
  "has_exclusions": "Boolean - Has exclusion records"
}
```

#### CompoundFull (extends CompoundBase)

```json
{
  "...CompoundBase fields...",
  "synonyms": "Array<String> - Alternative names",
  "inchi": "String - Full InChI",
  "m_plus_h": "Decimal - [M+H]+ mass",
  "m_plus_na": "Decimal - [M+Na]+ mass",
  "origin_reference": "Reference - First isolation reference",
  "origin_organism": "Organism - Source organism",
  "syntheses": "Array<Reference> - Synthesis references",
  "reassignments": "Array<Reassignment> - Structure revisions",
  "mol_structures": "Array<MolStructure> - Structure versions",
  "exclusions": "Array<Exclusion> - Exclusion records",
  "external_ids": "Array<ExternalDB> - Cross-references"
}
```

#### CompoundFullExtra (extends CompoundFull)

```json
{
  "...CompoundFull fields...",
  "classyfire": {
    "kingdom": "String",
    "superclass": "String",
    "class": "String",
    "subclass": "String",
    "direct_parent": "String"
  },
  "npclassifier": {
    "isglycoside": "Boolean",
    "class_results": "Array<String>",
    "superclass_results": "Array<String>",
    "pathway_results": "Array<String>"
  }
}
```

#### Organism

```json
{
  "id": "Integer - Organism ID",
  "type": "Enum - bacterium, fungi",
  "genus": "String - Genus name",
  "species": "String - Species name",
  "taxon": "TaxonCompound - Full taxonomic info"
}
```

#### TaxonBase

```json
{
  "id": "Integer - Taxon ID",
  "name": "String - Taxon name",
  "rank": "String - Taxonomic rank",
  "taxon_db": "Enum - lpsn, mycobank",
  "external_id": "String - External DB ID",
  "ncbi_id": "String - NCBI Taxonomy ID"
}
```

#### TaxonFull (extends TaxonBase)

```json
{
  "...TaxonBase fields...",
  "ancestors": "Array<TaxonBase> - Parent taxa",
  "children": "Array<TaxonBase> - Child taxa",
  "compound_count": "Integer - Number of compounds"
}
```

#### Reference

```json
{
  "doi": "String - Digital Object Identifier",
  "pmid": "String - PubMed ID",
  "authors": "String - Author list",
  "title": "String - Article title",
  "journal": "String - Journal name",
  "year": "Integer - Publication year",
  "volume": "String - Journal volume",
  "issue": "String - Journal issue",
  "pages": "String - Page range"
}
```

### 5.5 Advanced Query Schema

```json
{
  "operator": "Enum - eq, neq, gt, gte, lt, lte, contains, ncontains",
  "attribute": "String - QueryParamEnum value",
  "value": "Mixed - Search value"
}
```

**Conditional Expression:**
```json
{
  "conditional": "Enum - AND, OR, ANDNOT, ORNOT",
  "left_expression": "QueryExpr",
  "right_expression": "QueryExpr"
}
```

**Searchable Attributes (QueryParamEnum):**
- Basic: npaid, name, inchi, inchikey, smiles
- Mass: molecular_weight, exact_mass, m_plus_h, m_plus_na, molecular_formula
- Network: cluster_id, node_id
- Flags: has_synonyms, has_gnps, has_mibig, has_npmrd, has_syntheses, has_exclusions, has_reassignments
- NPClassifier: npclassifer_isglycoside, npclassifer_class, npclassifer_superclass, npclassifer_pathway
- ClassyFire: classyfire_kingdom, classyfire_superclass, classyfire_class, classyfire_subclass
- Taxonomy: origin_type, taxon_*
- Reference: reference_*

### 5.6 Download Formats

| Format | Contents |
|--------|----------|
| TSV | Names, structures, origins, metadata |
| Excel | Same as TSV |
| JSON | Full data + ClassyFire/NPClassifier annotations |
| SDF | Chemical structure data format |
| GraphML | Network files for Cytoscape |

### 5.7 Taxonomic Ranks

| Rank | Count (v3.0) |
|------|--------------|
| Domain | 3 |
| Kingdom | 2 |
| Phylum | 27 |
| Class | 65 |
| Order | 171 |
| Family | 427 |
| Genus | 1,178 |

---

## 6. FooDB

**URL:** https://foodb.ca
**API Documentation:** https://foodb.ca/api_doc
**Schema Reference:** https://foodb.ca/schema
**License:** CC BY-NC 4.0
**Current Version:** 1.0

### 6.1 Overview

FooDB is the world's largest database on food constituents, containing 100+ data fields per compound entry covering:
- Nomenclature and structure
- Physico-chemical properties
- Food sources and concentrations
- Sensory attributes (flavor, color, aroma, taste)
- Health effects

### 6.2 Database Schema

#### foods Table

| Column | Type | Description |
|--------|------|-------------|
| `food_id` | Integer | Primary identifier |
| `name` | String | Common food name |
| `name_scientific` | String | Scientific nomenclature |
| `description` | Text | Food description |
| `public_id` | String | FooDB public ID |

#### compounds Table

| Column | Type | Description |
|--------|------|-------------|
| `id` | Integer | Internal compound ID |
| `public_id` | String | FDB ID (e.g., FDB000004) |
| `name` | String | Compound name |
| `export` | Boolean | Active status (1=current, 2=inactive) |
| `description` | Text | Compound description |
| `status` | Integer | Quantification level (0-3) |

#### contents Table (Association)

| Column | Type | Description |
|--------|------|-------------|
| `food_id` | Integer | FK to foods table |
| `source_id` | Integer | FK to compound or nutrient |
| `source_type` | String | "Nutrient" or "Compound" |
| `orig_food_common_name` | String | Original food name |
| `orig_min` | Decimal | Minimum concentration |
| `orig_max` | Decimal | Maximum concentration |
| `orig_content` | Decimal | Average/direct value |
| `orig_unit` | String | Unit (mg/100g, etc.) |

#### nutrients Table

| Column | Type | Description |
|--------|------|-------------|
| `nutrient_id` | Integer | Nutrient identifier |
| `name` | String | Nutrient name |
| `description` | String | Nutrient information |
| `public_id` | String | FooDB nutrient ID |

### 6.3 API Endpoints

#### Authentication

- **Required:** API key via `api_key` parameter
- **Access:** Contact FooDB for credentials
- **Policy:** Unauthorized access results in permanent ban

#### Food Table Endpoint

**GET** `/foods`

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `food_name` | String | Yes | Search term |
| `page` | Integer | Yes | Pagination |
| `api_key` | String | Yes | API key |

**Response Fields:**
```json
{
  "food_id": "Integer",
  "public_id": "String - FooDB ID",
  "food_name": "String",
  "scientific_name": "String",
  "description": "String",
  "food_group": "String",
  "food_sub_group": "String",
  "updated_at": "String - ISO timestamp"
}
```

#### Compound Table Endpoint

**GET** `/compounds`

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `compound_name` | String | No | Search by name |
| `compound_id` | Integer/String | No | Internal or public ID |
| `cas_number` | String | No | CAS registry number |
| `inchikey` | String | No | Molecular InChIKey |
| `filters` | Object | No | max_weight, min_weight |
| `page` | Integer | Yes | Pagination |

**Response Fields:**
```json
{
  "id": "Integer",
  "public_id": "String - FDB ID",
  "name": "String",
  "molecular_formula": "String",
  "weight": "Decimal",
  "cas_number": "String",
  "pubchem_compound_id": "String",
  "hmdb_id": "String",
  "drugbank_id": "String",
  "chebi_id": "String",
  "smiles": "String",
  "inchi": "String",
  "inchikey": "String"
}
```

#### Content Table Endpoint

**GET** `/contents`

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `compound_id` | Integer/String | No* | Compound identifier |
| `food_id` | Integer/String | No* | Food identifier |
| `page` | Integer | Conditional | Required for single-param lookup |
| `page_size` | Integer | No | Results per page (default 10, max 100) |

*One of compound_id or food_id required

**Response Fields:**
```json
{
  "content_id": "Integer",
  "source_id": "Integer",
  "source_type": "String",
  "food_id": "Integer",
  "food_name": "String",
  "preparation_type": "String",
  "citation": "String",
  "citation_type": "String",
  "orig_content": "Decimal",
  "orig_min": "Decimal",
  "orig_max": "Decimal",
  "orig_unit": "String"
}
```

### 6.4 API Response Format

```json
{
  "status": 200,
  "value": ["Array of results"],
  "count": "Integer - Total results",
  "date": "String - Response timestamp",
  "api_version": "String - v1",
  "licence": "String - DEMO or license type"
}
```

### 6.5 Download Formats

| Format | Size | Description |
|--------|------|-------------|
| CSV | 952 MB | Structured data tables |
| XML | 6,438 MB | Full hierarchical data |
| JSON | 87 MB | Compact data format |
| MySQL Dump | 173 MB | Database import file |

### 6.6 Spectroscopic Data Downloads

| Type | Experimental | Predicted |
|------|--------------|-----------|
| LC-MS | 16 MB | 381 MB |
| MS-MS | 52 MB | 1,199 MB |
| NMR Spectra | 448 MB | - |
| NMR FID | 2,116 MB | - |

### 6.7 External Database Links

FooDB cross-references with:
- HMDB (Human Metabolome Database)
- PubChem
- ChEBI
- KEGG
- NCBI Taxonomy

### 6.8 Data Relationships

```
foods (1) ----< contents >---- (1) compounds
                   |
                   +---------- (1) nutrients

compounds ----< external_links >---- PubChem, HMDB, ChEBI, KEGG
```

---

## Summary: ID Systems Comparison

| Database | Primary ID | Format Example | External IDs |
|----------|-----------|----------------|--------------|
| DSLD | DSLD ID | Integer (25) | UPC/SKU, UNII |
| Dr. Duke's | FNID | Integer | CAS Number |
| LNHPD | lnhpd_id | Integer | NPN, DIN-HM |
| LOTUS | Wikidata Q | Q-number (Q123456) | InChIKey, OTL ID |
| NP Atlas | NPAID | NPA + digits (NPA012345) | GNPS, MIBiG, NPMRD |
| FooDB | public_id | FDB + digits (FDB000004) | PubChem, HMDB, ChEBI |

---

## References

1. NIH Office of Dietary Supplements - DSLD API Documentation: https://api.ods.od.nih.gov/dsld/v9/
2. USDA National Agricultural Library - Dr. Duke's Database: https://phytochem.nal.usda.gov
3. Health Canada LNHPD API: https://health-products.canada.ca/api/documentation/lnhpd-documentation-en.html
4. Rutz et al. (2022). The LOTUS initiative for open knowledge management in natural products research. eLife 11:e70780
5. van Santen et al. (2022). The Natural Products Atlas 2.0. Nucleic Acids Research 50(D1):D1317-D1323
6. FooDB Documentation: https://foodb.ca/about

---

*Document generated: 2026-01-18*
