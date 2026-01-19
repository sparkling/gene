# FDA and Regulatory Drug Database Data Models

This document provides comprehensive schema documentation for FDA and regulatory drug databases including OpenFDA, DailyMed, RxNorm, ClinicalTrials.gov, and the FDA Orange Book.

---

## Table of Contents

1. [OpenFDA](#1-openfda)
   - [Drug Event API (FAERS)](#11-drug-event-api-faers)
   - [Drug Label API](#12-drug-label-api)
   - [Drug NDC API](#13-drug-ndc-api)
   - [Drug Enforcement API](#14-drug-enforcement-api)
   - [Query Syntax](#15-openfda-query-syntax)
   - [OpenFDA Harmonized Fields](#16-openfda-harmonized-fields)
2. [DailyMed](#2-dailymed)
   - [SPL XML Schema](#21-spl-xml-schema)
   - [Key Sections](#22-key-sections)
   - [LOINC Codes](#23-loinc-codes)
   - [FTP Download Structure](#24-ftp-download-structure)
3. [RxNorm](#3-rxnorm)
   - [RRF File Format](#31-rrf-file-format)
   - [Concept Types (TTY)](#32-concept-types-tty)
   - [Relationship Types](#33-relationship-types)
   - [REST API](#34-rest-api)
   - [RxClass Integration](#35-rxclass-integration)
4. [ClinicalTrials.gov](#4-clinicaltrialsgov)
   - [Study Data Structure](#41-study-data-structure)
   - [Intervention Types](#42-intervention-types)
   - [Outcome Measures](#43-outcome-measures)
   - [API v2 Endpoints](#44-api-v2-endpoints)
5. [FDA Orange Book](#5-fda-orange-book)
   - [Data File Structure](#51-data-file-structure)
   - [Therapeutic Equivalence Codes](#52-therapeutic-equivalence-codes)
   - [Exclusivity Codes](#53-exclusivity-codes)
   - [Patent Data Format](#54-patent-data-format)
6. [ID System Reference](#6-id-system-reference)

---

## 1. OpenFDA

OpenFDA provides public access to FDA datasets through RESTful APIs. All endpoints return JSON data.

**Base URL**: `https://api.fda.gov`

### 1.1 Drug Event API (FAERS)

The Drug Adverse Event API returns data from the FDA Adverse Event Reporting System (FAERS), containing adverse event and medication error reports.

**Endpoint**: `https://api.fda.gov/drug/event.json`

**Data Source**: FAERS (2004-present, quarterly updates)

**Standard**: ICH E2b/M2 version 2.1

**Coding System**: MedDRA (Medical Dictionary for Regulatory Activities)

#### Schema Fields

##### Top-Level Report Fields

| Field | Type | Description |
|-------|------|-------------|
| `safetyreportid` | string | Unique identifier for the safety report (8-digit: 7 digits + checksum) |
| `safetyreportversion` | string | Version number of the safety report |
| `receivedate` | string | Date FDA received the report (YYYYMMDD) |
| `receivedateformat` | string | Format code for receivedate |
| `receiptdate` | string | Date the most recent information was received |
| `receiptdateformat` | string | Format code for receiptdate |
| `transmissiondate` | string | Date the report was transmitted |
| `transmissiondateformat` | string | Format code for transmissiondate |
| `reporttype` | string | Type of report: 1=Spontaneous, 2=Literature, 3=Study |
| `serious` | string | 1=Serious, 2=Not serious |
| `seriousnessdeath` | string | 1=Patient died |
| `seriousnesshospitalization` | string | 1=Patient hospitalized |
| `seriousnesslifethreatening` | string | 1=Life-threatening |
| `seriousnessdisabling` | string | 1=Disability resulted |
| `seriousnesscongenitalanomali` | string | 1=Congenital anomaly |
| `seriousnessother` | string | 1=Other serious outcome |
| `duplicate` | string | 1=Duplicate report |
| `companynumb` | string | Company's internal report number |
| `authoritynumb` | string | Regulatory authority's case number |
| `occurcountry` | string | Country where event occurred (ISO 3166-1 alpha-2) |
| `primarysourcecountry` | string | Country of primary source |
| `fulfillexpediteriteria` | string | Criteria for expedited reporting |

##### Patient Object (`patient`)

| Field | Type | Description |
|-------|------|-------------|
| `patient.patientonsetage` | string | Patient age at event onset |
| `patient.patientonsetageunit` | string | Unit: 800=Decade, 801=Year, 802=Month, 803=Week, 804=Day, 805=Hour |
| `patient.patientsex` | string | 0=Unknown, 1=Male, 2=Female |
| `patient.patientweight` | string | Patient weight in kg |
| `patient.patientdeathdate` | string | Date of death |
| `patient.patientdeathdateformat` | string | Format code for death date |
| `patient.reaction` | array | Array of reaction objects |
| `patient.drug` | array | Array of drug objects |

##### Reaction Object (`patient.reaction[]`)

| Field | Type | Description |
|-------|------|-------------|
| `reactionmeddrapt` | string | MedDRA preferred term for reaction |
| `reactionmeddraversionpt` | string | MedDRA version used |
| `reactionoutcome` | string | 1=Recovered, 2=Recovering, 3=Not recovered, 4=Recovered with sequelae, 5=Fatal, 6=Unknown |

##### Drug Object (`patient.drug[]`)

| Field | Type | Description |
|-------|------|-------------|
| `drugcharacterization` | string | 1=Suspect, 2=Concomitant, 3=Interacting |
| `medicinalproduct` | string | Drug name as reported |
| `drugbatchnumb` | string | Lot number |
| `drugauthorizationnumb` | string | NDA/ANDA number |
| `drugstructuredosagenumb` | string | Dose amount |
| `drugstructuredosageunit` | string | Dose unit code |
| `drugseparatedosagenumb` | string | Number of separate dosages |
| `drugintervaldosageunitnumb` | string | Interval between doses |
| `drugintervaldosagedefinition` | string | Interval unit |
| `drugcumulativedosagenumb` | string | Cumulative dose |
| `drugcumulativedosageunit` | string | Cumulative dose unit |
| `drugdosagetext` | string | Free-text dosage information |
| `drugdosageform` | string | Dosage form |
| `drugadministrationroute` | string | Route of administration code |
| `drugindication` | string | Indication for use |
| `drugstartdate` | string | Date drug started |
| `drugstartdateformat` | string | Format code |
| `drugenddate` | string | Date drug stopped |
| `drugenddateformat` | string | Format code |
| `drugtreatmentduration` | string | Duration of treatment |
| `drugtreatmentdurationunit` | string | Duration unit |
| `actiondrug` | string | 1=Withdrawn, 2=Dose reduced, 3=Dose increased, 4=Dose not changed, 5=Unknown, 6=Not applicable |
| `drugrecurreadministration` | string | 1=Yes, 2=No, 3=Unknown |
| `drugadditional` | string | 1=Yes, 2=No |
| `activesubstance` | object | Active ingredient information |
| `openfda` | object | OpenFDA harmonized fields |

##### Primary Source Object (`primarysource`)

| Field | Type | Description |
|-------|------|-------------|
| `qualification` | string | 1=Physician, 2=Pharmacist, 3=Other health professional, 4=Lawyer, 5=Consumer |
| `reportercountry` | string | Country code |

##### Sender Object (`sender`)

| Field | Type | Description |
|-------|------|-------------|
| `sendertype` | string | Type of sender |
| `senderorganization` | string | Organization name |

##### Receiver Object (`receiver`)

| Field | Type | Description |
|-------|------|-------------|
| `receivertype` | string | Type of receiver |
| `receiverorganization` | string | Organization name |

#### Example Record

```json
{
  "safetyreportid": "10003214",
  "safetyreportversion": "1",
  "receivedate": "20140101",
  "receivedateformat": "102",
  "serious": "1",
  "seriousnesshospitalization": "1",
  "occurcountry": "US",
  "patient": {
    "patientsex": "2",
    "patientonsetage": "45",
    "patientonsetageunit": "801",
    "reaction": [
      {
        "reactionmeddrapt": "Nausea",
        "reactionoutcome": "1"
      }
    ],
    "drug": [
      {
        "drugcharacterization": "1",
        "medicinalproduct": "ASPIRIN",
        "drugindication": "Pain",
        "openfda": {
          "brand_name": ["ASPIRIN"],
          "generic_name": ["ASPIRIN"],
          "rxcui": ["1191"]
        }
      }
    ]
  }
}
```

---

### 1.2 Drug Label API

Returns data from FDA-approved prescription and OTC drug product labels (package inserts).

**Endpoint**: `https://api.fda.gov/drug/label.json`

**Data Source**: Structured Product Labeling (SPL) submissions

#### Schema Fields

##### Metadata Fields

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Unique document identifier |
| `set_id` | string | Set ID (groups all versions of a label) |
| `version` | string | Label version number |
| `effective_time` | string | Date label became effective (YYYYMMDD) |

##### Clinical Information Sections

| Field | Type | Description |
|-------|------|-------------|
| `indications_and_usage` | array | Approved uses for the drug |
| `dosage_and_administration` | array | Dosing instructions |
| `dosage_forms_and_strengths` | array | Available formulations |
| `contraindications` | array | Conditions where drug should not be used |
| `warnings_and_cautions` | array | Safety warnings (PLR format) |
| `warnings` | array | Safety warnings (older format) |
| `precautions` | array | Precautionary information |
| `adverse_reactions` | array | Reported adverse events |
| `drug_interactions` | array | Drug-drug interactions |
| `use_in_specific_populations` | array | Pregnancy, nursing, pediatric, geriatric use |
| `overdosage` | array | Overdose management |

##### Pharmacology Sections

| Field | Type | Description |
|-------|------|-------------|
| `clinical_pharmacology` | array | Mechanism and pharmacokinetics |
| `mechanism_of_action` | array | How the drug works |
| `pharmacodynamics` | array | Drug effects on body |
| `pharmacokinetics` | array | Absorption, distribution, metabolism, excretion |

##### Additional Sections

| Field | Type | Description |
|-------|------|-------------|
| `description` | array | Drug description and formulation |
| `how_supplied` | array | Available packages and NDC codes |
| `storage_and_handling` | array | Storage requirements |
| `active_ingredient` | array | Active components |
| `inactive_ingredient` | array | Inactive components |
| `purpose` | array | Purpose (OTC drugs) |
| `boxed_warning` | array | Black box warnings |
| `recent_major_changes` | array | Recent significant label changes |
| `clinical_studies` | array | Clinical trial summaries |
| `nonclinical_toxicology` | array | Animal study data |
| `references` | array | Scientific references |
| `patient_medication_information` | array | Patient labeling |
| `spl_product_data_elements` | array | Product data |

##### OTC-Specific Sections

| Field | Type | Description |
|-------|------|-------------|
| `do_not_use` | array | Contraindications for OTC |
| `stop_use` | array | When to stop taking |
| `ask_doctor` | array | When to consult physician |
| `ask_doctor_or_pharmacist` | array | Consultation guidance |
| `keep_out_of_reach_of_children` | array | Safety warning |
| `directions` | array | Usage directions |
| `other_safety_information` | array | Additional safety info |
| `questions` | array | Contact information |

##### OpenFDA Fields

| Field | Type | Description |
|-------|------|-------------|
| `openfda.application_number` | array | NDA/ANDA/BLA number |
| `openfda.brand_name` | array | Brand name(s) |
| `openfda.generic_name` | array | Generic name(s) |
| `openfda.manufacturer_name` | array | Manufacturer name(s) |
| `openfda.product_ndc` | array | Product NDC codes |
| `openfda.product_type` | array | Product type |
| `openfda.route` | array | Route(s) of administration |
| `openfda.substance_name` | array | Active substance(s) |
| `openfda.rxcui` | array | RxNorm concept identifiers |
| `openfda.spl_id` | array | SPL document ID |
| `openfda.spl_set_id` | array | SPL set ID |
| `openfda.pharm_class_epc` | array | Established pharmacologic class |
| `openfda.pharm_class_moa` | array | Mechanism of action class |
| `openfda.pharm_class_pe` | array | Physiologic effect class |
| `openfda.pharm_class_cs` | array | Chemical structure class |
| `openfda.unii` | array | UNII codes |

#### Example Record

```json
{
  "id": "12345678-1234-1234-1234-123456789abc",
  "set_id": "abcdef12-3456-7890-abcd-ef1234567890",
  "version": "3",
  "effective_time": "20230615",
  "indications_and_usage": [
    "DRUG X is indicated for the treatment of moderate to severe pain."
  ],
  "dosage_and_administration": [
    "Adults: Take 1 tablet every 4-6 hours as needed. Do not exceed 6 tablets in 24 hours."
  ],
  "contraindications": [
    "Known hypersensitivity to the active ingredient."
  ],
  "adverse_reactions": [
    "The most common adverse reactions (>=5%) include: nausea, dizziness, headache."
  ],
  "openfda": {
    "application_number": ["NDA012345"],
    "brand_name": ["DRUG X"],
    "generic_name": ["ACETAMINOPHEN"],
    "manufacturer_name": ["EXAMPLE PHARMA INC"],
    "product_ndc": ["12345-6789"],
    "route": ["ORAL"],
    "rxcui": ["161"]
  }
}
```

---

### 1.3 Drug NDC API

Returns data from the National Drug Code Directory, containing information on drugs registered with the FDA.

**Endpoint**: `https://api.fda.gov/drug/ndc.json`

**Data Source**: NDC Directory (updated daily)

#### NDC Structure

The National Drug Code is a 10-digit, 3-segment identifier:
- **Segment 1 (Labeler)**: 4-5 digits - identifies the manufacturer/distributor
- **Segment 2 (Product)**: 3-4 digits - identifies the drug formulation
- **Segment 3 (Package)**: 1-2 digits - identifies the package size

**Configurations**: 4-4-2, 5-3-2, or 5-4-1

#### Schema Fields

##### Product Information

| Field | Type | Description |
|-------|------|-------------|
| `product_id` | string | Unique product identifier |
| `product_ndc` | string | Product NDC (labeler-product segments) |
| `spl_id` | string | SPL document identifier |
| `product_type` | string | HUMAN PRESCRIPTION DRUG, HUMAN OTC DRUG, etc. |
| `finished` | string | "true" = finished drug, "false" = active ingredient |
| `brand_name` | string | Proprietary name |
| `brand_name_base` | string | Base name without suffix |
| `brand_name_suffix` | string | Name suffix |
| `generic_name` | string | Non-proprietary name |
| `dosage_form` | string | TABLET, CAPSULE, SOLUTION, etc. |
| `route` | array | Route(s) of administration |
| `marketing_start_date` | string | Date marketing began (YYYYMMDD) |
| `marketing_end_date` | string | Date marketing ended |
| `marketing_category` | string | NDA, ANDA, BLA, OTC MONOGRAPH, etc. |
| `application_number` | string | FDA application number |
| `labeler_name` | string | Manufacturer/distributor name |
| `dea_schedule` | string | CI, CII, CIII, CIV, CV (controlled substances) |
| `listing_expiration_date` | string | When listing expires |

##### Active Ingredients Object (`active_ingredients[]`)

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Ingredient name |
| `strength` | string | Amount with units (e.g., "500 mg/1") |

##### Packaging Object (`packaging[]`)

| Field | Type | Description |
|-------|------|-------------|
| `package_ndc` | string | Full NDC including package segment |
| `description` | string | Package description |
| `marketing_start_date` | string | Package marketing start date |
| `marketing_end_date` | string | Package marketing end date |
| `sample` | boolean | True if sample package |

##### Pharmacologic Class (`pharm_class[]`)

| Field | Type | Description |
|-------|------|-------------|
| Value | string | Pharmacologic class(es) with type suffix: [EPC], [MoA], [PE], [CS] |

#### Example Record

```json
{
  "product_id": "12345-6789_12345678-1234-1234-1234-123456789012",
  "product_ndc": "12345-6789",
  "product_type": "HUMAN PRESCRIPTION DRUG",
  "brand_name": "METFORMIN HYDROCHLORIDE",
  "generic_name": "METFORMIN HYDROCHLORIDE",
  "dosage_form": "TABLET",
  "route": ["ORAL"],
  "marketing_category": "ANDA",
  "application_number": "ANDA076543",
  "labeler_name": "GENERIC PHARMA INC",
  "active_ingredients": [
    {
      "name": "METFORMIN HYDROCHLORIDE",
      "strength": "500 mg/1"
    }
  ],
  "packaging": [
    {
      "package_ndc": "12345-6789-01",
      "description": "100 TABLET in 1 BOTTLE",
      "marketing_start_date": "20100115"
    }
  ],
  "pharm_class": [
    "Biguanide [EPC]",
    "Biguanides [Chemical/Ingredient]"
  ],
  "openfda": {
    "rxcui": ["861004"],
    "unii": ["9100L32L2N"]
  }
}
```

---

### 1.4 Drug Enforcement API

Returns recall and enforcement data from the FDA Recall Enterprise System (RES).

**Endpoint**: `https://api.fda.gov/drug/enforcement.json`

**Data Source**: FDA Recall Enterprise System (2004-present, weekly updates)

#### Schema Fields

##### Recall Information

| Field | Type | Description |
|-------|------|-------------|
| `recall_number` | string | Unique recall identifier (e.g., D-1234-2023) |
| `event_id` | string | Event identifier |
| `status` | string | Ongoing, Completed, Terminated |
| `classification` | string | Class I (most serious), II, or III (least serious) |
| `recalling_firm` | string | Company initiating recall |
| `city` | string | City of recalling firm |
| `state` | string | State of recalling firm |
| `country` | string | Country of recalling firm |
| `address_1` | string | Street address line 1 |
| `address_2` | string | Street address line 2 |
| `postal_code` | string | Postal/ZIP code |
| `voluntary_mandated` | string | Voluntary, FDA Mandated |
| `initial_firm_notification` | string | How firm notified: Letter, Email, Press Release, etc. |
| `distribution_pattern` | string | Geographic distribution description |

##### Product Information

| Field | Type | Description |
|-------|------|-------------|
| `product_type` | string | Drugs, Devices, Food, etc. |
| `product_description` | string | Description of recalled product |
| `product_quantity` | string | Amount of product recalled |
| `code_info` | string | Lot numbers, expiration dates |
| `more_code_info` | string | Additional code information |
| `product_code` | string | Product code |

##### Recall Dates

| Field | Type | Description |
|-------|------|-------------|
| `recall_initiation_date` | string | When recall began (YYYYMMDD) |
| `center_classification_date` | string | When FDA classified recall |
| `report_date` | string | Date reported (YYYYMMDD) |
| `termination_date` | string | When recall terminated |

##### Reason and Actions

| Field | Type | Description |
|-------|------|-------------|
| `reason_for_recall` | string | Explanation of why recall occurred |

##### Classification Definitions

| Class | Definition |
|-------|------------|
| Class I | Dangerous or defective products that could cause serious health problems or death |
| Class II | Products that might cause temporary health problem or slight threat of serious nature |
| Class III | Products unlikely to cause any adverse health reaction but violate FDA regulations |

#### Example Record

```json
{
  "recall_number": "D-1234-2023",
  "event_id": "87654",
  "status": "Ongoing",
  "classification": "Class II",
  "recalling_firm": "EXAMPLE PHARMACEUTICALS",
  "city": "New York",
  "state": "NY",
  "country": "United States",
  "voluntary_mandated": "Voluntary: Firm Initiated",
  "product_type": "Drugs",
  "product_description": "ASPIRIN 325mg Tablets, 100 count bottle",
  "product_quantity": "50,000 bottles",
  "code_info": "Lot numbers: 12345A, 12345B. Expiration dates: 12/2024, 03/2025",
  "distribution_pattern": "Nationwide distribution in the United States",
  "reason_for_recall": "Failed dissolution specifications",
  "recall_initiation_date": "20230601",
  "report_date": "20230615",
  "openfda": {
    "brand_name": ["ASPIRIN"],
    "generic_name": ["ASPIRIN"]
  }
}
```

---

### 1.5 OpenFDA Query Syntax

#### Basic Query Structure

```
https://api.fda.gov/{endpoint}?search={field}:{term}&limit={n}
```

#### Search Parameter Syntax

| Syntax | Description | Example |
|--------|-------------|---------|
| `field:term` | Simple match | `patient.drug.openfda.brand_name:aspirin` |
| `field:"phrase"` | Exact phrase | `reason_for_recall:"failed dissolution"` |
| `field:term+AND+field:term` | Both conditions | `serious:1+AND+seriousnessdeath:1` |
| `field:term+field:term` | Either condition (OR) | `brand_name:aspirin+brand_name:ibuprofen` |
| `field:[min+TO+max]` | Range (inclusive) | `receivedate:[20200101+TO+20201231]` |
| `field:term*` | Wildcard | `brand_name:child*` |
| `_exists_:field` | Field has value | `_exists_:patient.patientdeathdate` |
| `_missing_:field` | Field is empty | `_missing_:companynumb` |

#### Query Parameters

| Parameter | Description | Max Value |
|-----------|-------------|-----------|
| `search` | Filter by field values | - |
| `count` | Count unique values in field | Returns top 1000 |
| `limit` | Number of records to return | 1000 |
| `skip` | Skip records (pagination) | 25000 |
| `sort` | Sort results | `field:asc` or `field:desc` |

#### Exact Match Suffix

Use `.exact` suffix for counting unique phrases:
```
count=patient.reaction.reactionmeddrapt.exact
```

#### Date Format

Dates use YYYYMMDD format:
```
search=receivedate:[20200101+TO+20201231]
```

#### Example Queries

```bash
# Find aspirin adverse events with hospitalization
https://api.fda.gov/drug/event.json?search=patient.drug.openfda.generic_name:aspirin+AND+seriousnesshospitalization:1&limit=10

# Count reactions for a drug
https://api.fda.gov/drug/event.json?search=patient.drug.openfda.brand_name:lipitor&count=patient.reaction.reactionmeddrapt.exact

# Search drug labels for diabetes
https://api.fda.gov/drug/label.json?search=indications_and_usage:diabetes&limit=25

# Find Class I recalls
https://api.fda.gov/drug/enforcement.json?search=classification:"Class+I"&limit=100
```

---

### 1.6 OpenFDA Harmonized Fields

OpenFDA adds standardized `openfda` fields to records across all endpoints.

| Field | Type | Description |
|-------|------|-------------|
| `application_number` | array | NDA, ANDA, or BLA number |
| `brand_name` | array | Proprietary name(s) |
| `generic_name` | array | Non-proprietary name(s) |
| `manufacturer_name` | array | Manufacturer/labeler name(s) |
| `nui` | array | NDF-RT NUI identifier(s) |
| `package_ndc` | array | Package-level NDC codes |
| `pharm_class_cs` | array | Chemical structure class |
| `pharm_class_epc` | array | Established pharmacologic class |
| `pharm_class_moa` | array | Mechanism of action class |
| `pharm_class_pe` | array | Physiologic effect class |
| `product_ndc` | array | Product-level NDC codes |
| `product_type` | array | Product type classification |
| `route` | array | Administration route(s) |
| `rxcui` | array | RxNorm concept identifier(s) |
| `spl_id` | array | SPL document identifier(s) |
| `spl_set_id` | array | SPL set identifier(s) |
| `substance_name` | array | Active substance name(s) |
| `unii` | array | UNII (Unique Ingredient Identifier) |

---

## 2. DailyMed

DailyMed is the official FDA drug label repository maintained by the National Library of Medicine.

**Website**: https://dailymed.nlm.nih.gov

### 2.1 SPL XML Schema

Structured Product Labeling (SPL) is an HL7 v3 standard for drug labeling in XML format.

#### Document Structure

```xml
<?xml version="1.0" encoding="UTF-8"?>
<document xmlns="urn:hl7-org:v3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <!-- Document header -->
  <id root="document-uuid"/>
  <code code="34391-3" codeSystem="2.16.840.1.113883.6.1" displayName="HUMAN PRESCRIPTION DRUG LABEL"/>
  <title>Drug Name - Package Insert</title>
  <effectiveTime value="YYYYMMDD"/>
  <setId root="set-uuid"/>
  <versionNumber value="n"/>

  <!-- Author/Organization -->
  <author>
    <assignedEntity>
      <representedOrganization>
        <id extension="labeler-code" root="1.3.6.1.4.1.519.1"/>
        <name>Manufacturer Name</name>
      </representedOrganization>
    </assignedEntity>
  </author>

  <!-- Document body -->
  <component>
    <structuredBody>
      <component>
        <section>
          <id root="section-uuid"/>
          <code code="LOINC-code" codeSystem="2.16.840.1.113883.6.1"/>
          <title>Section Title</title>
          <text>Section content...</text>
          <effectiveTime value="YYYYMMDD"/>
        </section>
      </component>
      <!-- Additional sections -->
    </structuredBody>
  </component>
</document>
```

#### Key XML Elements

| Element | Description |
|---------|-------------|
| `<document>` | Root element |
| `<id>` | Unique document identifier (UUID) |
| `<setId>` | Groups all versions of a label |
| `<versionNumber>` | Label version |
| `<effectiveTime>` | Date label became effective |
| `<code>` | LOINC code identifying section type |
| `<title>` | Section title |
| `<text>` | Human-readable content |
| `<author>` | Labeler/sponsor information |
| `<representedOrganization>` | Company information |
| `<manufacturedProduct>` | Product details |
| `<ingredient>` | Active/inactive ingredients |

#### Product Data Elements

```xml
<manufacturedProduct>
  <manufacturedMedicine>
    <code code="NDC-code" codeSystem="2.16.840.1.113883.6.69"/>
    <name>Drug Name</name>
    <formCode code="form-code" codeSystem="2.16.840.1.113883.3.26.1.1"/>
    <asEntityWithGeneric>
      <genericMedicine>
        <name>Generic Name</name>
      </genericMedicine>
    </asEntityWithGeneric>
    <ingredient classCode="ACTIB">
      <ingredientSubstance>
        <code code="UNII-code" codeSystem="2.16.840.1.113883.4.9"/>
        <name>Ingredient Name</name>
      </ingredientSubstance>
    </ingredient>
  </manufacturedMedicine>
</manufacturedProduct>
```

---

### 2.2 Key Sections

#### Prescription Drug Label Sections (PLR Format)

| Section | LOINC | Content |
|---------|-------|---------|
| Highlights of Prescribing Information | 48780-1 | Summary of key information |
| Boxed Warning | 34066-1 | Black box warning |
| Recent Major Changes | 43683-2 | Significant recent changes |
| Indications and Usage | 34067-9 | Approved uses |
| Dosage and Administration | 34068-7 | Dosing instructions |
| Dosage Forms and Strengths | 43678-2 | Available formulations |
| Contraindications | 34070-3 | When not to use |
| Warnings and Precautions | 43685-7 | Safety warnings |
| Adverse Reactions | 34084-4 | Side effects |
| Drug Interactions | 34073-7 | Drug-drug interactions |
| Use in Specific Populations | 43684-0 | Special population guidance |
| Clinical Pharmacology | 34090-1 | Mechanism and PK/PD |
| Nonclinical Toxicology | 43680-8 | Animal studies |
| Clinical Studies | 34092-7 | Human trial summaries |
| How Supplied/Storage and Handling | 34069-5 | Package information |
| Patient Counseling Information | 34076-0 | Patient guidance |

#### Clinical Pharmacology Subsections

| Subsection | LOINC | Content |
|------------|-------|---------|
| Mechanism of Action | 43679-0 | How drug works |
| Pharmacodynamics | 43681-6 | Effects on body |
| Pharmacokinetics | 43682-4 | ADME data |

---

### 2.3 LOINC Codes

Complete reference of LOINC codes for SPL sections:

| LOINC Code | Section Name |
|------------|--------------|
| 34066-1 | BOXED WARNING |
| 34067-9 | INDICATIONS AND USAGE |
| 34068-7 | DOSAGE AND ADMINISTRATION |
| 34069-5 | HOW SUPPLIED |
| 34070-3 | CONTRAINDICATIONS |
| 34071-1 | WARNINGS |
| 34072-9 | GENERAL PRECAUTIONS |
| 34073-7 | DRUG INTERACTIONS |
| 34074-5 | DRUG & OR LABORATORY TEST INTERACTIONS |
| 34075-2 | LABORATORY TESTS |
| 34076-0 | INFORMATION FOR PATIENTS |
| 34077-8 | TERATOGENIC EFFECTS |
| 34078-6 | NONTERATOGENIC EFFECTS |
| 34079-4 | LABOR AND DELIVERY |
| 34080-2 | NURSING MOTHERS |
| 34081-0 | PEDIATRIC USE |
| 34082-8 | GERIATRIC USE |
| 34083-6 | CARCINOGENESIS & MUTAGENESIS & IMPAIRMENT OF FERTILITY |
| 34084-4 | ADVERSE REACTIONS |
| 34085-1 | CONTROLLED SUBSTANCE |
| 34086-9 | ABUSE |
| 34087-7 | DEPENDENCE |
| 34088-5 | OVERDOSAGE |
| 34089-3 | DESCRIPTION |
| 34090-1 | CLINICAL PHARMACOLOGY |
| 34091-9 | ANIMAL PHARMACOLOGY & OR TOXICOLOGY |
| 34092-7 | CLINICAL STUDIES |
| 34093-5 | REFERENCES |
| 42229-5 | SPL UNCLASSIFIED SECTION |
| 43678-2 | DOSAGE FORMS AND STRENGTHS |
| 43679-0 | MECHANISM OF ACTION |
| 43680-8 | NONCLINICAL TOXICOLOGY |
| 43681-6 | PHARMACODYNAMICS |
| 43682-4 | PHARMACOKINETICS |
| 43683-2 | RECENT MAJOR CHANGES |
| 43684-0 | USE IN SPECIFIC POPULATIONS |
| 43685-7 | WARNINGS AND PRECAUTIONS |
| 48780-1 | HIGHLIGHTS OF PRESCRIBING INFORMATION |
| 51945-4 | PACKAGE LABEL.PRINCIPAL DISPLAY PANEL |

---

### 2.4 FTP Download Structure

**FTP URL**: ftp://public.nlm.nih.gov/nlmdata/.dailymed/

#### Available Downloads

| File Type | Description |
|-----------|-------------|
| `dm_spl_release_human_rx_part1.zip` - `partN.zip` | Human prescription drug labels (split files) |
| `dm_spl_release_human_otc_part1.zip` - `partN.zip` | Human OTC drug labels (split files) |
| `dm_spl_release_homeopathic.zip` | Homeopathic products |
| `dm_spl_release_animal.zip` | Animal drug labels |
| `dm_spl_release_remainder.zip` | Bulk ingredients, vaccines, devices |
| `dm_spl_daily_update_YYYYMMDD.zip` | Daily updates |
| `dm_spl_weekly_update_YYYYMMDD.zip` | Weekly updates |
| `dm_spl_monthly_update_YYYYMM.zip` | Monthly updates |

#### Index Files

| File | Description |
|------|-------------|
| `rxnorm_mappings.zip` | SPL Set ID to RxNorm mappings |
| `pharmacologic_class_indexing_spl_files.zip` | Pharmacologic class index |
| `product_info.zip` | Product information index |
| `application_numbers.zip` | NDA/ANDA number mappings |
| `rems_data_files.zip` | REMS information |

#### REST API Endpoints

**Base URL**: `https://dailymed.nlm.nih.gov/dailymed/services/v2/`

| Endpoint | Description |
|----------|-------------|
| `/spls.json` | List all SPLs |
| `/spls/{setId}.xml` | Get SPL by set ID |
| `/drugnames.json` | Search drug names |
| `/drugclasses.json` | Get pharmacologic classes |
| `/ndcs.json` | Get NDC information |
| `/rxcuis.json` | Get RxCUI mappings |

---

## 3. RxNorm

RxNorm is a normalized drug nomenclature system developed by NLM.

**Website**: https://www.nlm.nih.gov/research/umls/rxnorm/

### 3.1 RRF File Format

RxNorm uses the Rich Release Format (RRF) - pipe-delimited text files with UTF-8 encoding.

#### RXNCONSO.RRF (Concepts and Atoms)

Contains all drug names from source vocabularies and RxNorm normalized forms.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 0 | RXCUI | varchar(8) | RxNorm Concept Unique Identifier |
| 1 | LAT | varchar(3) | Language (ENG) |
| 2 | TS | varchar(1) | Term status |
| 3 | LUI | varchar(10) | Lexical unique identifier |
| 4 | STT | varchar(3) | String type |
| 5 | SUI | varchar(10) | String unique identifier |
| 6 | ISPREF | varchar(1) | Atom status (Y=preferred) |
| 7 | RXAUI | varchar(8) | RxNorm Atom Unique Identifier |
| 8 | SAUI | varchar(50) | Source atom identifier |
| 9 | SCUI | varchar(50) | Source concept identifier |
| 10 | SDUI | varchar(50) | Source descriptor identifier |
| 11 | SAB | varchar(20) | Source abbreviation |
| 12 | TTY | varchar(20) | Term type in source |
| 13 | CODE | varchar(50) | Source-specific code |
| 14 | STR | varchar(3000) | String (drug name) |
| 15 | SRL | varchar(10) | Source restriction level |
| 16 | SUPPRESS | varchar(1) | Suppress flag |
| 17 | CVF | varchar(50) | Content view flag |

#### RXNREL.RRF (Relationships)

Contains relationships between concepts.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 0 | RXCUI1 | varchar(8) | First concept identifier |
| 1 | RXAUI1 | varchar(8) | First atom identifier |
| 2 | STYPE1 | varchar(50) | Type of first identifier |
| 3 | REL | varchar(4) | Relationship type |
| 4 | RXCUI2 | varchar(8) | Second concept identifier |
| 5 | RXAUI2 | varchar(8) | Second atom identifier |
| 6 | STYPE2 | varchar(50) | Type of second identifier |
| 7 | RELA | varchar(100) | Relationship attribute |
| 8 | RUI | varchar(10) | Relationship identifier |
| 9 | SRUI | varchar(50) | Source relationship ID |
| 10 | SAB | varchar(20) | Source abbreviation |
| 11 | SL | varchar(1000) | Source of relationship label |
| 12 | DIR | varchar(1) | Directionality flag |
| 13 | RG | varchar(10) | Relationship group |
| 14 | SUPPRESS | varchar(1) | Suppress flag |
| 15 | CVF | varchar(50) | Content view flag |

#### RXNSAT.RRF (Attributes)

Contains attributes for concepts and atoms.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 0 | RXCUI | varchar(8) | Concept identifier |
| 1 | LUI | varchar(10) | Lexical unique identifier |
| 2 | SUI | varchar(10) | String unique identifier |
| 3 | RXAUI | varchar(8) | Atom identifier |
| 4 | STYPE | varchar(50) | Subject type |
| 5 | CODE | varchar(50) | Source-specific code |
| 6 | ATUI | varchar(10) | Attribute identifier |
| 7 | SATUI | varchar(50) | Source attribute identifier |
| 8 | ATN | varchar(100) | Attribute name |
| 9 | SAB | varchar(20) | Source abbreviation |
| 10 | ATV | varchar(4000) | Attribute value |
| 11 | SUPPRESS | varchar(1) | Suppress flag |
| 12 | CVF | varchar(50) | Content view flag |

#### RXNSTY.RRF (Semantic Types)

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 0 | RXCUI | varchar(8) | Concept identifier |
| 1 | TUI | varchar(4) | Semantic type identifier |
| 2 | STN | varchar(100) | Semantic type tree number |
| 3 | STY | varchar(50) | Semantic type name |
| 4 | ATUI | varchar(10) | Attribute identifier |
| 5 | CVF | varchar(50) | Content view flag |

#### RXNDOC.RRF (Documentation)

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 0 | DOCKEY | varchar(50) | Document key |
| 1 | VALUE | varchar(1000) | Value |
| 2 | TYPE | varchar(50) | Type of documentation |
| 3 | EXPL | varchar(1000) | Explanation |

---

### 3.2 Concept Types (TTY)

RxNorm Term Types define the semantic level of drug concepts.

#### Ingredient Level

| TTY | Name | Description | Example |
|-----|------|-------------|---------|
| IN | Ingredient | Single active ingredient | Acetaminophen |
| MIN | Multiple Ingredients | Combination of ingredients | Acetaminophen / Codeine |
| PIN | Precise Ingredient | Specific molecular form | Acetaminophen 4-Nitrophenyl Ether |

#### Clinical Drug Level

| TTY | Name | Description | Example |
|-----|------|-------------|---------|
| SCDC | Semantic Clinical Drug Component | Ingredient + strength | Acetaminophen 325 MG |
| SCDF | Semantic Clinical Drug Form | Ingredient + dose form | Acetaminophen Oral Tablet |
| SCDG | Semantic Clinical Drug Group | Ingredient + dose form group | Acetaminophen Oral Product |
| SCD | Semantic Clinical Drug | Ingredient + strength + form | Acetaminophen 325 MG Oral Tablet |

#### Branded Drug Level

| TTY | Name | Description | Example |
|-----|------|-------------|---------|
| BN | Brand Name | Proprietary name only | Tylenol |
| SBDC | Semantic Branded Drug Component | Brand + ingredient + strength | Tylenol Acetaminophen 325 MG |
| SBDF | Semantic Branded Drug Form | Brand + ingredient + form | Tylenol Acetaminophen Oral Tablet |
| SBDG | Semantic Branded Drug Group | Brand + dose form group | Tylenol Oral Product |
| SBD | Semantic Branded Drug | Brand + strength + form | Tylenol 325 MG Oral Tablet |

#### Pack Level

| TTY | Name | Description | Example |
|-----|------|-------------|---------|
| GPCK | Generic Pack | Multiple drugs in pack | {21 (Ethinyl Estradiol...) / 7 (Inert...)} Pack |
| BPCK | Branded Pack | Branded multi-drug pack | {21 (...) / 7 (...)} Pack [Ortho Tri-Cyclen] |

#### Form Level

| TTY | Name | Description | Example |
|-----|------|-------------|---------|
| DF | Dose Form | Physical form | Oral Tablet |
| DFG | Dose Form Group | Form group by route | Oral Product |

#### Other Types

| TTY | Name | Description |
|-----|------|-------------|
| SY | Synonym | Alternative name |
| TMSY | Tall Man Lettering Synonym | Safety-enhanced spelling |
| ET | Entry Term | Search term |
| PSN | Prescribable Name | Name for prescribing |

---

### 3.3 Relationship Types

#### REL (Relationship) Codes

| REL | Description |
|-----|-------------|
| AQ | Allowed qualifier |
| CHD | Has child |
| DEL | Deleted concept |
| PAR | Has parent |
| QB | Qualified by |
| RB | Has broader relationship |
| RL | Similar or alike |
| RN | Has narrower relationship |
| RO | Has other relationship |
| RQ | Related and possibly synonymous |
| RU | Related (unspecified) |
| SIB | Sibling |
| SY | Synonym |

#### RELA (Relationship Attribute) Codes

| RELA | Description |
|------|-------------|
| consists_of | Contains component |
| constitutes | Part of |
| contained_in | Contained in |
| contains | Contains |
| dose_form_of | Is dose form of |
| form_of | Is form of |
| has_dose_form | Has dose form |
| has_form | Has form |
| has_ingredient | Contains ingredient |
| has_ingredients | Contains multiple ingredients |
| has_part | Has part |
| has_precise_ingredient | Has precise ingredient |
| has_quantified_form | Has quantified form |
| has_tradename | Has brand name |
| ingredient_of | Is ingredient of |
| ingredients_of | Are ingredients of |
| inverse_isa | Inverse is-a |
| isa | Is a |
| part_of | Is part of |
| precise_ingredient_of | Is precise ingredient of |
| quantified_form_of | Is quantified form of |
| tradename_of | Is brand name of |

---

### 3.4 REST API

**Base URL**: `https://rxnav.nlm.nih.gov/REST/`

#### Core Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/rxcui.json?name={name}` | GET | Find RXCUI by drug name |
| `/rxcui/{rxcui}.json` | GET | Get concept properties |
| `/rxcui/{rxcui}/allrelated.json` | GET | Get all related concepts |
| `/rxcui/{rxcui}/related.json?tty={tty}` | GET | Get related concepts by TTY |
| `/rxcui/{rxcui}/related.json?rela={rela}` | GET | Get related concepts by relationship |
| `/rxcui/{rxcui}/ndcs.json` | GET | Get NDCs for concept |
| `/rxcui/{rxcui}/filter.json` | GET | Filter by properties |
| `/drugs.json?name={name}` | GET | Get drug products by name |
| `/approximateTerm.json?term={term}` | GET | Fuzzy match drug name |
| `/spellingsuggestions.json?name={name}` | GET | Get spelling suggestions |
| `/version.json` | GET | Get RxNorm version |
| `/allconcepts.json?tty={tty}` | GET | Get all concepts by TTY |
| `/allstatus.json?status={status}` | GET | Get concepts by status |

#### NDC-Related Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/ndcstatus.json?ndc={ndc}` | GET | Get NDC status |
| `/ndcproperties.json?id={ndc}` | GET | Get NDC properties |
| `/rxcui.json?idtype=NDC&id={ndc}` | GET | Find RXCUI by NDC |
| `/relatedndc.json?ndc={ndc}` | GET | Get related NDCs |

#### Response Format

```json
{
  "rxnormdata": {
    "idGroup": {
      "rxnormId": [
        {
          "rxcui": "161"
        }
      ]
    }
  }
}
```

#### NDC Status Values

| Status | Description |
|--------|-------------|
| ACTIVE | NDC is active in current RxNorm |
| OBSOLETE | NDC was previously active, now inactive |
| ALIEN | NDC exists only in non-RxNorm source |
| UNKNOWN | NDC not recognized |

---

### 3.5 RxClass Integration

RxClass links drug classes to RxNorm drug concepts.

**Base URL**: `https://rxnav.nlm.nih.gov/REST/rxclass/`

#### Classification Systems

| System | Description | Source |
|--------|-------------|--------|
| ATC | Anatomical Therapeutic Chemical | WHO |
| EPC | Established Pharmacologic Class | FDA |
| MoA | Mechanism of Action | MED-RT |
| PE | Physiologic Effect | MED-RT |
| TC | Therapeutic Category | MED-RT |
| VA | VA Drug Class | VA |
| DISEASE | Disease relationship | MED-RT |

#### Endpoints

| Endpoint | Description |
|----------|-------------|
| `/class/byDrugName.json?drugName={name}` | Get classes for drug |
| `/class/byRxcui.json?rxcui={rxcui}` | Get classes by RXCUI |
| `/classMembers.json?classId={id}&relaSource={src}` | Get drugs in class |
| `/classTree.json?classId={id}` | Get class hierarchy |
| `/allClasses.json?classTypes={types}` | List all classes |

#### Example Response

```json
{
  "rxclassMinConceptList": {
    "rxclassMinConcept": [
      {
        "classId": "N0000175503",
        "className": "Proton Pump Inhibitors",
        "classType": "EPC"
      },
      {
        "classId": "N0000175696",
        "className": "Decreased Gastric Acid Secretion",
        "classType": "PE"
      }
    ]
  }
}
```

---

## 4. ClinicalTrials.gov

ClinicalTrials.gov is the U.S. registry of clinical studies.

**Website**: https://clinicaltrials.gov

**API Version**: 2.0 (REST API with OpenAPI 3.0 specification)

### 4.1 Study Data Structure

#### Protocol Section

```json
{
  "protocolSection": {
    "identificationModule": {
      "nctId": "NCT12345678",
      "orgStudyIdInfo": {
        "id": "SPONSOR-STUDY-001"
      },
      "organization": {
        "fullName": "Organization Name",
        "class": "INDUSTRY"
      },
      "briefTitle": "Study Short Title",
      "officialTitle": "Complete Official Study Title",
      "acronym": "STUDY"
    },
    "statusModule": {
      "statusVerifiedDate": "2024-01",
      "overallStatus": "RECRUITING",
      "expandedAccessInfo": {
        "hasExpandedAccess": false
      },
      "startDateStruct": {
        "date": "2024-01-15",
        "type": "ACTUAL"
      },
      "primaryCompletionDateStruct": {
        "date": "2025-12-31",
        "type": "ESTIMATED"
      },
      "completionDateStruct": {
        "date": "2026-06-30",
        "type": "ESTIMATED"
      },
      "studyFirstSubmitDate": "2023-12-01",
      "studyFirstPostDateStruct": {
        "date": "2023-12-15",
        "type": "ACTUAL"
      },
      "lastUpdatePostDateStruct": {
        "date": "2024-06-01",
        "type": "ACTUAL"
      }
    },
    "sponsorCollaboratorsModule": {
      "responsibleParty": {
        "type": "SPONSOR"
      },
      "leadSponsor": {
        "name": "Sponsor Name",
        "class": "INDUSTRY"
      },
      "collaborators": [
        {
          "name": "Collaborator Name",
          "class": "OTHER"
        }
      ]
    },
    "oversightModule": {
      "isFdaRegulatedDrug": true,
      "isFdaRegulatedDevice": false,
      "isUnapprovedDevice": false
    },
    "descriptionModule": {
      "briefSummary": "Brief description of the study...",
      "detailedDescription": "Detailed description..."
    },
    "conditionsModule": {
      "conditions": ["Condition 1", "Condition 2"],
      "keywords": ["keyword1", "keyword2"]
    },
    "designModule": {
      "studyType": "INTERVENTIONAL",
      "phases": ["PHASE3"],
      "designInfo": {
        "allocation": "RANDOMIZED",
        "interventionModel": "PARALLEL",
        "interventionModelDescription": "Description...",
        "primaryPurpose": "TREATMENT",
        "maskingInfo": {
          "masking": "DOUBLE",
          "maskingDescription": "Description...",
          "whoMasked": ["PARTICIPANT", "INVESTIGATOR"]
        }
      },
      "enrollmentInfo": {
        "count": 500,
        "type": "ESTIMATED"
      }
    },
    "armsInterventionsModule": {
      "armGroups": [...],
      "interventions": [...]
    },
    "outcomesModule": {
      "primaryOutcomes": [...],
      "secondaryOutcomes": [...]
    },
    "eligibilityModule": {
      "eligibilityCriteria": "Inclusion/exclusion criteria text...",
      "healthyVolunteers": false,
      "sex": "ALL",
      "minimumAge": "18 Years",
      "maximumAge": "75 Years",
      "stdAges": ["ADULT", "OLDER_ADULT"]
    },
    "contactsLocationsModule": {
      "centralContacts": [...],
      "overallOfficials": [...],
      "locations": [...]
    },
    "referencesModule": {
      "references": [...],
      "seeAlsoLinks": [...]
    }
  }
}
```

#### Results Section

```json
{
  "resultsSection": {
    "participantFlowModule": {
      "preAssignmentDetails": "Details...",
      "recruitmentDetails": "Details...",
      "groups": [...],
      "periods": [...]
    },
    "baselineCharacteristicsModule": {
      "populationDescription": "Description...",
      "groups": [...],
      "denoms": [...],
      "measures": [...]
    },
    "outcomeMeasuresModule": {
      "outcomeMeasures": [...]
    },
    "adverseEventsModule": {
      "frequencyThreshold": "5",
      "timeFrame": "Duration...",
      "description": "Description...",
      "eventGroups": [...],
      "seriousEvents": [...],
      "otherEvents": [...]
    },
    "moreInfoModule": {
      "certainAgreement": {...},
      "pointOfContact": {...}
    }
  }
}
```

---

### 4.2 Intervention Types

| Type | Description |
|------|-------------|
| DRUG | Pharmaceutical intervention |
| DEVICE | Medical device |
| BIOLOGICAL | Biological/vaccine |
| PROCEDURE | Surgical or other procedure |
| RADIATION | Radiation therapy |
| BEHAVIORAL | Behavioral intervention |
| GENETIC | Gene transfer/therapy |
| DIETARY_SUPPLEMENT | Dietary supplement |
| COMBINATION_PRODUCT | Drug/device combination |
| DIAGNOSTIC_TEST | Diagnostic procedure |
| OTHER | Other intervention type |

#### Intervention Object

```json
{
  "interventions": [
    {
      "type": "DRUG",
      "name": "Drug Name",
      "description": "Description of the intervention",
      "armGroupLabels": ["Arm 1", "Arm 2"],
      "otherNames": ["Alternative Name 1"]
    }
  ]
}
```

---

### 4.3 Outcome Measures

#### Primary Outcome

```json
{
  "primaryOutcomes": [
    {
      "measure": "Change in Primary Endpoint",
      "description": "Detailed description of outcome measure",
      "timeFrame": "Baseline to Week 12"
    }
  ]
}
```

#### Secondary Outcome

```json
{
  "secondaryOutcomes": [
    {
      "measure": "Secondary Endpoint Name",
      "description": "Description",
      "timeFrame": "Week 12"
    }
  ]
}
```

#### Outcome Results (Results Section)

```json
{
  "outcomeMeasures": [
    {
      "type": "PRIMARY",
      "title": "Outcome Title",
      "description": "Description...",
      "populationDescription": "Population...",
      "reportingStatus": "POSTED",
      "paramType": "MEAN",
      "dispersionType": "STANDARD_DEVIATION",
      "unitOfMeasure": "units",
      "timeFrame": "Time frame",
      "groups": [...],
      "denoms": [...],
      "classes": [...]
    }
  ]
}
```

---

### 4.4 API v2 Endpoints

**Base URL**: `https://clinicaltrials.gov/api/v2/`

#### Study Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/studies` | GET | Search studies |
| `/studies/{nctId}` | GET | Get single study |
| `/studies/metadata` | GET | Get field metadata |
| `/stats/size` | GET | Get database statistics |

#### Query Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `query.cond` | string | Condition/disease search |
| `query.term` | string | General search term |
| `query.intr` | string | Intervention search |
| `query.titles` | string | Title search |
| `query.outc` | string | Outcome measure search |
| `query.spons` | string | Sponsor search |
| `query.lead` | string | Lead sponsor search |
| `query.id` | string | Study ID search |
| `query.patient` | string | Patient-level search |
| `filter.overallStatus` | array | Status filter |
| `filter.geo` | string | Geographic filter |
| `filter.ids` | array | NCT ID filter |
| `filter.advanced` | string | Advanced filter expression |
| `aggFilters` | string | Aggregation filters |
| `fields` | array | Fields to return |
| `sort` | array | Sort order |
| `pageSize` | integer | Results per page (max 1000) |
| `pageToken` | string | Pagination token |
| `format` | string | Response format (json, csv) |

#### Example Queries

```bash
# Search for diabetes drug trials
https://clinicaltrials.gov/api/v2/studies?query.cond=diabetes&query.intr=DRUG&filter.overallStatus=RECRUITING

# Get specific study
https://clinicaltrials.gov/api/v2/studies/NCT12345678

# Search with field selection
https://clinicaltrials.gov/api/v2/studies?query.term=cancer&fields=NCTId,BriefTitle,OverallStatus&pageSize=50
```

#### Overall Status Values

| Status | Description |
|--------|-------------|
| NOT_YET_RECRUITING | Has not started recruiting |
| RECRUITING | Currently recruiting |
| ENROLLING_BY_INVITATION | Enrolling by invitation only |
| ACTIVE_NOT_RECRUITING | Active, not recruiting |
| SUSPENDED | Temporarily suspended |
| TERMINATED | Stopped early |
| COMPLETED | Study completed |
| WITHDRAWN | Withdrawn before enrollment |
| UNKNOWN | Status unknown |

#### Study Phase Values

| Phase | Description |
|-------|-------------|
| EARLY_PHASE1 | Early Phase 1 |
| PHASE1 | Phase 1 |
| PHASE2 | Phase 2 |
| PHASE3 | Phase 3 |
| PHASE4 | Phase 4 (post-marketing) |
| NA | Not Applicable |

---

## 5. FDA Orange Book

The Orange Book lists approved drug products with therapeutic equivalence evaluations.

**Website**: https://www.fda.gov/drugs/drug-approvals-and-databases/orange-book-data-files

### 5.1 Data File Structure

Orange Book data is provided as tilde-delimited (~) ASCII text files.

#### Products.txt

| Field | Type | Description |
|-------|------|-------------|
| Ingredient | string | Active ingredient(s), alphabetical, semicolon-separated |
| DF;Route | string | Dosage form; Route of administration |
| Trade_Name | string | Proprietary name |
| Applicant | string | Applicant short name (max 20 chars) |
| Strength | string | Product strength |
| Appl_Type | char | N=NDA, A=ANDA |
| Appl_No | string | Application number (6 digits) |
| Product_No | string | Product number (3 digits) |
| TE_Code | string | Therapeutic equivalence code |
| Approval_Date | string | Approval date (Mmm DD, YYYY) |
| RLD | string | Y=Reference Listed Drug |
| RS | string | Y=Reference Standard |
| Type | string | RX, OTC, or DISCN (discontinued) |
| Applicant_Full_Name | string | Full company name |

#### Patent.txt

| Field | Type | Description |
|-------|------|-------------|
| Appl_Type | char | N=NDA, A=ANDA |
| Appl_No | string | Application number (6 digits) |
| Product_No | string | Product number (3 digits) |
| Patent_No | string | Patent number (7-11 digits) |
| Patent_Expire_Date_Text | string | Expiration date (Mmm DD, YYYY) |
| Drug_Substance_Flag | char | Y=Claims drug substance |
| Drug_Product_Flag | char | Y=Claims drug product |
| Patent_Use_Code | string | Use code (if method-of-use patent) |
| Delist_Flag | char | Y=Delist requested |
| Submission_Date | string | Patent submission date |

#### Exclusivity.txt

| Field | Type | Description |
|-------|------|-------------|
| Appl_Type | char | N=NDA, A=ANDA |
| Appl_No | string | Application number (6 digits) |
| Product_No | string | Product number (3 digits) |
| Exclusivity_Code | string | Exclusivity type code |
| Exclusivity_Date | string | Expiration date (Mmm DD, YYYY) |

---

### 5.2 Therapeutic Equivalence Codes

#### A-Rated (Therapeutically Equivalent - Substitutable)

| Code | Description |
|------|-------------|
| AA | Conventional dosage forms not presenting bioequivalence problems |
| AB | Products meeting bioequivalence requirements |
| AB1, AB2, AB3... | Bioequivalent within same group but not between groups |
| AN | Aerosol or nasal solutions/sprays |
| AO | Injectable oil solutions |
| AP | Injectable aqueous solutions |
| AT | Topical products |

#### B-Rated (Not Therapeutically Equivalent)

| Code | Description |
|------|-------------|
| BC | Extended-release dosage forms (oral) |
| BD | Active ingredients with documented bioequivalence problems |
| BE | Delayed-release oral products |
| BN | Aerosol or nasal sprays with problems |
| BP | Active ingredients with potential bioequivalence problems |
| BR | Suppositories or enemas |
| BS | Drug standard deficiencies |
| BT | Topical products with bioequivalence issues |
| BX | Insufficient data for determination |
| B* | Requires further FDA investigation |

---

### 5.3 Exclusivity Codes

#### Primary Exclusivity Types

| Code | Duration | Description |
|------|----------|-------------|
| NCE | 5 years | New Chemical Entity |
| ODE | 7 years | Orphan Drug Exclusivity |
| NCI | 3 years | New Clinical Investigation |
| PC | 180 days | Paragraph IV Patent Challenge (Generic) |
| PED | 6 months | Pediatric Exclusivity (added to existing) |
| GAIN | 5 years | Qualified Infectious Disease Product (added to existing) |

#### Exclusivity Code Variants

| Code | Description |
|------|-------------|
| I-### | Product-specific exclusivity with number identifier |
| D-### | Drug substance exclusivity with number identifier |
| M-### | Miscellaneous exclusivity with number identifier |
| ODE-### | Orphan Drug with indication number |
| NC | New combination |
| NC-### | New combination with identifier |
| NP | New dosage form, formulation, or delivery system |
| NP-### | New product with identifier |

---

### 5.4 Patent Data Format

#### Patent Types

| Flag | Description |
|------|-------------|
| Drug Substance Flag = Y | Patent claims the active ingredient |
| Drug Product Flag = Y | Patent claims the formulation/product |
| Patent Use Code present | Method-of-use patent (claims specific indication) |

#### Patent Use Codes

Use codes identify specific approved uses covered by method-of-use patents.

| Code Format | Description |
|-------------|-------------|
| U-### | Use code identifier (e.g., U-001) |

#### Example Record

```
Products.txt:
ACETAMINOPHEN~TABLET;ORAL~TYLENOL~MCNEIL~325MG~N~016004~001~AB~Apr 1, 1976~Y~Y~RX~McNeil Consumer Healthcare

Patent.txt:
N~016004~001~1234567~Dec 31, 2025~Y~N~~N~Jan 15, 2020

Exclusivity.txt:
N~020560~001~NCE~Mar 15, 2028
N~020560~001~PED~Sep 15, 2028
```

---

## 6. ID System Reference

### Drug Identifiers

| ID Type | Format | Description | Source |
|---------|--------|-------------|--------|
| NDC | 10 digits (4-4-2, 5-3-2, 5-4-1) | National Drug Code | FDA |
| NDA | NXXXXXX (6 digits) | New Drug Application | FDA |
| ANDA | AXXXXXX (6 digits) | Abbreviated NDA (generic) | FDA |
| BLA | XXXXXX (6 digits) | Biologics License Application | FDA |
| RXCUI | 1-7 digits | RxNorm Concept Unique ID | NLM |
| RXAUI | 8 digits | RxNorm Atom Unique ID | NLM |
| UNII | 10 alphanumeric | Unique Ingredient Identifier | FDA |
| SPL ID | UUID | SPL Document ID | NLM |
| SPL Set ID | UUID | SPL Set ID (groups versions) | NLM |
| NCT ID | NCT + 8 digits | ClinicalTrials.gov ID | NLM |

### Clinical Coding Systems

| System | Description | Use |
|--------|-------------|-----|
| MedDRA | Medical Dictionary for Regulatory Activities | Adverse events |
| LOINC | Logical Observation Identifiers Names and Codes | Lab/document sections |
| SNOMED CT | Systematized Nomenclature of Medicine | Clinical terms |
| ICD-10 | International Classification of Diseases | Diagnoses |
| ATC | Anatomical Therapeutic Chemical | Drug classification |

### Organization Identifiers

| ID Type | Description |
|---------|-------------|
| DUNS | Data Universal Numbering System (FDA registration) |
| Labeler Code | 4-5 digit NDC labeler segment |
| FDA Establishment ID | Facility registration number |

---

## Sources and References

### OpenFDA
- [Drug Adverse Event API](https://open.fda.gov/apis/drug/event/)
- [Drug Label API](https://open.fda.gov/apis/drug/label/)
- [Drug NDC API](https://open.fda.gov/apis/drug/ndc/)
- [Drug Enforcement API](https://open.fda.gov/apis/drug/enforcement/)
- [OpenFDA Query Syntax](https://open.fda.gov/apis/query-syntax/)
- [OpenFDA Advanced Syntax](https://open.fda.gov/apis/advanced-syntax/)
- [OpenFDA Fields](https://open.fda.gov/apis/openfda-fields/)
- [GitHub: FDA/openfda](https://github.com/FDA/openfda)

### DailyMed
- [DailyMed SPL Resources](https://dailymed.nlm.nih.gov/dailymed/spl-resources.cfm)
- [FDA SPL Resources](https://www.fda.gov/industry/fda-data-standards-advisory-board/structured-product-labeling-resources)
- [Section Headings (LOINC)](https://www.fda.gov/industry/structured-product-labeling-resources/section-headings-loinc)

### RxNorm
- [RxNorm Overview](https://www.nlm.nih.gov/research/umls/rxnorm/overview.html)
- [RxNorm Technical Documentation](https://www.nlm.nih.gov/research/umls/rxnorm/docs/techdoc.html)
- [RxNorm Files](https://www.nlm.nih.gov/research/umls/rxnorm/docs/rxnormfiles.html)
- [RxNorm API Documentation](https://lhncbc.nlm.nih.gov/RxNav/APIs/RxNormAPIs.html)
- [RxClass Overview](https://lhncbc.nlm.nih.gov/RxNav/applications/RxClassIntro.html)

### ClinicalTrials.gov
- [API Documentation](https://clinicaltrials.gov/data-api/api)
- [Study Data Structure](https://clinicaltrials.gov/data-api/about-api/study-data-structure)
- [API Version 2.0 Release](https://www.nlm.nih.gov/pubs/techbull/ma24/ma24_clinicaltrials_api.html)

### FDA Orange Book
- [Orange Book Data Files](https://www.fda.gov/drugs/drug-approvals-and-databases/orange-book-data-files)
- [Orange Book Preface](https://www.fda.gov/drugs/development-approval-process-drugs/orange-book-preface)
- [Orange Book FAQ](https://www.fda.gov/drugs/drug-approvals-and-databases/frequently-asked-questions-orange-book)
- [Patents and Exclusivity FAQ](https://www.fda.gov/drugs/development-approval-process-drugs/frequently-asked-questions-patents-and-exclusivity)
- [NBER Orange Book Dataset Guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC10731339/)

---

*Document Version: 1.0*
*Last Updated: January 2025*
