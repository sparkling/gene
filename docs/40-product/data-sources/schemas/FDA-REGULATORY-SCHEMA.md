# FDA Regulatory Data Schema

**Document ID:** FDA-REGULATORY-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**License:** Public Domain (US Government Work)

---

## TL;DR

This document covers all major US FDA and NLM regulatory drug databases: OpenFDA (4 drug APIs), DailyMed (SPL labels with 133+ LOINC sections), RxNorm (normalized drug nomenclature with 5 RRF files), ClinicalTrials.gov (API v2), and FDA Orange Book (therapeutic equivalence). All data is Public Domain as US Government work. Key cross-references: NDC, RXCUI, UNII, NDA/ANDA/BLA numbers, and NCT IDs.

---

## Database Overview

| Database | Maintainer | Records | Update Frequency | Format |
|----------|------------|---------|------------------|--------|
| **OpenFDA Drug Event** | FDA | 25M+ reports | Quarterly | JSON API |
| **OpenFDA Drug Label** | FDA | 150K+ labels | Weekly | JSON API |
| **OpenFDA Drug NDC** | FDA | 350K+ products | Daily | JSON API |
| **OpenFDA Enforcement** | FDA | 50K+ recalls | Weekly | JSON API |
| **DailyMed** | NLM | 150K+ SPL docs | Daily | XML/REST |
| **RxNorm** | NLM | 1.2M+ concepts | Monthly | RRF/REST |
| **ClinicalTrials.gov** | NLM | 500K+ studies | Daily | JSON API |
| **Orange Book** | FDA | 35K+ products | Monthly | Delimited |

---

## 1. OpenFDA APIs

**Base URL:** `https://api.fda.gov`
**Authentication:** Optional API key for higher rate limits
**Rate Limits:** 240 requests/minute without key, 120,000/day with key

### 1.1 Drug Event API (FAERS)

**Endpoint:** `https://api.fda.gov/drug/event.json`
**Data Source:** FDA Adverse Event Reporting System (2004-present)
**Standard:** ICH E2b/M2 version 2.1
**Coding System:** MedDRA (Medical Dictionary for Regulatory Activities)

#### Top-Level Report Fields

| Field | Type | Description |
|-------|------|-------------|
| `safetyreportid` | string | Unique 8-digit identifier (7 digits + checksum) |
| `safetyreportversion` | string | Version number of the report |
| `receivedate` | string | Date FDA received report (YYYYMMDD) |
| `receivedateformat` | string | Format code for date (102 = CCYYMMDD) |
| `receiptdate` | string | Date most recent information received |
| `transmissiondate` | string | Date report transmitted |
| `reporttype` | string | 1=Spontaneous, 2=Literature, 3=Study |
| `serious` | string | 1=Serious, 2=Not serious |
| `seriousnessdeath` | string | 1=Patient died |
| `seriousnesshospitalization` | string | 1=Patient hospitalized |
| `seriousnesslifethreatening` | string | 1=Life-threatening condition |
| `seriousnessdisabling` | string | 1=Resulted in disability |
| `seriousnesscongenitalanomali` | string | 1=Congenital anomaly |
| `seriousnessother` | string | 1=Other serious outcome |
| `duplicate` | string | 1=Duplicate report |
| `companynumb` | string | Manufacturer's internal case number |
| `authoritynumb` | string | Regulatory authority case number |
| `occurcountry` | string | Country where event occurred (ISO 3166-1 alpha-2) |
| `primarysourcecountry` | string | Country of primary source |
| `fulfillexpediteriteria` | string | Expedited reporting criteria met |

#### Patient Object (`patient`)

| Field | Type | Description |
|-------|------|-------------|
| `patient.patientonsetage` | string | Patient age at event onset |
| `patient.patientonsetageunit` | string | 800=Decade, 801=Year, 802=Month, 803=Week, 804=Day, 805=Hour |
| `patient.patientsex` | string | 0=Unknown, 1=Male, 2=Female |
| `patient.patientweight` | string | Weight in kg |
| `patient.patientdeathdate` | string | Date of death (YYYYMMDD) |
| `patient.patientdeathdateformat` | string | Format code for death date |
| `patient.reaction` | array | Array of reaction objects |
| `patient.drug` | array | Array of drug objects |

#### Reaction Object (`patient.reaction[]`)

| Field | Type | Description |
|-------|------|-------------|
| `reactionmeddrapt` | string | MedDRA preferred term for reaction |
| `reactionmeddraversionpt` | string | MedDRA version used |
| `reactionoutcome` | string | 1=Recovered, 2=Recovering, 3=Not recovered, 4=Recovered with sequelae, 5=Fatal, 6=Unknown |

#### Drug Object (`patient.drug[]`)

| Field | Type | Description |
|-------|------|-------------|
| `drugcharacterization` | string | 1=Suspect, 2=Concomitant, 3=Interacting |
| `medicinalproduct` | string | Drug name as reported |
| `drugbatchnumb` | string | Lot/batch number |
| `drugauthorizationnumb` | string | NDA/ANDA number |
| `drugstructuredosagenumb` | string | Dose amount |
| `drugstructuredosageunit` | string | Dose unit code |
| `drugseparatedosagenumb` | string | Number of separate dosages |
| `drugintervaldosageunitnumb` | string | Interval between doses |
| `drugintervaldosagedefinition` | string | Interval unit |
| `drugcumulativedosagenumb` | string | Cumulative dose amount |
| `drugcumulativedosageunit` | string | Cumulative dose unit |
| `drugdosagetext` | string | Free-text dosage information |
| `drugdosageform` | string | Dosage form (tablet, capsule, etc.) |
| `drugadministrationroute` | string | Route of administration code |
| `drugindication` | string | Indication for use |
| `drugstartdate` | string | Date drug therapy started |
| `drugenddate` | string | Date drug therapy ended |
| `drugtreatmentduration` | string | Duration of treatment |
| `drugtreatmentdurationunit` | string | Duration unit |
| `actiondrug` | string | 1=Withdrawn, 2=Dose reduced, 3=Dose increased, 4=Unchanged, 5=Unknown, 6=Not applicable |
| `drugrecurreadministration` | string | 1=Yes, 2=No, 3=Unknown (rechallenge) |
| `activesubstance` | object | Active ingredient information |
| `openfda` | object | OpenFDA harmonized fields |

#### Primary Source Object (`primarysource`)

| Field | Type | Description |
|-------|------|-------------|
| `qualification` | string | 1=Physician, 2=Pharmacist, 3=Other health professional, 4=Lawyer, 5=Consumer |
| `reportercountry` | string | Reporter's country code |

#### Sample JSON Response

```json
{
  "meta": {
    "disclaimer": "...",
    "terms": "...",
    "results": {"skip": 0, "limit": 1, "total": 25432109}
  },
  "results": [{
    "safetyreportid": "10003214",
    "safetyreportversion": "1",
    "receivedate": "20240115",
    "receivedateformat": "102",
    "serious": "1",
    "seriousnesshospitalization": "1",
    "occurcountry": "US",
    "patient": {
      "patientsex": "2",
      "patientonsetage": "45",
      "patientonsetageunit": "801",
      "reaction": [{
        "reactionmeddrapt": "Nausea",
        "reactionmeddraversionpt": "25.1",
        "reactionoutcome": "1"
      }],
      "drug": [{
        "drugcharacterization": "1",
        "medicinalproduct": "METFORMIN",
        "drugindication": "Type 2 diabetes mellitus",
        "drugdosagetext": "500 MG, ORAL, TWICE DAILY",
        "openfda": {
          "brand_name": ["GLUCOPHAGE"],
          "generic_name": ["METFORMIN HYDROCHLORIDE"],
          "rxcui": ["861004"],
          "application_number": ["NDA020357"]
        }
      }]
    },
    "primarysource": {
      "qualification": "1",
      "reportercountry": "US"
    }
  }]
}
```

---

### 1.2 Drug Label API

**Endpoint:** `https://api.fda.gov/drug/label.json`
**Data Source:** Structured Product Labeling (SPL) submissions

#### Metadata Fields

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Unique SPL document identifier (UUID) |
| `set_id` | string | Set ID grouping all versions of a label (UUID) |
| `version` | string | Label version number |
| `effective_time` | string | Date label became effective (YYYYMMDD) |

#### Clinical Information Sections

| Field | Type | Description |
|-------|------|-------------|
| `indications_and_usage` | array | Approved therapeutic uses |
| `dosage_and_administration` | array | Dosing instructions |
| `dosage_forms_and_strengths` | array | Available formulations |
| `contraindications` | array | Conditions where drug should not be used |
| `warnings_and_cautions` | array | Safety warnings (PLR format) |
| `warnings` | array | Safety warnings (older format) |
| `precautions` | array | Precautionary information |
| `adverse_reactions` | array | Reported side effects |
| `drug_interactions` | array | Drug-drug interactions |
| `use_in_specific_populations` | array | Pregnancy, nursing, pediatric, geriatric use |
| `overdosage` | array | Overdose management information |
| `boxed_warning` | array | Black box warnings |

#### Pharmacology Sections

| Field | Type | Description |
|-------|------|-------------|
| `clinical_pharmacology` | array | Mechanism and pharmacokinetics overview |
| `mechanism_of_action` | array | How the drug works |
| `pharmacodynamics` | array | Drug effects on body |
| `pharmacokinetics` | array | ADME (absorption, distribution, metabolism, excretion) |

#### Product Information Sections

| Field | Type | Description |
|-------|------|-------------|
| `description` | array | Drug description and formulation |
| `how_supplied` | array | Available packages and NDC codes |
| `storage_and_handling` | array | Storage requirements |
| `active_ingredient` | array | Active components |
| `inactive_ingredient` | array | Inactive components |
| `clinical_studies` | array | Clinical trial summaries |
| `nonclinical_toxicology` | array | Animal study data |
| `references` | array | Scientific references |

#### OTC-Specific Sections

| Field | Type | Description |
|-------|------|-------------|
| `purpose` | array | Therapeutic purpose |
| `do_not_use` | array | Contraindications for OTC |
| `stop_use` | array | When to stop taking |
| `ask_doctor` | array | When to consult physician |
| `ask_doctor_or_pharmacist` | array | Consultation guidance |
| `keep_out_of_reach_of_children` | array | Child safety warning |
| `directions` | array | Usage directions |
| `other_safety_information` | array | Additional safety info |
| `questions` | array | Contact information |

#### OpenFDA Harmonized Fields

| Field | Type | Description |
|-------|------|-------------|
| `openfda.application_number` | array | NDA/ANDA/BLA number |
| `openfda.brand_name` | array | Proprietary name(s) |
| `openfda.generic_name` | array | Non-proprietary name(s) |
| `openfda.manufacturer_name` | array | Manufacturer/labeler name(s) |
| `openfda.product_ndc` | array | Product NDC codes |
| `openfda.product_type` | array | Product type classification |
| `openfda.route` | array | Route(s) of administration |
| `openfda.substance_name` | array | Active substance name(s) |
| `openfda.rxcui` | array | RxNorm concept identifiers |
| `openfda.spl_id` | array | SPL document ID |
| `openfda.spl_set_id` | array | SPL set ID |
| `openfda.pharm_class_epc` | array | Established pharmacologic class |
| `openfda.pharm_class_moa` | array | Mechanism of action class |
| `openfda.pharm_class_pe` | array | Physiologic effect class |
| `openfda.pharm_class_cs` | array | Chemical structure class |
| `openfda.unii` | array | UNII codes |

#### Sample JSON Response

```json
{
  "results": [{
    "id": "12345678-1234-1234-1234-123456789abc",
    "set_id": "abcdef12-3456-7890-abcd-ef1234567890",
    "version": "3",
    "effective_time": "20250615",
    "indications_and_usage": [
      "METFORMIN HYDROCHLORIDE TABLETS are indicated as an adjunct to diet and exercise to improve glycemic control in adults with type 2 diabetes mellitus."
    ],
    "dosage_and_administration": [
      "Starting dose: 500 mg twice daily or 850 mg once daily with meals. Maximum dose: 2550 mg daily in divided doses."
    ],
    "contraindications": [
      "Severe renal impairment (eGFR below 30 mL/min/1.73 m2). Known hypersensitivity to metformin."
    ],
    "boxed_warning": [
      "LACTIC ACIDOSIS: Metformin can cause lactic acidosis, a rare but serious complication..."
    ],
    "adverse_reactions": [
      "Most common adverse reactions (incidence >= 5%): diarrhea, nausea, vomiting, flatulence, asthenia, headache."
    ],
    "openfda": {
      "application_number": ["NDA020357"],
      "brand_name": ["GLUCOPHAGE"],
      "generic_name": ["METFORMIN HYDROCHLORIDE"],
      "manufacturer_name": ["BRISTOL-MYERS SQUIBB"],
      "product_ndc": ["0087-6060"],
      "product_type": ["HUMAN PRESCRIPTION DRUG"],
      "route": ["ORAL"],
      "rxcui": ["861004"],
      "unii": ["9100L32L2N"],
      "pharm_class_epc": ["Biguanide [EPC]"]
    }
  }]
}
```

---

### 1.3 Drug NDC API

**Endpoint:** `https://api.fda.gov/drug/ndc.json`
**Data Source:** NDC Directory (updated daily)

#### NDC Structure

The National Drug Code is a 10-digit, 3-segment identifier:
- **Segment 1 (Labeler):** 4-5 digits - manufacturer/distributor
- **Segment 2 (Product):** 3-4 digits - drug formulation
- **Segment 3 (Package):** 1-2 digits - package size

**Configurations:** 4-4-2, 5-3-2, or 5-4-1

#### Product Information Fields

| Field | Type | Description |
|-------|------|-------------|
| `product_id` | string | Unique product identifier |
| `product_ndc` | string | Product NDC (labeler-product segments) |
| `spl_id` | string | SPL document identifier |
| `product_type` | string | HUMAN PRESCRIPTION DRUG, HUMAN OTC DRUG, etc. |
| `finished` | string | "true"=finished drug, "false"=active ingredient |
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

#### Active Ingredients Object (`active_ingredients[]`)

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Ingredient name |
| `strength` | string | Amount with units (e.g., "500 mg/1") |

#### Packaging Object (`packaging[]`)

| Field | Type | Description |
|-------|------|-------------|
| `package_ndc` | string | Full NDC including package segment |
| `description` | string | Package description |
| `marketing_start_date` | string | Package marketing start date |
| `marketing_end_date` | string | Package marketing end date |
| `sample` | boolean | True if sample package |

#### Pharmacologic Class (`pharm_class[]`)

Pharmacologic classes with type suffix: `[EPC]`, `[MoA]`, `[PE]`, `[CS]`

#### Sample JSON Response

```json
{
  "results": [{
    "product_id": "12345-6789_12345678-1234",
    "product_ndc": "12345-6789",
    "product_type": "HUMAN PRESCRIPTION DRUG",
    "brand_name": "METFORMIN HYDROCHLORIDE",
    "generic_name": "METFORMIN HYDROCHLORIDE",
    "dosage_form": "TABLET",
    "route": ["ORAL"],
    "marketing_category": "ANDA",
    "application_number": "ANDA076543",
    "labeler_name": "GENERIC PHARMA INC",
    "active_ingredients": [{
      "name": "METFORMIN HYDROCHLORIDE",
      "strength": "500 mg/1"
    }],
    "packaging": [{
      "package_ndc": "12345-6789-01",
      "description": "100 TABLET in 1 BOTTLE",
      "marketing_start_date": "20100115"
    }],
    "pharm_class": [
      "Biguanide [EPC]",
      "Biguanides [Chemical/Ingredient]"
    ],
    "openfda": {
      "rxcui": ["861004"],
      "unii": ["9100L32L2N"]
    }
  }]
}
```

---

### 1.4 Drug Enforcement API

**Endpoint:** `https://api.fda.gov/drug/enforcement.json`
**Data Source:** FDA Recall Enterprise System (RES) (2004-present, weekly updates)

#### Recall Information Fields

| Field | Type | Description |
|-------|------|-------------|
| `recall_number` | string | Unique recall identifier (e.g., D-1234-2024) |
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
| `initial_firm_notification` | string | Letter, Email, Press Release, etc. |
| `distribution_pattern` | string | Geographic distribution description |

#### Product Information Fields

| Field | Type | Description |
|-------|------|-------------|
| `product_type` | string | Drugs, Devices, Food, etc. |
| `product_description` | string | Description of recalled product |
| `product_quantity` | string | Amount of product recalled |
| `code_info` | string | Lot numbers, expiration dates |
| `more_code_info` | string | Additional code information |
| `product_code` | string | Product code |

#### Recall Dates

| Field | Type | Description |
|-------|------|-------------|
| `recall_initiation_date` | string | When recall began (YYYYMMDD) |
| `center_classification_date` | string | When FDA classified recall |
| `report_date` | string | Date reported (YYYYMMDD) |
| `termination_date` | string | When recall terminated |

#### Reason Field

| Field | Type | Description |
|-------|------|-------------|
| `reason_for_recall` | string | Explanation of why recall occurred |

#### Classification Definitions

| Class | Definition | Example Issues |
|-------|------------|----------------|
| **Class I** | Dangerous/defective products that could cause serious health problems or death | Contamination, mislabeling of potent drugs |
| **Class II** | Might cause temporary health problem or slight threat of serious nature | Subpotent, dissolution failures |
| **Class III** | Unlikely to cause adverse health reaction but violate FDA regulations | Label errors, minor GMP violations |

#### Sample JSON Response

```json
{
  "results": [{
    "recall_number": "D-1234-2024",
    "event_id": "87654",
    "status": "Ongoing",
    "classification": "Class II",
    "recalling_firm": "EXAMPLE PHARMACEUTICALS",
    "city": "New York",
    "state": "NY",
    "country": "United States",
    "voluntary_mandated": "Voluntary: Firm Initiated",
    "product_type": "Drugs",
    "product_description": "METFORMIN HCl ER Tablets, 500mg, 90 count bottle",
    "product_quantity": "50,000 bottles",
    "code_info": "Lot numbers: MET2024A, MET2024B. Expiration: 12/2025, 03/2026",
    "distribution_pattern": "Nationwide distribution in the United States",
    "reason_for_recall": "Failed dissolution specifications at 6-month stability",
    "recall_initiation_date": "20240601",
    "report_date": "20240615",
    "openfda": {
      "brand_name": ["METFORMIN HCL ER"],
      "generic_name": ["METFORMIN HYDROCHLORIDE"],
      "application_number": ["ANDA076543"]
    }
  }]
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
| `field:[min+TO+max]` | Range (inclusive) | `receivedate:[20230101+TO+20231231]` |
| `field:term*` | Wildcard | `brand_name:metfor*` |
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

#### Example Queries

```bash
# Find metformin adverse events with hospitalization
https://api.fda.gov/drug/event.json?search=patient.drug.openfda.generic_name:metformin+AND+seriousnesshospitalization:1&limit=10

# Count reactions for a drug
https://api.fda.gov/drug/event.json?search=patient.drug.openfda.brand_name:glucophage&count=patient.reaction.reactionmeddrapt.exact

# Search drug labels for diabetes
https://api.fda.gov/drug/label.json?search=indications_and_usage:diabetes&limit=25

# Find Class I recalls in 2024
https://api.fda.gov/drug/enforcement.json?search=classification:"Class+I"+AND+report_date:[20240101+TO+20241231]&limit=100

# Find NDC products by labeler
https://api.fda.gov/drug/ndc.json?search=labeler_name:"PFIZER"&limit=100
```

---

### 1.6 OpenFDA Harmonized Fields Reference

These `openfda` fields are added to records across all endpoints:

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
| `unii` | array | Unique Ingredient Identifier |

---

## 2. DailyMed

**Website:** https://dailymed.nlm.nih.gov
**Maintainer:** National Library of Medicine (NLM)
**Format:** SPL XML (HL7 v3 standard)

### 2.1 SPL XML Schema

Structured Product Labeling (SPL) is an HL7 v3 standard for drug labeling.

#### Document Structure

```xml
<?xml version="1.0" encoding="UTF-8"?>
<document xmlns="urn:hl7-org:v3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <!-- Document header -->
  <id root="document-uuid"/>
  <code code="34391-3" codeSystem="2.16.840.1.113883.6.1"
        displayName="HUMAN PRESCRIPTION DRUG LABEL"/>
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
        <activeMoiety>
          <activeMoiety>
            <name>Active Moiety Name</name>
          </activeMoiety>
        </activeMoiety>
      </ingredientSubstance>
    </ingredient>
  </manufacturedMedicine>
</manufacturedProduct>
```

---

### 2.2 Key Sections

#### Prescription Drug Label Sections (PLR Format)

| Section | LOINC | Description |
|---------|-------|-------------|
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

### 2.3 LOINC Section Codes (133+ Codes)

#### Core Label Sections

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

#### PLR (Physician Labeling Rule) Sections

| LOINC Code | Section Name |
|------------|--------------|
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

#### Package Label and Display

| LOINC Code | Section Name |
|------------|--------------|
| 51945-4 | PACKAGE LABEL.PRINCIPAL DISPLAY PANEL |
| 53404-0 | ROUTE OF ADMINISTRATION |
| 55105-1 | OTC - ACTIVE INGREDIENT |
| 55106-9 | OTC - PURPOSE |
| 34071-1 | WARNINGS |
| 50565-1 | OTC - ASK DOCTOR |
| 50566-9 | OTC - ASK DOCTOR/PHARMACIST |
| 50567-7 | OTC - DO NOT USE |
| 50568-5 | OTC - KEEP OUT OF REACH OF CHILDREN |
| 50569-3 | OTC - STOP USE |
| 34068-7 | OTC - DIRECTIONS |
| 53413-1 | OTC - QUESTIONS |

---

### 2.4 FTP Download Structure

**FTP URL:** `ftp://public.nlm.nih.gov/nlmdata/.dailymed/`

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

**Base URL:** `https://dailymed.nlm.nih.gov/dailymed/services/v2/`

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/spls.json` | GET | List all SPLs |
| `/spls/{setId}.xml` | GET | Get SPL by set ID |
| `/drugnames.json` | GET | Search drug names |
| `/drugclasses.json` | GET | Get pharmacologic classes |
| `/ndcs.json` | GET | Get NDC information |
| `/rxcuis.json` | GET | Get RxCUI mappings |

---

## 3. RxNorm

**Website:** https://www.nlm.nih.gov/research/umls/rxnorm/
**Maintainer:** National Library of Medicine (NLM)
**Format:** RRF (Rich Release Format) - pipe-delimited UTF-8

### 3.1 RRF File Format

RxNorm provides 5 core RRF files:

#### RXNCONSO.RRF (Concepts and Atoms)

All drug names from source vocabularies and RxNorm normalized forms.

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

**Sample Record:**
```
861004|ENG|P|L12345678|PF|S12345678|Y|12345678|123456|861004||RXNORM|SCD|861004|metformin hydrochloride 500 MG Oral Tablet|0|N||
```

#### RXNREL.RRF (Relationships)

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
| IN | Ingredient | Single active ingredient | Metformin |
| MIN | Multiple Ingredients | Combination of ingredients | Metformin / Sitagliptin |
| PIN | Precise Ingredient | Specific molecular form | Metformin Hydrochloride |

#### Clinical Drug Level

| TTY | Name | Description | Example |
|-----|------|-------------|---------|
| SCDC | Semantic Clinical Drug Component | Ingredient + strength | Metformin 500 MG |
| SCDF | Semantic Clinical Drug Form | Ingredient + dose form | Metformin Oral Tablet |
| SCDG | Semantic Clinical Drug Group | Ingredient + dose form group | Metformin Oral Product |
| SCD | Semantic Clinical Drug | Ingredient + strength + form | Metformin 500 MG Oral Tablet |

#### Branded Drug Level

| TTY | Name | Description | Example |
|-----|------|-------------|---------|
| BN | Brand Name | Proprietary name only | Glucophage |
| SBDC | Semantic Branded Drug Component | Brand + ingredient + strength | Glucophage Metformin 500 MG |
| SBDF | Semantic Branded Drug Form | Brand + ingredient + form | Glucophage Metformin Oral Tablet |
| SBDG | Semantic Branded Drug Group | Brand + dose form group | Glucophage Oral Product |
| SBD | Semantic Branded Drug | Brand + strength + form | Glucophage 500 MG Oral Tablet |

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

**Base URL:** `https://rxnav.nlm.nih.gov/REST/`

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

#### NDC-Related Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/ndcstatus.json?ndc={ndc}` | GET | Get NDC status |
| `/ndcproperties.json?id={ndc}` | GET | Get NDC properties |
| `/rxcui.json?idtype=NDC&id={ndc}` | GET | Find RXCUI by NDC |
| `/relatedndc.json?ndc={ndc}` | GET | Get related NDCs |

#### Sample API Response

```json
{
  "idGroup": {
    "name": "metformin",
    "rxnormId": ["6809"]
  }
}
```

```json
{
  "rxnormdata": {
    "conceptGroup": [{
      "tty": "SCD",
      "conceptProperties": [{
        "rxcui": "861004",
        "name": "metformin hydrochloride 500 MG Oral Tablet",
        "synonym": "",
        "tty": "SCD",
        "language": "ENG",
        "suppress": "N",
        "umlscui": ""
      }]
    }]
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

**Base URL:** `https://rxnav.nlm.nih.gov/REST/rxclass/`

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

#### Sample Response

```json
{
  "rxclassMinConceptList": {
    "rxclassMinConcept": [{
      "classId": "N0000175503",
      "className": "Biguanides",
      "classType": "EPC"
    }, {
      "classId": "A10BA02",
      "className": "metformin",
      "classType": "ATC"
    }]
  }
}
```

---

## 4. ClinicalTrials.gov

**Website:** https://clinicaltrials.gov
**Maintainer:** National Library of Medicine (NLM)
**API Version:** 2.0 (REST API with OpenAPI 3.0 specification)

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
        "fullName": "Pharmaceutical Company Inc",
        "class": "INDUSTRY"
      },
      "briefTitle": "Study of Drug X in Type 2 Diabetes",
      "officialTitle": "A Phase 3, Randomized, Double-Blind, Placebo-Controlled Study...",
      "acronym": "DRUG-X-DM"
    },
    "statusModule": {
      "statusVerifiedDate": "2025-01",
      "overallStatus": "RECRUITING",
      "expandedAccessInfo": {
        "hasExpandedAccess": false
      },
      "startDateStruct": {
        "date": "2024-01-15",
        "type": "ACTUAL"
      },
      "primaryCompletionDateStruct": {
        "date": "2026-12-31",
        "type": "ESTIMATED"
      },
      "completionDateStruct": {
        "date": "2027-06-30",
        "type": "ESTIMATED"
      },
      "studyFirstSubmitDate": "2023-12-01",
      "studyFirstPostDateStruct": {
        "date": "2023-12-15",
        "type": "ACTUAL"
      }
    },
    "sponsorCollaboratorsModule": {
      "responsibleParty": {
        "type": "SPONSOR"
      },
      "leadSponsor": {
        "name": "Pharmaceutical Company Inc",
        "class": "INDUSTRY"
      },
      "collaborators": [{
        "name": "Academic Medical Center",
        "class": "OTHER"
      }]
    },
    "oversightModule": {
      "isFdaRegulatedDrug": true,
      "isFdaRegulatedDevice": false,
      "isUnapprovedDevice": false
    },
    "descriptionModule": {
      "briefSummary": "This study evaluates the efficacy and safety...",
      "detailedDescription": "Detailed description of the study design..."
    },
    "conditionsModule": {
      "conditions": ["Type 2 Diabetes Mellitus", "Hyperglycemia"],
      "keywords": ["diabetes", "glycemic control", "HbA1c"]
    },
    "designModule": {
      "studyType": "INTERVENTIONAL",
      "phases": ["PHASE3"],
      "designInfo": {
        "allocation": "RANDOMIZED",
        "interventionModel": "PARALLEL",
        "primaryPurpose": "TREATMENT",
        "maskingInfo": {
          "masking": "QUADRUPLE",
          "whoMasked": ["PARTICIPANT", "CARE_PROVIDER", "INVESTIGATOR", "OUTCOMES_ASSESSOR"]
        }
      },
      "enrollmentInfo": {
        "count": 500,
        "type": "ESTIMATED"
      }
    },
    "armsInterventionsModule": {
      "armGroups": [{
        "label": "Drug X 10mg",
        "type": "EXPERIMENTAL",
        "description": "Drug X 10mg once daily"
      }, {
        "label": "Placebo",
        "type": "PLACEBO_COMPARATOR",
        "description": "Matching placebo once daily"
      }],
      "interventions": [{
        "type": "DRUG",
        "name": "Drug X",
        "description": "Oral tablet, 10mg once daily",
        "armGroupLabels": ["Drug X 10mg"],
        "otherNames": ["Investigational Drug X"]
      }]
    },
    "outcomesModule": {
      "primaryOutcomes": [{
        "measure": "Change in HbA1c from Baseline to Week 24",
        "description": "Change in glycated hemoglobin",
        "timeFrame": "Baseline to Week 24"
      }],
      "secondaryOutcomes": [{
        "measure": "Proportion of Subjects Achieving HbA1c < 7%",
        "description": "Percentage of participants reaching target",
        "timeFrame": "Week 24"
      }]
    },
    "eligibilityModule": {
      "eligibilityCriteria": "Inclusion Criteria:\n- Adults aged 18-75\n- Diagnosed with T2DM\n- HbA1c 7.5-10.5%\n\nExclusion Criteria:\n- Type 1 diabetes\n- eGFR < 45 mL/min/1.73m2",
      "healthyVolunteers": false,
      "sex": "ALL",
      "minimumAge": "18 Years",
      "maximumAge": "75 Years",
      "stdAges": ["ADULT", "OLDER_ADULT"]
    },
    "contactsLocationsModule": {
      "centralContacts": [{
        "name": "Clinical Trial Contact",
        "role": "CONTACT",
        "phone": "1-800-555-0123",
        "email": "trials@pharma.com"
      }],
      "locations": [{
        "facility": "Research Hospital",
        "city": "Boston",
        "state": "Massachusetts",
        "country": "United States",
        "status": "RECRUITING"
      }]
    }
  }
}
```

#### Results Section

```json
{
  "resultsSection": {
    "participantFlowModule": {
      "preAssignmentDetails": "500 subjects enrolled, 480 randomized",
      "recruitmentDetails": "Recruited from 50 sites across the US",
      "groups": [{
        "id": "FG000",
        "title": "Drug X 10mg",
        "description": "Drug X 10mg once daily for 24 weeks"
      }],
      "periods": [{
        "title": "Overall Study",
        "milestones": [{
          "type": "STARTED",
          "achievements": [{"groupId": "FG000", "numSubjects": "240"}]
        }]
      }]
    },
    "baselineCharacteristicsModule": {
      "populationDescription": "All randomized subjects",
      "groups": [{
        "id": "BG000",
        "title": "Drug X 10mg",
        "description": "Treatment arm"
      }],
      "measures": [{
        "title": "Age",
        "paramType": "MEAN",
        "dispersionType": "STANDARD_DEVIATION",
        "unitOfMeasure": "years",
        "classes": [{
          "categories": [{
            "measurements": [{"groupId": "BG000", "value": "55.2", "spread": "10.1"}]
          }]
        }]
      }]
    },
    "outcomeMeasuresModule": {
      "outcomeMeasures": [{
        "type": "PRIMARY",
        "title": "Change in HbA1c from Baseline to Week 24",
        "description": "Change in glycated hemoglobin",
        "populationDescription": "Intent-to-treat population",
        "reportingStatus": "POSTED",
        "paramType": "LEAST_SQUARES_MEAN",
        "dispersionType": "95% Confidence Interval",
        "unitOfMeasure": "percentage points",
        "timeFrame": "Baseline to Week 24",
        "groups": [{
          "id": "OG000",
          "title": "Drug X 10mg",
          "description": "Treatment arm"
        }],
        "classes": [{
          "categories": [{
            "measurements": [{
              "groupId": "OG000",
              "value": "-1.2",
              "lowerLimit": "-1.4",
              "upperLimit": "-1.0"
            }]
          }]
        }],
        "analyses": [{
          "groupIds": ["OG000", "OG001"],
          "statisticalMethod": "ANCOVA",
          "pValue": "<0.001"
        }]
      }]
    },
    "adverseEventsModule": {
      "frequencyThreshold": "5",
      "timeFrame": "24 weeks",
      "description": "Treatment-emergent adverse events",
      "eventGroups": [{
        "id": "EG000",
        "title": "Drug X 10mg",
        "seriousNumAffected": 5,
        "seriousNumAtRisk": 240,
        "otherNumAffected": 45,
        "otherNumAtRisk": 240
      }],
      "seriousEvents": [{
        "term": "Myocardial infarction",
        "organSystem": "Cardiac disorders",
        "assessmentType": "SYSTEMATIC_ASSESSMENT",
        "stats": [{"groupId": "EG000", "numEvents": 2, "numAffected": 2, "numAtRisk": 240}]
      }],
      "otherEvents": [{
        "term": "Nausea",
        "organSystem": "Gastrointestinal disorders",
        "assessmentType": "NON_SYSTEMATIC_ASSESSMENT",
        "stats": [{"groupId": "EG000", "numEvents": 15, "numAffected": 12, "numAtRisk": 240}]
      }]
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

---

### 4.3 Outcome Measures

#### Primary Outcome Structure

```json
{
  "primaryOutcomes": [{
    "measure": "Change in Primary Endpoint",
    "description": "Detailed description of outcome measure",
    "timeFrame": "Baseline to Week 12"
  }]
}
```

#### Secondary Outcome Structure

```json
{
  "secondaryOutcomes": [{
    "measure": "Secondary Endpoint Name",
    "description": "Description",
    "timeFrame": "Week 12"
  }]
}
```

#### Result Outcome Structure

```json
{
  "outcomeMeasures": [{
    "type": "PRIMARY",
    "title": "Outcome Title",
    "description": "Description...",
    "populationDescription": "Analysis population...",
    "reportingStatus": "POSTED",
    "paramType": "MEAN",
    "dispersionType": "STANDARD_DEVIATION",
    "unitOfMeasure": "units",
    "timeFrame": "Time frame",
    "groups": [{
      "id": "OG000",
      "title": "Treatment Arm",
      "description": "Description"
    }],
    "classes": [{
      "categories": [{
        "measurements": [{
          "groupId": "OG000",
          "value": "10.5",
          "spread": "2.3"
        }]
      }]
    }]
  }]
}
```

---

### 4.4 API v2 Endpoints

**Base URL:** `https://clinicaltrials.gov/api/v2/`

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
https://clinicaltrials.gov/api/v2/studies?query.term=metformin&fields=NCTId,BriefTitle,OverallStatus,Phase&pageSize=50

# Find Phase 3 trials for specific condition
https://clinicaltrials.gov/api/v2/studies?query.cond=type+2+diabetes&aggFilters=phase:3&filter.overallStatus=RECRUITING,ACTIVE_NOT_RECRUITING
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
| EARLY_PHASE1 | Early Phase 1 (formerly Phase 0) |
| PHASE1 | Phase 1 |
| PHASE2 | Phase 2 |
| PHASE3 | Phase 3 |
| PHASE4 | Phase 4 (post-marketing) |
| NA | Not Applicable (non-drug studies) |

---

## 5. FDA Orange Book

**Website:** https://www.fda.gov/drugs/drug-approvals-and-databases/orange-book-data-files
**Maintainer:** FDA Center for Drug Evaluation and Research (CDER)
**Format:** Tilde-delimited (~) ASCII text files

### 5.1 Data File Structure

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

**Sample Record:**
```
METFORMIN HYDROCHLORIDE~TABLET;ORAL~GLUCOPHAGE~BMS~500MG~N~020357~001~AB~Mar 03, 1995~Y~Y~RX~Bristol-Myers Squibb Company
```

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

**Sample Record:**
```
N~020357~001~5614492~Dec 31, 2012~Y~N~~N~Jan 15, 1996
```

#### Exclusivity.txt

| Field | Type | Description |
|-------|------|-------------|
| Appl_Type | char | N=NDA, A=ANDA |
| Appl_No | string | Application number (6 digits) |
| Product_No | string | Product number (3 digits) |
| Exclusivity_Code | string | Exclusivity type code |
| Exclusivity_Date | string | Expiration date (Mmm DD, YYYY) |

**Sample Record:**
```
N~020357~001~NCE~Mar 03, 2000
```

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
| GAIN | 5 years | Qualified Infectious Disease Product (added) |

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

---

## 6. ID System Cross-Reference

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
| MedDRA | Medical Dictionary for Regulatory Activities | Adverse events coding |
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

## 7. Integration Notes for Knowledge Base

### Cross-Reference Mapping Strategy

| Source | Target | Join Key |
|--------|--------|----------|
| OpenFDA Label | RxNorm | openfda.rxcui |
| OpenFDA Label | DailyMed | openfda.spl_set_id |
| OpenFDA NDC | Orange Book | application_number |
| DailyMed | RxNorm | rxnorm_mappings index |
| ClinicalTrials.gov | OpenFDA | Drug name fuzzy match |
| Orange Book | OpenFDA | Appl_No to application_number |

### Recommended ID Resolution Order

1. **RXCUI** - Most reliable for drug concept matching
2. **NDC** - Product/package level identification
3. **UNII** - Active ingredient matching
4. **Application Number (NDA/ANDA/BLA)** - Regulatory approval linkage
5. **SPL Set ID** - Label version tracking

### Data Quality Considerations

- **OpenFDA FAERS:** Contains duplicates; use safetyreportid + version for deduplication
- **Drug Names:** Highly variable; use RxNorm normalization
- **NDC Codes:** May be formatted with/without hyphens; normalize to 11-digit format
- **Dates:** Multiple formats (YYYYMMDD, Mmm DD, YYYY); standardize to ISO 8601

---

## 8. Licensing and Terms

### Public Domain Status

All data covered in this document is **Public Domain** as work of the US Government:

| Database | License | Attribution |
|----------|---------|-------------|
| OpenFDA | Public Domain | FDA disclaimer required |
| DailyMed | Public Domain | NLM citation appreciated |
| RxNorm | Public Domain (UMLS terms apply to some sources) | NLM citation |
| ClinicalTrials.gov | Public Domain | NLM citation |
| Orange Book | Public Domain | FDA citation |

### API Terms of Service

- **OpenFDA:** Rate limits apply (240/min without key, 120K/day with key)
- **RxNorm:** Unlimited non-commercial; commercial use requires UMLS license review
- **ClinicalTrials.gov:** No explicit limits; reasonable use expected
- **DailyMed:** Bulk downloads preferred over API for large extractions

---

## 9. References

### OpenFDA
- [Drug Adverse Event API](https://open.fda.gov/apis/drug/event/)
- [Drug Label API](https://open.fda.gov/apis/drug/label/)
- [Drug NDC API](https://open.fda.gov/apis/drug/ndc/)
- [Drug Enforcement API](https://open.fda.gov/apis/drug/enforcement/)
- [OpenFDA Query Syntax](https://open.fda.gov/apis/query-syntax/)
- [GitHub: FDA/openfda](https://github.com/FDA/openfda)

### DailyMed
- [DailyMed SPL Resources](https://dailymed.nlm.nih.gov/dailymed/spl-resources.cfm)
- [FDA SPL Resources](https://www.fda.gov/industry/fda-data-standards-advisory-board/structured-product-labeling-resources)
- [Section Headings (LOINC)](https://www.fda.gov/industry/structured-product-labeling-resources/section-headings-loinc)

### RxNorm
- [RxNorm Overview](https://www.nlm.nih.gov/research/umls/rxnorm/overview.html)
- [RxNorm Technical Documentation](https://www.nlm.nih.gov/research/umls/rxnorm/docs/techdoc.html)
- [RxNorm API Documentation](https://lhncbc.nlm.nih.gov/RxNav/APIs/RxNormAPIs.html)
- [RxClass Overview](https://lhncbc.nlm.nih.gov/RxNav/applications/RxClassIntro.html)

### ClinicalTrials.gov
- [API Documentation](https://clinicaltrials.gov/data-api/api)
- [Study Data Structure](https://clinicaltrials.gov/data-api/about-api/study-data-structure)

### FDA Orange Book
- [Orange Book Data Files](https://www.fda.gov/drugs/drug-approvals-and-databases/orange-book-data-files)
- [Orange Book Preface](https://www.fda.gov/drugs/development-approval-process-drugs/orange-book-preface)
- [Patents and Exclusivity FAQ](https://www.fda.gov/drugs/development-approval-process-drugs/frequently-asked-questions-patents-and-exclusivity)

---

## 10. Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial comprehensive FDA regulatory schema documentation |

---

*Document Version: 1.0*
*Last Updated: January 2026*
