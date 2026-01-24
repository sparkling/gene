# 8.4 Regulatory & Legal - Data Dictionary

**Subcategory ID:** 8.4
**Subcategory Name:** Regulatory & Legal
**Data Sources:** ClinicalTrials.gov, FDA OpenFDA

## Overview

This subcategory integrates regulatory and legal data related to clinical trials and FDA drug information, including adverse event reports, drug labeling, and recall information.

---

## Unified Fields

### Core Record Information

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `record_id` | string | Yes | Primary identifier for the regulatory record | `NCT00000001`, `NDA020357` |
| `record_type` | string | Yes | Type of regulatory record | `clinical_trial`, `adverse_event`, `drug_label`, `recall` |
| `status` | string | No | Current status of the record | `RECRUITING`, `COMPLETED`, `current` |
| `title` | string | Yes | Official title of the study or product | `Study of Drug X in Type 2 Diabetes` |
| `sponsor` | string | No | Organization sponsoring the trial or product | `Pfizer Inc` |

**Source Mappings - Core Fields:**

| Field | ClinicalTrials.gov | FDA OpenFDA |
|-------|-------------------|-------------|
| `record_id` | nctId | safetyreportid / spl_id / recall_number |
| `status` | overallStatus | status |
| `title` | briefTitle / officialTitle | brand_name |
| `sponsor` | leadSponsor.name | manufacturer_name |

---

### Date Information

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `start_date` | date | No | Date of study initiation or marketing start | `2023-06-01` |
| `completion_date` | date | No | Expected or actual completion date | `2026-06-30` |

**Source Mappings - Dates:**

| Field | ClinicalTrials.gov | FDA OpenFDA |
|-------|-------------------|-------------|
| `start_date` | startDateStruct.date | marketing_start_date |
| `completion_date` | completionDateStruct.date | marketing_end_date |

---

### Conditions and Interventions

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `conditions` | array[string] | No | Medical conditions or diseases addressed | `["Type 2 Diabetes Mellitus", "Hyperglycemia"]` |
| `interventions` | array[object] | No | Drugs, devices, or procedures being tested | See structure below |

**Intervention Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `type` | string | Intervention type (DRUG, DEVICE, BIOLOGICAL, etc.) |
| `name` | string | Intervention name |
| `description` | string | Detailed description |

**Source Mappings - Conditions/Interventions:**

| Field | ClinicalTrials.gov | FDA OpenFDA |
|-------|-------------------|-------------|
| `conditions` | conditions | indications_and_usage |
| `interventions` | interventions | medicinalproduct |

---

## ClinicalTrials.gov Specific Fields

### Study Identification

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `clinical_trial_data.nct_id` | string | NCT number (National Clinical Trial identifier) | `NCT00000001` |
| `clinical_trial_data.org_study_id` | string | Sponsor-assigned study identifier | `PROTOCOL-001` |
| `clinical_trial_data.secondary_ids` | array[string] | Secondary identifiers (EUDRACT, etc.) | `["2020-001234-56"]` |
| `clinical_trial_data.acronym` | string | Study acronym | `CHAMPION` |

### Sponsor Information

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `clinical_trial_data.organization_class` | string | Sponsor type | `INDUSTRY`, `NIH`, `OTHER`, `FED`, `NETWORK` |

### Study Status

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `clinical_trial_data.overall_status` | string | Recruitment status | `NOT_YET_RECRUITING`, `RECRUITING`, `ENROLLING_BY_INVITATION`, `ACTIVE_NOT_RECRUITING`, `COMPLETED`, `SUSPENDED`, `TERMINATED`, `WITHDRAWN` |
| `clinical_trial_data.expanded_access` | boolean | Whether expanded access is available | `true` |

### Study Design

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `clinical_trial_data.study_type` | string | Study type | `INTERVENTIONAL`, `OBSERVATIONAL`, `EXPANDED_ACCESS` |
| `clinical_trial_data.phases` | array[string] | Trial phases | `["PHASE1", "PHASE2", "PHASE3", "PHASE4", "NA"]` |

**Design Object Structure:**

| Property | Data Type | Description | Example Values |
|----------|-----------|-------------|----------------|
| `allocation` | string | Allocation method | `RANDOMIZED`, `NON_RANDOMIZED`, `NA` |
| `intervention_model` | string | Study design | `PARALLEL`, `CROSSOVER`, `SEQUENTIAL`, `FACTORIAL`, `SINGLE_GROUP` |
| `masking` | string | Blinding level | `NONE`, `SINGLE`, `DOUBLE`, `TRIPLE`, `QUADRUPLE` |
| `who_masked` | array[string] | Who is blinded | `["PARTICIPANT", "INVESTIGATOR", "OUTCOMES_ASSESSOR", "CARE_PROVIDER"]` |
| `primary_purpose` | string | Purpose | `TREATMENT`, `PREVENTION`, `DIAGNOSTIC`, `SUPPORTIVE_CARE`, `SCREENING`, `HEALTH_SERVICES_RESEARCH`, `BASIC_SCIENCE`, `DEVICE_FEASIBILITY`, `OTHER` |

### Enrollment and Arms

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `clinical_trial_data.enrollment_count` | integer | Target or actual enrollment | `500` |
| `clinical_trial_data.arm_groups` | array[object] | Treatment arms | See structure below |

**Arm Group Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `label` | string | Arm label |
| `type` | string | Arm type (EXPERIMENTAL, ACTIVE_COMPARATOR, PLACEBO_COMPARATOR, etc.) |
| `description` | string | Arm description |

### Outcomes

| Field | Data Type | Description |
|-------|-----------|-------------|
| `clinical_trial_data.outcomes.primary` | array[object] | Primary endpoint measures |
| `clinical_trial_data.outcomes.secondary` | array[object] | Secondary endpoint measures |

**Outcome Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `measure` | string | Outcome measure name |
| `description` | string | Detailed description |
| `time_frame` | string | Time frame for measurement |

### Eligibility

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `clinical_trial_data.eligibility.criteria` | string | Inclusion/exclusion criteria text | `Inclusion: Age >= 18...` |
| `clinical_trial_data.eligibility.minimum_age` | string | Minimum eligible age | `18 Years` |
| `clinical_trial_data.eligibility.maximum_age` | string | Maximum eligible age | `65 Years` |
| `clinical_trial_data.eligibility.sex` | string | Eligible sex | `ALL`, `MALE`, `FEMALE` |
| `clinical_trial_data.eligibility.healthy_volunteers` | boolean | Whether healthy volunteers accepted | `false` |

### Locations

| Field | Data Type | Description |
|-------|-----------|-------------|
| `clinical_trial_data.locations` | array[object] | Study site locations |

**Location Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `facility` | string | Facility name |
| `city` | string | City name |
| `country` | string | Country name |
| `geo_point` | object | Geographic coordinates `{lat, lon}` |

### Results Data

| Field | Data Type | Description |
|-------|-----------|-------------|
| `clinical_trial_data.results.participant_flow` | object | Participant flow through study |
| `clinical_trial_data.results.baseline_characteristics` | object | Baseline demographics |
| `clinical_trial_data.results.outcome_measures` | array | Outcome data with statistics |
| `clinical_trial_data.results.adverse_events` | object | Adverse event summary |

### Regulatory Oversight

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `clinical_trial_data.fda_regulated_drug` | boolean | Whether FDA-regulated drug study | `true` |
| `clinical_trial_data.fda_regulated_device` | boolean | Whether FDA-regulated device study | `false` |

---

## FDA OpenFDA Specific Fields

### Adverse Event Data (FAERS)

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `fda_adverse_event_data.safety_report_id` | string | FAERS unique report identifier | `10003456` |
| `fda_adverse_event_data.safety_report_version` | string | Report version number | `1` |
| `fda_adverse_event_data.receive_date` | string | Date FDA received report (YYYYMMDD) | `20240115` |
| `fda_adverse_event_data.report_type` | string | Report type code | `1` (Spontaneous), `2` (Literature), `3` (Study) |

**Seriousness Details:**

| Field | Data Type | Description |
|-------|-----------|-------------|
| `fda_adverse_event_data.serious` | boolean | Whether serious adverse event |
| `fda_adverse_event_data.seriousness_details.death` | boolean | Whether resulted in death |
| `fda_adverse_event_data.seriousness_details.hospitalization` | boolean | Whether required hospitalization |
| `fda_adverse_event_data.seriousness_details.life_threatening` | boolean | Whether life-threatening |
| `fda_adverse_event_data.seriousness_details.disabling` | boolean | Whether caused disability |

**Patient Demographics:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `age` | string | Patient age at event onset |
| `sex` | string | Patient sex (0=Unknown, 1=Male, 2=Female) |
| `weight` | string | Patient weight in kg |

**Reaction Details:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `meddra_pt` | string | MedDRA preferred term for reaction |
| `outcome` | string | Reaction outcome (1-6 scale) |

**Drug Information:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `characterization` | string | Drug role (1=Suspect, 2=Concomitant, 3=Interacting) |
| `medicinal_product` | string | Drug name as reported |
| `authorization_numb` | string | NDA/ANDA number |
| `dosage_text` | string | Free-text dosage information |
| `indication` | string | Indication for use |
| `administration_route` | string | Route of administration code |
| `action_drug` | string | Action taken (withdrawn, dose reduced, etc.) |

| Field | Data Type | Description |
|-------|-----------|-------------|
| `fda_adverse_event_data.occur_country` | string | Country where event occurred (ISO 3166-1) |
| `fda_adverse_event_data.reporter_qualification` | string | Reporter type (physician, pharmacist, consumer) |

---

### Drug Label Data (SPL)

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `fda_drug_label_data.spl_id` | string | SPL document UUID | `a1b2c3d4-e5f6-7890-abcd-ef1234567890` |
| `fda_drug_label_data.set_id` | string | SPL set UUID (groups versions) | `b2c3d4e5-f6a7-8901-bcde-f23456789012` |
| `fda_drug_label_data.effective_time` | string | Label effective date (YYYYMMDD) | `20240101` |

**Label Content Sections:**

| Field | Data Type | Description |
|-------|-----------|-------------|
| `fda_drug_label_data.indications_and_usage` | array[string] | Approved therapeutic uses |
| `fda_drug_label_data.dosage_and_administration` | array[string] | Dosing instructions |
| `fda_drug_label_data.contraindications` | array[string] | Contraindication warnings |
| `fda_drug_label_data.warnings_and_cautions` | array[string] | Safety warnings |
| `fda_drug_label_data.adverse_reactions` | array[string] | Reported side effects |
| `fda_drug_label_data.drug_interactions` | array[string] | Drug-drug interactions |
| `fda_drug_label_data.boxed_warning` | array[string] | Black box warning text |
| `fda_drug_label_data.clinical_pharmacology` | array[string] | Mechanism and PK/PD |

**Product Identification:**

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `fda_drug_label_data.application_number` | string | NDA/ANDA/BLA number | `NDA020357` |
| `fda_drug_label_data.brand_name` | array[string] | Proprietary name(s) | `["GLUCOPHAGE"]` |
| `fda_drug_label_data.generic_name` | array[string] | Non-proprietary name(s) | `["METFORMIN HYDROCHLORIDE"]` |
| `fda_drug_label_data.product_ndc` | string | NDC product code (labeler-product) | `0087-6060` |
| `fda_drug_label_data.package_ndc` | string | NDC with package segment | `0087-6060-01` |
| `fda_drug_label_data.product_type` | string | Product type | `HUMAN PRESCRIPTION DRUG`, `HUMAN OTC DRUG` |
| `fda_drug_label_data.dosage_form` | string | Physical form | `TABLET`, `CAPSULE`, `INJECTION` |
| `fda_drug_label_data.route` | array[string] | Route(s) of administration | `["ORAL"]` |

**Active Ingredients:**

| Field | Data Type | Description |
|-------|-----------|-------------|
| `fda_drug_label_data.active_ingredients` | array[object] | Active ingredients with strength |

**Active Ingredient Object:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `name` | string | Ingredient name |
| `strength` | string | Ingredient strength |

**Manufacturer and Classification:**

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `fda_drug_label_data.labeler_name` | string | Manufacturer/distributor | `Bristol-Myers Squibb Company` |
| `fda_drug_label_data.dea_schedule` | string | DEA schedule | `CI`, `CII`, `CIII`, `CIV`, `CV` |
| `fda_drug_label_data.marketing_category` | string | Category | `NDA`, `ANDA`, `BLA`, `OTC MONOGRAPH FINAL` |

**Pharmacologic Classes:**

| Field | Data Type | Description |
|-------|-----------|-------------|
| `fda_drug_label_data.pharmacologic_classes.epc` | array[string] | Established pharmacologic class |
| `fda_drug_label_data.pharmacologic_classes.moa` | array[string] | Mechanism of action class |
| `fda_drug_label_data.pharmacologic_classes.pe` | array[string] | Physiologic effect class |

**External Identifiers:**

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `fda_drug_label_data.rxcui` | array[string] | RxNorm concept identifier(s) | `["6809", "860975"]` |
| `fda_drug_label_data.unii` | array[string] | Unique Ingredient Identifier(s) | `["9100L32L2N"]` |

---

### Recall Data

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `fda_recall_data.recall_number` | string | Unique recall identifier | `D-1234-2024` |
| `fda_recall_data.recall_classification` | string | Recall class | `Class I` (most serious), `Class II`, `Class III` |
| `fda_recall_data.recalling_firm` | string | Company initiating recall | `ABC Pharmaceuticals Inc` |
| `fda_recall_data.reason_for_recall` | string | Recall reason explanation | `Contamination with foreign particles` |
| `fda_recall_data.distribution_pattern` | string | Geographic distribution | `Nationwide`, `CA, NY, TX` |
| `fda_recall_data.voluntary_mandated` | string | Recall initiation type | `Voluntary`, `FDA Mandated` |
| `fda_recall_data.code_info` | string | Lot numbers and expiration dates | `Lot ABC123, Exp 12/2024` |

**Recall Classification Definitions:**

| Class | Definition |
|-------|------------|
| Class I | Dangerous or defective products that predictably could cause serious health problems or death |
| Class II | Products that might cause a temporary health problem, or pose only a slight threat of a serious nature |
| Class III | Products that are unlikely to cause any adverse health reaction, but violate FDA labeling or manufacturing laws |

---

## Source Mappings Summary

### ClinicalTrials.gov Field Mappings

| Unified Field | ClinicalTrials.gov Source |
|---------------|--------------------------|
| `record_id` | nctId |
| `status` | overallStatus |
| `title` | briefTitle / officialTitle |
| `sponsor` | leadSponsor.name |
| `start_date` | startDateStruct.date |
| `completion_date` | completionDateStruct.date |
| `conditions` | conditions |
| `interventions` | interventions |
| `clinical_trial_data.org_study_id` | orgStudyId |
| `clinical_trial_data.secondary_ids` | secondaryIds |
| `clinical_trial_data.acronym` | acronym |
| `clinical_trial_data.organization_class` | organization.class |
| `clinical_trial_data.expanded_access` | expandedAccess |
| `clinical_trial_data.study_type` | studyType |
| `clinical_trial_data.phases` | phases |
| `clinical_trial_data.design` | designInfo |
| `clinical_trial_data.enrollment_count` | enrollmentCount |
| `clinical_trial_data.arm_groups` | armGroups |
| `clinical_trial_data.outcomes.primary` | primaryOutcomes |
| `clinical_trial_data.outcomes.secondary` | secondaryOutcomes |
| `clinical_trial_data.eligibility` | eligibility |
| `clinical_trial_data.locations` | locations |
| `clinical_trial_data.results` | results |
| `clinical_trial_data.fda_regulated_drug` | fdaRegulatedDrug |
| `clinical_trial_data.fda_regulated_device` | fdaRegulatedDevice |

### FDA OpenFDA Field Mappings

| Unified Field | FDA OpenFDA Source |
|---------------|-------------------|
| `record_id` | safetyreportid / spl_id / recall_number |
| `fda_adverse_event_data.safety_report_version` | safetyreportversion |
| `fda_adverse_event_data.receive_date` | receivedate |
| `fda_adverse_event_data.report_type` | reporttype |
| `fda_adverse_event_data.serious` | serious |
| `fda_adverse_event_data.seriousness_details.death` | seriousnessdeath |
| `fda_adverse_event_data.occur_country` | occurcountry |
| `fda_adverse_event_data.patient` | patient |
| `fda_adverse_event_data.reaction` | reaction |
| `fda_adverse_event_data.drug` | drug |
| `fda_drug_label_data.spl_id` | spl_id |
| `fda_drug_label_data.set_id` | set_id |
| `fda_drug_label_data.effective_time` | effective_time |
| `fda_drug_label_data.indications_and_usage` | indications_and_usage |
| `fda_drug_label_data.brand_name` | brand_name |
| `fda_drug_label_data.generic_name` | generic_name |
| `fda_drug_label_data.product_ndc` | product_ndc |
| `fda_recall_data.recall_number` | recall_number |
| `fda_recall_data.recall_classification` | classification |
| `fda_recall_data.recalling_firm` | recalling_firm |
| `fda_recall_data.reason_for_recall` | reason_for_recall |

---

## Data Source Metadata

| Field | Data Type | Description |
|-------|-----------|-------------|
| `_source.primary_source` | string | Name of the primary data source |
| `_source.source_id` | string | Original ID in the source system |
| `_source.extraction_date` | date | Date when data was extracted |
| `_source.source_version` | string | Version of the source data/API |

---

## API Endpoints

### ClinicalTrials.gov

- **REST API v2:** `https://clinicaltrials.gov/api/v2/studies`
- **Parameters:** `query.term`, `filter.overallStatus`, `filter.geo`, `pageSize`
- **Documentation:** https://clinicaltrials.gov/data-api/api

### FDA OpenFDA

- **Drug Events:** `https://api.fda.gov/drug/event.json`
- **Drug Labels:** `https://api.fda.gov/drug/label.json`
- **Drug Recalls:** `https://api.fda.gov/drug/enforcement.json`
- **Parameters:** `search`, `count`, `limit`, `skip`
- **Documentation:** https://open.fda.gov/apis/
