# DailyMed - Data Dictionary

## Overview

This data dictionary documents the schema for DailyMed FDA drug labeling database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dailymed |
| **Name** | DailyMed |
| **Parent** | 2.2.pharmaceuticals |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### SPL Document

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| set_id | uuid | 1:1 | Yes | Unique document identifier | 12345678-1234-1234-1234-123456789012 |
| version_number | integer | 1:1 | Yes | Label version | 3 |
| effective_time | date | 1:1 | Yes | Label effective date | 2024-01-01 |
| document_type | string | 1:1 | Yes | Document category | Human rx, OTC, vaccine |
| title | string | 1:1 | Yes | Drug product title | TYLENOL (acetaminophen) tablet |
| marketing_status | string | 1:1 | No | Marketing status | RX, OTC, DISCONTINUED |

### Product Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| product_id | string | 1:1 | Yes | Internal identifier | PROD-12345 |
| name | string | 1:1 | Yes | Product name | Tylenol |
| nonproprietary_name | string | 1:1 | No | Generic name | Acetaminophen |
| dose_form | string | 1:1 | No | Dosage form | Tablet, Capsule, Injection |
| route | string | 1:1 | No | Administration route | Oral, Intravenous, Topical |
| marketing_category | string | 1:1 | No | Regulatory category | NDA, ANDA, BLA |
| application_number | string | 1:1 | No | NDA/ANDA/BLA number | NDA020702 |

### NDC Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ndc | string | 1:1 | Yes | National Drug Code | 50580-176-01 |
| package_description | string | 1:1 | No | Package type | Bottle, Blister pack |
| package_quantity | integer | 1:1 | No | Count per package | 100 |

### Active Ingredient

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ingredient_id | string | 1:1 | Yes | Internal identifier | ING-12345 |
| substance_name | string | 1:1 | Yes | Active ingredient name | ACETAMINOPHEN |
| strength | string | 1:1 | No | Amount per unit | 500 mg |
| strength_unit | string | 1:1 | No | Unit of strength | mg, ml |
| unii | string | 1:1 | No | FDA Unique Ingredient ID | 362O9ITL9D |

### Label Section

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| section_id | string | 1:1 | Yes | Internal identifier | SEC-12345 |
| section_code | string | 1:1 | Yes | LOINC code | 34067-9 |
| section_name | string | 1:1 | Yes | Section title | Indications & Usage |
| text | xml | 1:1 | Yes | Section content | XML content |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Set ID | UUID format | 12345678-1234-1234-1234-123456789012 | Document family ID |
| NDC | 10-11 digits | 50580-176-01 | National Drug Code |
| UNII | 10 characters | 362O9ITL9D | FDA substance ID |
| RXCUI | Integer | 198440 | RxNorm concept ID |
| LOINC | XX-X format | 34067-9 | Section code |

---

## Enumerations

### SPL Section Codes (LOINC)

| Code | Section Name |
|------|--------------|
| 34067-9 | Indications & Usage |
| 34068-7 | Dosage & Administration |
| 34070-3 | Contraindications |
| 34071-1 | Warnings and Precautions |
| 34084-4 | Adverse Reactions |
| 34073-7 | Drug Interactions |
| 43684-0 | Use in Specific Populations |
| 34089-3 | Description |
| 34090-1 | Clinical Pharmacology |
| 34069-5 | How Supplied |
| 34076-0 | Patient Counseling Information |
| 42229-5 | SPL Unclassified Section |

### Marketing Categories

| Category | Description |
|----------|-------------|
| NDA | New Drug Application |
| ANDA | Abbreviated NDA (generic) |
| BLA | Biologics License Application |
| OTC | Over-the-counter monograph |
| UNAPPROVED | Unapproved drug |

### Document Types

| Type | Description |
|------|-------------|
| HUMAN PRESCRIPTION DRUG LABEL | Rx products |
| HUMAN OTC DRUG LABEL | OTC products |
| VACCINE LABEL | Vaccines |
| PLASMA DERIVATIVE LABEL | Blood products |
| CELLULAR THERAPY LABEL | Cell therapies |

### Marketing Status

| Status | Description |
|--------|-------------|
| RX | Prescription |
| OTC | Over-the-counter |
| DISCONTINUED | No longer marketed |

---

## Entity Relationships

### SPL Document to Products
- **Cardinality:** 1:N
- **Description:** One label document describes multiple products
- **Key Fields:** set_id, product_id

### Product to NDC Codes
- **Cardinality:** 1:N
- **Description:** One product can have multiple NDC codes (packages)
- **Key Fields:** product_id, ndc

### Product to Active Ingredients
- **Cardinality:** 1:N
- **Description:** Products can have multiple active ingredients
- **Key Fields:** product_id, ingredient_id

### SPL Document to Sections
- **Cardinality:** 1:N
- **Description:** Labels contain multiple content sections
- **Key Fields:** set_id, section_id

### Document to Application
- **Cardinality:** N:1
- **Description:** Multiple labels may share application number
- **Key Fields:** set_id, application_number

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| SPL | Structured Product Labeling | HL7 XML format |
| NDC | National Drug Code | 10-11 digit identifier |
| NDA | New Drug Application | Brand approval |
| ANDA | Abbreviated New Drug Application | Generic approval |
| BLA | Biologics License Application | Biologic approval |
| LOINC | Logical Observation Identifiers Names and Codes | Section codes |
| UNII | Unique Ingredient Identifier | FDA substance ID |
| NLM | National Library of Medicine | DailyMed host |
| FDA | Food and Drug Administration | Regulatory body |
| OTC | Over-The-Counter | Non-prescription |
| RXCUI | RxNorm Concept Unique Identifier | Drug terminology |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| RxNorm | RXCUI | Drug terminology |
| FDA Orange Book | Application # | Approval status |
| DrugBank | DrugBank ID | Drug information |
| PubChem | CID | Chemical data |
| MeSH | MeSH ID | Medical terms |
| NCI Thesaurus | NCIt ID | Drug classification |

---

## Data Quality Notes

1. **Continuous Updates:** Labels updated daily as FDA receives new submissions
2. **Version Control:** Each label has version number for change tracking
3. **XML Format:** SPL uses HL7 CDA standard for structured content
4. **NDC Completeness:** Not all products have NDC codes listed
5. **UNII Mapping:** Active ingredients linked via FDA substance registry
6. **Public Domain:** DailyMed content is in public domain (no license restrictions)
