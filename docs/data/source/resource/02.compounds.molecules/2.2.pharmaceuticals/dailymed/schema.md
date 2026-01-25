---
id: schema-dailymed
title: "DailyMed Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, fda, drug-labels, spl, regulatory]
---

# DailyMed - FDA Drug Labeling Schema

**Document ID:** SCHEMA-DAILYMED
**Version:** 1.0
**Source Version:** Current (continuously updated)

---

## TL;DR

DailyMed provides FDA-approved drug labeling in Structured Product Labeling (SPL) XML format. The schema organizes drug labels into hierarchical sections (indications, dosage, warnings, etc.) with standardized identifiers (NDC, Set ID) enabling programmatic extraction of clinical drug information.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Drug Labels | 140,000+ | SPL documents |
| Active Labels | 100,000+ | Current products |
| Drug Products | 45,000+ | Unique drugs |
| NDC Numbers | 200,000+ | Package identifiers |
| Update Frequency | Daily | NLM |

---

## Entity Relationship Overview

```
SPL Documents (1) ←→ (many) Products (many) ←→ (many) NDC Codes
       ↓                         ↓
   Set ID (UUID)            Drug Names/Strengths

SPL Documents (1) ←→ (many) Sections
                               ↓
                     Label content (XML)

SPL Documents (1) ←→ (1) Application (NDA/ANDA)
```

---

## Core Tables/Entities

### spl_documents

**Description:** Root elements representing drug label documents.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| set_id | uuid | Yes | Unique document identifier |
| version_number | integer | Yes | Label version |
| effective_time | date | Yes | Label effective date |
| document_type | string | Yes | Human rx, OTC, vaccine, etc. |
| title | string | Yes | Drug product title |
| marketing_status | string | No | RX, OTC, DISCONTINUED |

### products

**Description:** Drug products described in labels.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| product_id | string | Yes | Internal identifier |
| set_id | uuid | Yes | Foreign key to spl_documents |
| name | string | Yes | Product name |
| nonproprietary_name | string | No | Generic name |
| dose_form | string | No | Tablet, injection, etc. |
| route | string | No | Oral, intravenous, etc. |
| marketing_category | string | No | NDA, ANDA, BLA |
| application_number | string | No | NDA/ANDA/BLA number |

### ndc_codes

**Description:** National Drug Code identifiers for products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ndc | string | Yes | NDC code (11 digits) |
| product_id | string | Yes | Foreign key to products |
| package_description | string | No | Package type |
| package_quantity | integer | No | Count per package |

### active_ingredients

**Description:** Active substances in drug products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ingredient_id | string | Yes | Internal identifier |
| product_id | string | Yes | Foreign key to products |
| substance_name | string | Yes | Active ingredient name |
| strength | string | No | Amount per unit |
| strength_unit | string | No | mg, ml, etc. |
| unii | string | No | FDA Unique Ingredient ID |

### spl_sections

**Description:** Content sections within drug labels.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| section_id | string | Yes | Internal identifier |
| set_id | uuid | Yes | Foreign key to spl_documents |
| section_code | string | Yes | LOINC code |
| section_name | string | Yes | Section title |
| text | xml | Yes | Section content |

---

## SPL Section Codes (LOINC)

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

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /drugnames.json | GET | List all drug names |
| /rxcui/{rxcui} | GET | Get label by RxNorm ID |
| /ndc/{ndc} | GET | Get label by NDC |
| /setid/{set_id} | GET | Get label by Set ID |
| /spellingsuggestions.json | GET | Drug name suggestions |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | SPL XML (HL7 standard) |
| Alternative | PDF (human-readable) |
| API Response | JSON, XML |
| Encoding | UTF-8 |

---

## Sample Record

```xml
<document>
  <id root="12345678-1234-1234-1234-123456789012"/>
  <effectiveTime value="20240101"/>
  <title>TYLENOL (acetaminophen) tablet</title>
  <component>
    <structuredBody>
      <component>
        <section>
          <code code="34067-9" displayName="INDICATIONS &amp; USAGE"/>
          <text>
            For temporary relief of minor aches and pains...
          </text>
        </section>
      </component>
    </structuredBody>
  </component>
</document>
```

---

## Glossary

| Term | Definition |
|------|------------|
| Set ID | Unique UUID identifying a label family |
| NDC | National Drug Code (10-11 digit identifier) |
| SPL | Structured Product Labeling (XML format) |
| LOINC | Logical Observation Identifiers Names and Codes |
| UNII | Unique Ingredient Identifier |
| NDA | New Drug Application |
| ANDA | Abbreviated New Drug Application |

---

## References

1. DailyMed: https://dailymed.nlm.nih.gov/
2. SPL Resources: https://www.fda.gov/industry/fda-resources-data-standards/structured-product-labeling-resources
3. API Documentation: https://dailymed.nlm.nih.gov/dailymed/app-support-web-services.cfm
