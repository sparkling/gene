---
id: schema-orange.book
title: "FDA Orange Book Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, fda, approved-drugs, patents, therapeutic-equivalence]
---

# FDA Orange Book Schema

**Document ID:** SCHEMA-ORANGE-BOOK
**Version:** 1.0
**Source Version:** Current (monthly updates)

---

## TL;DR

The FDA Orange Book contains approved drug products with therapeutic equivalence evaluations, patent information, and exclusivity data. The schema organizes products by application number (NDA/ANDA) with related patent listings and exclusivity records, supporting generic drug development and regulatory intelligence.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Approved Products | 34,000+ | Products table |
| Active Ingredients | 4,000+ | Distinct ingredients |
| Patent Listings | 45,000+ | Patents table |
| Exclusivity Records | 8,000+ | Exclusivity table |
| Update Frequency | Monthly | FDA |

---

## Entity Relationship Overview

```
Products (1) ←→ (many) Patents
    ↓                    ↓
Application #      Patent Number, Expiration

Products (1) ←→ (many) Exclusivity
    ↓                    ↓
Application #      Exclusivity Code, Expiration

Products (1) ←→ (1) TE_Code
                     ↓
            Therapeutic Equivalence Rating
```

---

## Core Tables/Entities

### products

**Description:** FDA-approved drug products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| appl_no | string | Yes | NDA/ANDA/BLA number |
| product_no | string | Yes | Product number within application |
| form | string | Yes | Dosage form |
| strength | string | Yes | Strength/concentration |
| reference_drug | string | No | Y/N if reference listed drug |
| drug_name | string | Yes | Trade name |
| active_ingredient | string | Yes | Active substance |
| reference_standard | string | No | Y/N if RS designated |
| te_code | string | No | Therapeutic equivalence code |
| approval_date | date | No | FDA approval date |
| rld | string | No | Reference Listed Drug flag |
| rs | string | No | Reference Standard flag |
| type | string | No | RX, OTC, DISCN |
| applicant | string | No | Company name |
| applicant_full_name | string | No | Full applicant name |

### patents

**Description:** Patent information for drug products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| appl_no | string | Yes | Foreign key to products |
| product_no | string | Yes | Foreign key to products |
| patent_no | string | Yes | Patent number |
| patent_expire_date | date | Yes | Expiration date |
| drug_substance_flag | string | No | Y/N covers drug substance |
| drug_product_flag | string | No | Y/N covers drug product |
| patent_use_code | string | No | Patent use code |
| delist_flag | string | No | Delist request flag |

### exclusivity

**Description:** Market exclusivity periods for products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| appl_no | string | Yes | Foreign key to products |
| product_no | string | Yes | Foreign key to products |
| exclusivity_code | string | Yes | Exclusivity type code |
| exclusivity_date | date | Yes | Expiration date |

---

## Therapeutic Equivalence Codes

| Code | Category | Description |
|------|----------|-------------|
| AA | A-rated | No bioequivalence issues |
| AB | A-rated | Bioequivalence demonstrated |
| AN | A-rated | Aerosol solutions |
| AO | A-rated | Injectable oil solutions |
| AP | A-rated | Injectable aqueous solutions |
| AT | A-rated | Topical products |
| BC | B-rated | Extended-release, not equivalent |
| BD | B-rated | Active ingredient studies needed |
| BE | B-rated | Delayed-release enteric, studies needed |
| BN | B-rated | Extended-release, no determination |
| BP | B-rated | Active ingredient potential problem |
| BR | B-rated | Suppositories/enemas, studies needed |
| BS | B-rated | Standard studies deficient |
| BT | B-rated | Topical, studies needed |
| BX | B-rated | Insufficient data |

---

## Exclusivity Codes

| Code | Description | Duration |
|------|-------------|----------|
| NCE | New Chemical Entity | 5 years |
| NP | New Product | 3 years |
| I | Orphan Drug Exclusivity | 7 years |
| ODE | Orphan Drug Exclusivity | 7 years |
| PED | Pediatric Exclusivity | 6 months |
| M | New formulation/indication | 3 years |
| PC | Patent Challenge | 180 days |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | ZIP archive with text files |
| Alternative | Web search interface |
| Encoding | ASCII |
| Delimiter | ~ (tilde) separated |

---

## Sample Record

```
APPL_NO~PRODUCT_NO~FORM~STRENGTH~REFERENCE_DRUG~DRUG_NAME~ACTIVE_INGREDIENT~TE_CODE
NDA020702~001~TABLET;ORAL~10MG~Y~LIPITOR~ATORVASTATIN CALCIUM~AB

PATENT_NO~APPL_NO~PRODUCT_NO~PATENT_EXPIRE_DATE~DRUG_SUBSTANCE~DRUG_PRODUCT
7456789~NDA020702~001~20251231~Y~N

EXCLUSIVITY_CODE~APPL_NO~PRODUCT_NO~EXCLUSIVITY_DATE
NCE~NDA020702~001~20080101
```

---

## Glossary

| Term | Definition |
|------|------------|
| NDA | New Drug Application |
| ANDA | Abbreviated New Drug Application |
| BLA | Biologics License Application |
| TE Code | Therapeutic Equivalence evaluation code |
| RLD | Reference Listed Drug |
| NCE | New Chemical Entity |
| Exclusivity | Period of market protection |

---

## References

1. FDA Orange Book: https://www.accessdata.fda.gov/scripts/cder/ob/
2. Orange Book Preface: https://www.fda.gov/drugs/development-approval-process-drugs/orange-book-preface
3. Data Files: https://www.fda.gov/drugs/drug-approvals-and-databases/orange-book-data-files
