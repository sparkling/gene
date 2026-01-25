# FDA Orange Book - Data Dictionary

## Overview

This data dictionary documents the schema for FDA Orange Book (Approved Drug Products with Therapeutic Equivalence Evaluations).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | orange.book |
| **Name** | FDA Orange Book |
| **Parent** | 2.2.pharmaceuticals |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Product Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| appl_no | string | 1:1 | Yes | NDA/ANDA/BLA number | NDA020702 |
| product_no | string | 1:1 | Yes | Product number within application | 001 |
| drug_name | string | 1:1 | Yes | Trade name | LIPITOR |
| active_ingredient | string | 1:1 | Yes | Active substance | ATORVASTATIN CALCIUM |
| form | string | 1:1 | Yes | Dosage form | TABLET;ORAL |
| strength | string | 1:1 | Yes | Strength/concentration | 10MG |
| te_code | string | 1:1 | No | Therapeutic equivalence code | AB |
| reference_drug | string | 1:1 | No | Reference listed drug flag | Y/N |
| rld | string | 1:1 | No | RLD designation | Y/N |
| rs | string | 1:1 | No | Reference standard flag | Y/N |
| approval_date | date | 1:1 | No | FDA approval date | 1996-12-17 |
| type | string | 1:1 | No | Product type | RX, OTC, DISCN |
| applicant | string | 1:1 | No | Company name | PFIZER INC |

### Patent Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| patent_no | string | 1:1 | Yes | Patent number | 7456789 |
| patent_expire_date | date | 1:1 | Yes | Expiration date | 2025-12-31 |
| drug_substance_flag | string | 1:1 | No | Covers drug substance | Y/N |
| drug_product_flag | string | 1:1 | No | Covers drug product | Y/N |
| patent_use_code | string | 1:1 | No | Patent use code | U-1234 |
| delist_flag | string | 1:1 | No | Delist request | Y/N |

### Exclusivity Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| exclusivity_code | string | 1:1 | Yes | Exclusivity type | NCE, NP, I, PED |
| exclusivity_date | date | 1:1 | Yes | Expiration date | 2028-01-15 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Application Number | NDA/ANDA/BLA + digits | NDA020702 | FDA application ID |
| Product Number | 3 digits | 001 | Product within application |
| Patent Number | 7 digits | 7456789 | US patent number |
| TE Code | 2-letter code | AB | Therapeutic equivalence |

---

## Enumerations

### Therapeutic Equivalence Codes

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
| BE | B-rated | Delayed-release, studies needed |
| BN | B-rated | Extended-release, no determination |
| BP | B-rated | Active ingredient potential problem |
| BR | B-rated | Suppositories/enemas, studies needed |
| BS | B-rated | Standard studies deficient |
| BT | B-rated | Topical, studies needed |
| BX | B-rated | Insufficient data |

### Exclusivity Codes

| Code | Description | Duration |
|------|-------------|----------|
| NCE | New Chemical Entity | 5 years |
| NP | New Product (clinical studies) | 3 years |
| I | Orphan Drug Exclusivity | 7 years |
| ODE | Orphan Drug Exclusivity | 7 years |
| PED | Pediatric Exclusivity | 6 months |
| M | New formulation/indication | 3 years |
| PC | Patent Challenge (first generic) | 180 days |
| CGT | Competitive Generic Therapy | 180 days |

### Application Types

| Type | Description |
|------|-------------|
| NDA | New Drug Application (brand) |
| ANDA | Abbreviated NDA (generic) |
| BLA | Biologics License Application |
| N | New molecular entity |
| S | Supplement to approved NDA |

### Product Types

| Type | Description |
|------|-------------|
| RX | Prescription drug |
| OTC | Over-the-counter |
| DISCN | Discontinued |

---

## Entity Relationships

### Product to Patents
- **Cardinality:** 1:N
- **Description:** Products may have multiple patents listed
- **Key Fields:** appl_no, product_no, patent_no

### Product to Exclusivity
- **Cardinality:** 1:N
- **Description:** Products may have multiple exclusivity periods
- **Key Fields:** appl_no, product_no, exclusivity_code

### Product to TE Code
- **Cardinality:** 1:1
- **Description:** Each product has one therapeutic equivalence rating
- **Key Fields:** appl_no, product_no, te_code

### Application to Products
- **Cardinality:** 1:N
- **Description:** One application may cover multiple products
- **Key Fields:** appl_no

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NDA | New Drug Application | Brand drug approval |
| ANDA | Abbreviated New Drug Application | Generic approval |
| BLA | Biologics License Application | Biologic approval |
| TE | Therapeutic Equivalence | Substitutability rating |
| RLD | Reference Listed Drug | Original product |
| RS | Reference Standard | Quality benchmark |
| NCE | New Chemical Entity | Novel compound |
| ODE | Orphan Drug Exclusivity | Rare disease protection |
| PED | Pediatric Exclusivity | Pediatric study extension |
| FDA | Food and Drug Administration | Regulatory agency |
| OTC | Over-The-Counter | Non-prescription |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| DailyMed | Application # | Label information |
| RxNorm | RXCUI | Drug terminology |
| DrugBank | DrugBank ID | Drug details |
| USPTO | Patent # | Patent details |
| NDC Directory | NDC | Package codes |

---

## Data Quality Notes

1. **Monthly Updates:** Orange Book updated monthly by FDA
2. **TE Ratings:** Only drugs with multiple sources have TE codes
3. **Patent Accuracy:** Patent listings are applicant-submitted
4. **RLD Designation:** Critical for generic drug development (ANDA)
5. **Exclusivity Tracking:** Essential for generic entry planning
6. **File Format:** Tilde (~) delimited text files
7. **Public Access:** Data freely available from FDA
