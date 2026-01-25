# ConsumerLab - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | consumerlab |
| **Name** | ConsumerLab |
| **Total Fields** | 40+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| review_id | string | 1:1 | Yes | Unique review identifier | `CL-2024-VIT-D-001` |
| product_name | string | 1:1 | Yes | Full product name | `Vitamin D3 5000 IU` |
| brand | string | 1:1 | Yes | Brand name | `Nature Made` |
| upc | string | 1:1 | No | Universal Product Code | `012345678901` |

### Testing Results

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| test_date | date | 1:1 | No | Date of testing | `2024-01-15` |
| approval_status | string | 1:1 | Yes | Overall status | `Approved` |
| quality_seal | boolean | 1:1 | No | CL seal awarded | `true` |
| percent_claim | number | 1:1 | No | Measured vs labeled | `102.4` |

### Heavy Metal Testing

| Field Name | Data Type | Cardinality | Required | Description | Allowed Values |
|------------|-----------|-------------|----------|-------------|----------------|
| lead_mcg_daily | number | 1:1 | No | Lead per daily dose | <0.5 (USP limit) |
| arsenic_mcg_daily | number | 1:1 | No | Arsenic per daily dose | <10.0 (USP limit) |
| cadmium_mcg_daily | number | 1:1 | No | Cadmium per daily dose | <4.1 (USP limit) |
| mercury_mcg_daily | number | 1:1 | No | Mercury per daily dose | <2.0 (USP limit) |

---

## Enumerations

### Approval Status

| Value | Description |
|-------|-------------|
| Approved | Meets all quality criteria |
| Not Approved | Failed one or more tests |
| Voluntarily Withdrawn | Removed by manufacturer |
| Not Tested | Awaiting evaluation |

### Test Status

| Value | Description |
|-------|-------------|
| Pass | Meets or exceeds criteria |
| Fail | Below acceptable threshold |
| ND | Not detected |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CL | ConsumerLab | Testing organization |
| USP | United States Pharmacopeia | Quality standards |
| UPC | Universal Product Code | Barcode identifier |
| IU | International Unit | Vitamin measurement |
| ND | Not Detected | Below detection limit |
| mcg | Microgram | Weight unit (1/1000 mg) |

---

## Data Quality Notes

1. **Proprietary data**: Access requires paid subscription
2. **Testing methodology**: Uses USP standards for limits
3. **Potency range**: Acceptable is 90-110% of label claim
4. **Disintegration**: Must occur within 30 minutes per USP

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
