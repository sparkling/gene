# Natural Medicines - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | natural-medicines |
| **Name** | Natural Medicines Comprehensive Database |
| **Total Fields** | 50+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| monograph_id | string | 1:1 | Yes | Unique ID | `NM-001234` |
| name | string | 1:1 | Yes | Primary name | `Turmeric` |
| scientific_names | array | 1:N | No | Latin names | `Curcuma longa` |
| common_names | array | 1:N | No | Alternative names | `Curcumin` |
| family | string | 1:1 | No | Plant family | `Zingiberaceae` |

### Safety & Efficacy

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| safety_grade | string | 1:1 | Yes | Overall safety | `Likely Safe` |
| effectiveness_ratings | array | 1:N | No | By condition | See enumerations |
| evidence_level | string | 1:1 | No | Evidence quality | `A`, `B`, `C`, `D` |

### Interactions

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| drug | string | 1:N | No | Interacting drug | `Warfarin` |
| severity | string | 1:1 | No | Interaction severity | `Major` |
| mechanism | string | 1:1 | No | Interaction mechanism | `CYP3A4 induction` |

---

## Enumerations

### Safety Grades

| Value | Description |
|-------|-------------|
| Likely Safe | No significant concerns when used appropriately |
| Possibly Safe | Limited data, no significant concerns reported |
| Possibly Unsafe | Reports of adverse effects; use with caution |
| Likely Unsafe | Significant adverse effects documented |
| Unsafe | Do not use; serious adverse effects |
| Insufficient Information | Not enough safety data available |

### Effectiveness Ratings

| Value | Evidence Level | Description |
|-------|----------------|-------------|
| Effective | A | Strong evidence from well-designed trials |
| Likely Effective | A-B | Good evidence from clinical trials |
| Possibly Effective | B | Some evidence, more research needed |
| Possibly Ineffective | B | Evidence suggests lack of benefit |
| Likely Ineffective | A-B | Good evidence of no benefit |
| Ineffective | A | Strong evidence of no benefit |
| Insufficient Evidence | C-D | Not enough data to evaluate |

### Interaction Severity

| Value | Action Required |
|-------|-----------------|
| Contraindicated | Do not use together |
| Major | Use alternative or monitor closely |
| Moderate | Monitor and consider alternatives |
| Minor | Be aware; usually no action needed |
| Unknown | Use caution |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NM | Natural Medicines | Database name |
| TRC | Therapeutic Research Center | Publisher |
| CYP | Cytochrome P450 | Drug metabolism enzymes |
| RCT | Randomized Controlled Trial | Study design |
| UNII | Unique Ingredient Identifier | FDA identifier |
| EHR | Electronic Health Record | Integration target |
| INR | International Normalized Ratio | Warfarin monitoring |

---

## Data Quality Notes

1. **Proprietary data**: Requires paid subscription
2. **Evidence-based**: Systematic review of clinical literature
3. **Regular updates**: Continuously updated monographs
4. **EHR integration**: Available for healthcare systems

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
