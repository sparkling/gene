# HMP Gut Microbiome - Data Dictionary

## Overview

This data dictionary documents HMP gut microbiome sample records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | hmp |
| **Name** | Human Microbiome Project |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| sample_id | string | Yes | Unique sample ID |
| rand_subject_id | string | No | HIPAA-compliant subject ID |
| body_site | string | Yes | stool, cecum, colon |
| visit_number | integer | No | Visit sequence number |
| study | string | No | HMP1, IBDMDB, T2D, MOMS-PI |

---

## Studies

| Study | Focus | Description |
|-------|-------|-------------|
| HMP1 | Healthy cohort | Baseline characterization |
| IBDMDB | IBD | Inflammatory bowel disease |
| T2D | Prediabetes | Type 2 diabetes |
| MOMS-PI | Pregnancy | Preterm birth risk |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| HMP | Human Microbiome Project |
| iHMP | Integrative HMP |
| OSDF | Open Science Data Framework |
| MIxS | Minimum Information about any Sequence |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
