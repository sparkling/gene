# HMP Body Sites - Data Dictionary

## Overview

This data dictionary documents HMP multi-body-site sample records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | hmp |
| **Name** | Human Microbiome Project (Body Sites) |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| sample_id | string | Yes | Unique sample identifier |
| rand_subject_id | string | No | HIPAA-compliant subject ID |
| body_site | string | Yes | Specific anatomical location |
| supersite | string | Yes | Major body region (5 categories) |
| fma_body_site | string | No | UBERON ontology term |
| study | string | No | HMP1, IBDMDB, T2D, MOMS-PI |

---

## Supersites (Major Body Regions)

| Supersite | Body Sites | Sample Count |
|-----------|------------|--------------|
| Airways | Anterior nares, Throat | ~2,000 |
| Gastrointestinal | Stool | ~3,500 |
| Oral | Tongue, Buccal, Gingiva, Palate, Saliva, Plaque | ~4,000 |
| Skin | Retroauricular crease, Antecubital fossa | ~1,000 |
| Urogenital | Vagina, Posterior fornix | ~500 |

---

## Body Site Ontology Mapping

| Body Site | UBERON ID |
|-----------|-----------|
| Stool | UBERON:0001988 |
| Tongue dorsum | UBERON:0009471 |
| Buccal mucosa | UBERON:0006956 |
| Anterior nares | UBERON:0001707 |
| Supragingival plaque | UBERON:0016485 |
| Posterior fornix | UBERON:0012247 |

---

## HMP Studies

| Study | Focus | Data Types |
|-------|-------|------------|
| HMP1 | Healthy cohort | 16S, WGS |
| IBDMDB | Inflammatory bowel disease | Multi-omics |
| T2D | Type 2 diabetes | Multi-omics |
| MOMS-PI | Pregnancy/preterm birth | Multi-omics |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| HMP | Human Microbiome Project |
| iHMP | Integrative Human Microbiome Project |
| OSDF | Open Science Data Framework |
| MIxS | Minimum Information about any Sequence |
| FMA | Foundational Model of Anatomy |
| UBERON | Uber Anatomy Ontology |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
