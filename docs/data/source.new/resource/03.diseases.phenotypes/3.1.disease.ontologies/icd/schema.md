---
id: schema-icd
title: "International Classification of Diseases (ICD) Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, ontology, classification, diseases, who, clinical-coding, icd-10, icd-11]
---

# International Classification of Diseases (ICD) Schema Documentation

**Document ID:** SCHEMA-ICD

---

## TL;DR

ICD is the WHO global standard for diagnostic health information, providing a common language for recording, reporting, and monitoring diseases. ICD-10 contains ~70,000 codes while ICD-11 contains ~80,000 entities with improved digital functionality.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| ICD-10 Codes | ~70,000 | WHO |
| ICD-11 Entities | ~80,000 | WHO |
| Chapters (ICD-10) | 22 | WHO |
| Countries Using | 194 | WHO |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | XML, JSON, Tabular |
| API | Yes (ICD API - registration required) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| ICD-10-CM | `[A-Z][0-9]{2}(\.[0-9]{1,4})?` | E11.9 (Type 2 diabetes) |
| ICD-10 | `[A-Z][0-9]{2}(\.[0-9])?` | I21.0 (MI) |
| ICD-11 | `[A-Z0-9]{4,6}` | 5A11 (Type 2 diabetes) |
| Block Code | `[A-Z][0-9]{2}-[A-Z][0-9]{2}` | E10-E14 |

---

## References

See [Overview](./_index.md) for full details.
