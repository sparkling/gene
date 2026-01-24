---
id: methodology-curation-framework
title: "Data Sources Pruning Analysis"
type: methodology
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [methodology, curation, data-quality, pruning]
---

**Parent:** [Methodology](_index.md)

# Data Sources Pruning Analysis

**Created:** January 2026
**Updated:** January 2026 (added original research evaluation)
**Purpose:** Critical evaluation of ALL databases (original ~40 + new 200+) for removal

---

## Executive Summary

| Pruning Reason | New DBs | Original DBs | Total |
|----------------|---------|--------------|-------|
| **Commercial/Paid License** | 18 | 0 | 18 |
| **Controlled Access** | 15 | 0 | 15 |
| **Legal/Scraping Risk** | 12 | 0 | 12 |
| **Out of Scope** | 22 | 2 | 24 |
| **Redundant/Superseded** | 14 | 6 | 20 |
| **Defunct/Unavailable** | 6 | 1 | 7 |
| **Low Value** | 16 | 2 | 18 |
| **TOTAL PRUNABLE** | **103** | **11** | **114** |
| **KEEP (New)** | ~97 | - | - |
| **KEEP (Original)** | - | ~29 | - |
| **TOTAL KEEP** | - | - | **~126** |

---

[Full content from refactor.plan.md continues...]

---

## Download

| Resource | Method | URL |
|----------|--------|-----|
| **Curation checklist** | This document | N/A |
| **Database inventory** | Internal | See curation spreadsheet |
| **Quality metrics** | Generated | Pipeline output |

**Access Requirements:** Internal documentation; external sources vary.

## Data Format

| Format | Description |
|--------|-------------|
| Documentation | Markdown |
| Inventory | CSV, TSV |
| Metrics | JSON |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `database_name` | string | Source database | "gnomAD" |
| `tier` | integer | Quality tier (1-3) | 1 |
| `status` | string | Curation status | "keep" |
| `reason` | string | Decision rationale | "Primary source" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `replaces` | Database | 1:N |
| `supplements` | Database | N:M |

## Sample Data

### Example Curation Decision
```json
{
  "database": "legacy_db",
  "tier": 3,
  "status": "prune",
  "reason": "redundant",
  "replacement": "gnomAD v4.1",
  "reviewed_date": "2026-01-15"
}
```

### Sample Query Result
| database | tier | status | reason |
|----------|------|--------|--------|
| gnomAD | 1 | keep | Primary |
| old_variant_db | 3 | prune | Redundant |

## License

| Resource | License | Notes |
|----------|---------|-------|
| Framework | Internal | Project documentation |
| External sources | Varies | See individual databases |

## Data Set Size

| Metric | Value |
|--------|-------|
| Databases evaluated | 240 |
| Databases retained | ~126 |
| Databases pruned | ~114 |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Data Curation | Process of organizing, validating, and maintaining data quality | Removing duplicates |
| Pruning | Removing data sources that don't meet quality or relevance criteria | Excluding defunct DBs |
| Data Source Evaluation | Assessment of database quality, licensing, and value | Tier 1/2/3 ranking |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Controlled Access | Data requiring application or approval for access | Restricted databases |
| Commercial License | Paid license required for data use | Enterprise pricing |
| Redundant Database | Database whose data is available in other sources | Superseded sources |
| Defunct Database | Database no longer maintained or available | Dead links |
| Out of Scope | Database not relevant to project requirements | Non-genomics data |
| Low Value | Database providing minimal unique information | Limited data |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CC BY | Creative Commons Attribution | Open license |
| CC BY-NC | Creative Commons Attribution Non-Commercial | Restricted license |
| DB | Database | Data source |
| MVP | Minimum Viable Product | Initial release |

---

*Generated: January 2026*
