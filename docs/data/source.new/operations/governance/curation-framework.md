---
id: governance-curation-framework
title: "Data Sources Pruning Analysis"
type: governance
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [governance, curation, data-quality, pruning]
---

**Parent:** [Governance](./_index.md)

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
