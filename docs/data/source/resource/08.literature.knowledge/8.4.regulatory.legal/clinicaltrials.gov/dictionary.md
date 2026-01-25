# ClinicalTrials.gov - Data Dictionary

## Overview

This data dictionary documents ClinicalTrials.gov study records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | clinicaltrials.gov |
| **Name** | ClinicalTrials.gov |
| **Total Fields** | 50+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identifiers

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| nctId | string | Yes | NCT identifier | `NCT00000001` |
| briefTitle | string | Yes | Short title | Study title |
| officialTitle | string | No | Full title | Complete title |

### Status

| Value | Description |
|-------|-------------|
| NOT_YET_RECRUITING | Approved but not started |
| RECRUITING | Actively enrolling |
| ACTIVE_NOT_RECRUITING | Ongoing, closed to enrollment |
| COMPLETED | Finished |
| TERMINATED | Ended early |
| WITHDRAWN | Never started |

### Phases

| Phase | Description |
|-------|-------------|
| EARLY_PHASE1 | Exploratory |
| PHASE1 | Safety/dosing |
| PHASE2 | Efficacy signals |
| PHASE3 | Confirmatory |
| PHASE4 | Post-market |
| NA | Not applicable |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| NCT | National Clinical Trial |
| IRB | Institutional Review Board |
| DSMB | Data Safety Monitoring Board |
| ITT | Intention to Treat |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
