# FDA OpenFDA - Data Dictionary

## Overview

This data dictionary documents FDA OpenFDA drug adverse event records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | fda.openfda |
| **Name** | FDA OpenFDA |
| **Total Fields** | 40+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Report Identifiers

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| safetyreportid | string | Unique 8-digit report ID |
| safetyreportversion | string | Report version |
| receivedate | string | Date received (YYYYMMDD) |

### Seriousness

| Field | Value | Meaning |
|-------|-------|---------|
| serious | 1 | Serious event |
| serious | 2 | Non-serious |
| seriousnessdeath | 1 | Death occurred |
| seriousnesshospitalization | 1 | Hospitalization |

### Drug Characterization

| Value | Meaning |
|-------|---------|
| 1 | Suspect drug |
| 2 | Concomitant drug |
| 3 | Interacting drug |

### Reaction Outcome

| Value | Meaning |
|-------|---------|
| 1 | Recovered |
| 2 | Recovering |
| 3 | Not recovered |
| 4 | Recovered with sequelae |
| 5 | Fatal |
| 6 | Unknown |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| FAERS | FDA Adverse Event Reporting System |
| MedDRA | Medical Dictionary for Regulatory Activities |
| NDC | National Drug Code |
| RXCUI | RxNorm Concept Unique Identifier |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
