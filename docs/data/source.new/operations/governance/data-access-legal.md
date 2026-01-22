---
id: governance-data-access-legal
title: "Data Access: Legal Considerations & Governance"
type: governance
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [governance, legal, compliance, data-access]
---

**Parent:** [Governance](./_index.md)

# Data Access: Legal Considerations & Governance

**Last Updated:** January 2026
**Purpose:** Actionable guidance for legitimate data access and legal risk mitigation

---

## Executive Summary

| Risk Level | Approach | Examples |
|------------|----------|----------|
| **LOW** | Direct access via APIs or partnerships | Clinical Trials.gov, AACT, PROMIS, Open Humans |
| **MEDIUM** | Restricted access with consent/licensing | Wearable APIs, Health Connect, Fitabase |
| **HIGH** | Avoid or partnership-only | Reddit, ConsumerLab, Labdoor, HealthUnlocked |

**Key Principle:** Use official APIs, respect Terms of Service, and obtain proper consent.

---

[Full content from data-access-legal.md continues...]

---

## Download

| Resource | Method | Notes |
|----------|--------|-------|
| **This guide** | Internal document | N/A |
| **License inventory** | Internal spreadsheet | Tracks all database licenses |
| **API agreements** | Per provider | Registration required |

**Access Requirements:** This is internal guidance; external API access varies by provider.

## Data Format

| Format | Description |
|--------|-------------|
| Documentation | Markdown |
| License tracking | CSV |
| Compliance records | JSON |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `source` | string | Data provider | "PubMed" |
| `license_type` | string | License category | "public_domain" |
| `risk_level` | string | Legal risk (low/medium/high) | "low" |
| `requires_consent` | boolean | User consent needed | false |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `requires` | Agreement | 1:N |
| `permits` | Use Case | 1:N |

## Sample Data

### Example License Record
```json
{
  "source": "ClinicalTrials.gov",
  "license": "public_domain",
  "risk_level": "low",
  "commercial_use": true,
  "requires_consent": false,
  "api_available": true
}
```

### Sample Query Result
| source | license | risk | commercial |
|--------|---------|------|------------|
| PubMed | Public domain | Low | Yes |
| DrugBank | Academic | Medium | No (license required) |

## License

| Resource | License | Notes |
|----------|---------|-------|
| This guide | Internal | Project documentation |
| Data sources | Varies | See individual assessments |

## Data Set Size

| Metric | Value |
|--------|-------|
| Sources evaluated | 150+ |
| Low risk | ~60% |
| Medium risk | ~30% |
| High risk | ~10% |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Terms of Service (ToS) | Legal agreement governing use of a service or data | API usage restrictions |
| Data Licensing | Legal framework defining how data can be used | CC BY, academic license |
| Risk Mitigation | Strategies to reduce legal and compliance exposure | Using official APIs |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| API Access | Programmatic data retrieval via official interface | REST, GraphQL |
| Web Scraping | Automated extraction of data from websites | Legal risk |
| Consent | Permission from users or data owners for data use | GDPR, HIPAA |
| Partnership Access | Data access through formal business agreements | Enterprise licenses |
| Public Domain | Data with no copyright restrictions | CC0 |
| Controlled Access | Data requiring application/approval for access | dbGaP |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | Official data access |
| CC BY | Creative Commons Attribution | Requires attribution |
| CC BY-NC | Creative Commons Attribution Non-Commercial | No commercial use |
| CC BY-SA | Creative Commons Attribution Share-Alike | Derivative works |
| CC0 | Creative Commons Zero (Public Domain) | No restrictions |
| GDPR | General Data Protection Regulation | EU privacy law |
| HIPAA | Health Insurance Portability and Accountability Act | US health data law |
| ToS | Terms of Service | Usage agreement |

---

**Next Steps:** Reference this document when planning data integration. For specific sources, confirm with legal team before implementation.
